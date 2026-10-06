"""Compare inherited CDS test codes to pinned taxonomic context without recoding."""
from collections import Counter
import csv
import gzip
import json
from pathlib import Path
import tarfile


def fields(line):
    pieces = line.decode('utf-8').rstrip('\n').split('|')
    assert pieces.pop().strip() == ''
    return [p.strip() for p in pieces]


def read_taxonomy(archive):
    nodes, merged, codes = {}, {}, {}
    with tarfile.open(archive) as tar:
        for line in tar.extractfile('nodes.dmp'):
            row = fields(line)
            assert len(row) >= 13
            taxid, parent = int(row[0]), int(row[1])
            assert taxid not in nodes and row[7] in ('0', '1') and row[9] in ('0', '1')
            nodes[taxid] = (parent, int(row[6]), int(row[7]), int(row[8]), int(row[9]))
        for line in tar.extractfile('merged.dmp'):
            old, new = map(int, fields(line))
            assert old not in merged
            merged[old] = new
        for line in tar.extractfile('gencode.dmp'):
            row = fields(line)
            assert len(row) == 5
            code = int(row[0])
            assert code not in codes
            codes[code] = dict(abbreviation=row[1], name=row[2], amino_acids=row[3], starts=row[4])
    return nodes, merged, codes


def canonical(taxid, merged):
    taxid, seen = int(taxid), set()
    while taxid in merged:
        assert taxid not in seen, 'merged taxonomy cycle'
        seen.add(taxid)
        taxid = merged[taxid]
    return taxid


def ancestry(taxid, nodes):
    seen, path = set(), []
    while True:
        assert taxid not in seen and taxid in nodes, 'taxonomy cycle or missing parent'
        seen.add(taxid)
        path.append(taxid)
        parent = nodes[taxid][0]
        if parent == taxid:
            assert taxid == 1
            return path
        taxid = parent


def inherited_code(taxid, nodes, compartment):
    index = 1 if compartment == 'nuclear' else 3
    assert compartment in ('nuclear', 'mitochondrial')
    path, seen = [], set()
    while True:
        assert taxid in nodes and taxid not in seen, 'code inheritance cycle or missing node'
        seen.add(taxid)
        node = nodes[taxid]
        path.append(dict(taxid=taxid, recorded_code=node[index], inherited=bool(node[index + 1])))
        if not node[index + 1]:
            code = node[index]
            return dict(code=code, source_taxid=taxid, inheritance_path=path,
                        recorded_values_agree=all(r['recorded_code'] == code for r in path))
        taxid = node[0]


def relationship(code, context):
    if code in (None, ''):
        return 'inherited_test_code_not_recorded'
    code = int(code)
    if code == 0:
        return 'inherited_test_code_unspecified'
    nuclear, mitochondrial = context['nuclear']['code'], context['mitochondrial']['code']
    if code == nuclear == mitochondrial:
        return 'matches_snapshot_nuclear_and_mitochondrial_codes_compartment_unverified'
    if code == nuclear:
        return 'matches_snapshot_nuclear_code_compartment_unverified'
    if code == mitochondrial:
        return 'matches_snapshot_mitochondrial_code_requires_compartment_verification'
    return 'differs_from_snapshot_nuclear_and_mitochondrial_codes_requires_review'


def target_context(row, context):
    evidence = row['translation_evidence']
    return dict(ordinal=row['ordinal'], cds_id=row['cds_id'], protein_id=row['protein_id'],
        original_target_dna_sha256=row['original_target_dna_sha256'], original_joined_status=row['joined_status'],
        inherited_translation_status=evidence['inherited_status'], inherited_test_code=evidence['translation_table'],
        inherited_code_provenance=evidence['code_provenance'],
        snapshot_code_relationship=relationship(evidence['translation_table'], context),
        taxon_context_id=context['taxon_id'], genetic_code_admission=False)


def product_context(row, context):
    originals = [target_context(t, context) for t in row['original_target_evidence']]
    supplement = row['separate_derived_evidence']
    derived = None
    if supplement is not None:
        assert supplement['kind'] in ('separate_annotation_boundary_codon_projection', 'separate_strict_annotation_derived_cds')
        derived = dict(kind=supplement['kind'], test_code=1,
                       source='table_1_test_in_original_separate_derived_audit',
                       exact_derived_translation=supplement['exact_derived_translation'],
                       snapshot_code_relationship=relationship(1, context),
                       shared_original_genome_dependency=supplement['shared_original_genome_dependency'])
    relations = {r['snapshot_code_relationship'] for r in originals}
    if len(relations) > 1:
        classification = 'multiple_original_code_context_relationships_require_review'
    elif relations:
        classification = next(iter(relations))
    elif derived is not None:
        classification = 'separate_derived_test:' + derived['snapshot_code_relationship']
    else:
        classification = 'no_original_or_separate_derived_test_code'
    model = row['representative_availability']
    return dict(protein_id=row['protein_id'], sequence_sha256=row['sequence_sha256'],
        selected_representative=row['selected_representative'], original_joined_status=row['joined_status'],
        original_target_code_contexts=originals, separate_derived_code_context=derived,
        availability=model['availability'] if model is not None else 'alternative_not_assigned',
        snapshot_code_classification=classification, taxon_context_id=context['taxon_id'],
        original_target_count=len(originals), genetic_code_admission=False, biological_codon_eligibility=False)


def rows(path):
    with gzip.open(path, 'rt') as f:
        for line in f:
            yield json.loads(line)


def audit_taxon(entry, output):
    root = Path(output) / entry['taxon_id']
    root.mkdir(exist_ok=False)
    targets, products, cross = Counter(), Counter(), Counter()
    target_path, product_path = root / 'target_code_context.jsonl.gz', root / 'product_code_context.jsonl.gz'
    with gzip.open(target_path, 'wt') as f:
        for source in rows(entry['targets']):
            row = target_context(source, entry['context'])
            f.write(json.dumps(row, separators=(',', ':')) + '\n')
            targets[row['snapshot_code_relationship']] += 1
    selected, derived = 0, 0
    with gzip.open(product_path, 'wt') as f:
        for source in rows(entry['products']):
            row = product_context(source, entry['context'])
            f.write(json.dumps(row, separators=(',', ':')) + '\n')
            classification = row['snapshot_code_classification']
            products[classification] += 1
            cross[row['original_joined_status'] + ':' + row['availability'] + ':' + classification] += 1
            selected += int(row['selected_representative'])
            derived += int(row['separate_derived_code_context'] is not None)
    assert sum(targets.values()) == entry['target_records']
    assert sum(products.values()) == entry['source_products']
    assert selected == entry['selected_representatives']
    return dict(taxon_id=entry['taxon_id'], target_records=sum(targets.values()), source_products=sum(products.values()),
                selected_representatives=selected, separate_derived_evidence_products=derived,
                target_relationship_counts=dict(targets), product_classification_counts=dict(products),
                coding_model_code_cross_counts=dict(cross), artifact_paths=[str(target_path), str(product_path)])
