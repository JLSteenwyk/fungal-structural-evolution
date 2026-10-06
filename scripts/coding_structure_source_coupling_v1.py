"""Link inherited coding evidence and representative model availability without admission."""
from collections import Counter, defaultdict
import csv
import gzip
import json
from pathlib import Path
import sqlite3


EXACT = 'original_target_genome_agrees_and_inherited_strict_translation_exact'
REVIEW = 'original_target_review_or_no_strict_translation_agreement'
DERIVED = 'no_original_target_with_separate_derived_translation_evidence'
MISSING = 'no_original_target_or_exact_derived_translation_evidence'


def table(path):
    with Path(path).open() as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def json_rows(path):
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            yield json.loads(line)


def evidence(entry):
    kind = entry['mapping_mode']
    primary, supplements = {}, {}
    if kind == 'ncbi_protein_gff':
        for ordinal, row in enumerate(table(entry['translation_table']), 1):
            assert row['taxon_id'] == entry['taxon_id']
            if row['status'] == 'protein_without_cds_record':
                continue
            assert row['cds_id'] not in primary
            primary[row['cds_id']] = dict(inherited_status=row['status'], protein_id=row['protein_id'] or None,
                translation_table=row['translation_table'] or None, code_provenance=row['code_source'] or None,
                terminal_stop=row['terminal_stop'] or None, annotation_flags=row['gff_flags'],
                original_header=row['cds_header'], source_kind='ncbi_unmodified_publisher_cds',
                evidence_path=entry['translation_table'], evidence_ordinal=ordinal,
                translation_independently_recomputed_here=False)
    elif kind == 'transcript_gff':
        for ordinal, row in enumerate(table(entry['translation_table']), 1):
            if row['taxon_id'] != entry['taxon_id']:
                continue
            pid = row['protein_id']
            assert pid not in primary
            primary[pid] = dict(inherited_status=row['status'], protein_id=pid, translation_table='1',
                code_provenance='table_1_test_in_original_publisher_audit', terminal_stop=row['terminal_stop_in_cds'],
                annotation_flags=None, source_kind='external_unmodified_publisher_cds',
                evidence_path=entry['translation_table'], evidence_ordinal=ordinal,
                original_dna_sha256=row['cds_sequence_sha256'], original_protein_sha256=row['protein_sequence_sha256'],
                translation_independently_recomputed_here=False)
        for ordinal, row in enumerate(table(entry['supplementary_table']), 1):
            if row['taxon_id'] == entry['taxon_id']:
                assert row['protein_id'] not in supplements
                supplements[row['protein_id']] = dict(kind='separate_annotation_boundary_codon_projection',
                    evidence_path=entry['supplementary_table'], evidence_ordinal=ordinal,
                    original_row=row, shared_original_genome_dependency=True,
                    exact_derived_translation=row['projection_status'] == 'exact_genome_linked_codon_translation')
    elif kind == 'verified_orf_coordinates':
        for ordinal, row in enumerate(table(entry['translation_table']), 1):
            assert row['protein_id'] not in primary
            primary[row['protein_id']] = dict(inherited_status=row['status'], protein_id=row['protein_id'],
                translation_table='1', code_provenance='table_1_test_in_original_dependent_orf_audit',
                terminal_stop=None, annotation_flags=None, source_kind='dependent_genome_derived_orf_cds',
                evidence_path=entry['translation_table'], evidence_ordinal=ordinal,
                shared_original_genome_dependency=True, translation_independently_recomputed_here=False)
    else:
        assert kind == 'creolimax_gtf'
        for ordinal, row in enumerate(table(entry['supplementary_table']), 1):
            assert row['protein_id'] not in supplements
            supplements[row['protein_id']] = dict(kind='separate_strict_annotation_derived_cds',
                evidence_path=entry['supplementary_table'], evidence_ordinal=ordinal,
                original_row=row, shared_original_genome_dependency=True,
                exact_derived_translation=row['status'] == 'exact_translation')
    return primary, supplements


def target_record(target, inherited):
    assert target['product_id'] == inherited['protein_id']
    if 'original_header' in inherited:
        assert target['description'] == inherited['original_header']
    if 'original_dna_sha256' in inherited:
        assert target['sequence_sha256'] == inherited['original_dna_sha256']
    exact = target['status'] == 'one_exact_unmodified_genomic_cds_candidate' and inherited['inherited_status'] == 'exact_translation'
    return dict(ordinal=target['ordinal'], cds_id=target['cds_id'], protein_id=target['product_id'],
        original_target_dna_sha256=target['sequence_sha256'], original_target_length=target['sequence_length'],
        genomic_status=target['status'], genomic_candidate_count=target['candidate_count'],
        matching_candidate_indices=target['matching_candidate_indices'],
        translation_evidence=inherited, joined_status=EXACT if exact else REVIEW,
        biological_codon_eligibility=False)


def product_record(product, targets, supplement, availability):
    selected = product['selected_representative']
    if selected:
        assert availability is not None and availability['protein_id'] == product['protein_id']
        assert availability['sequence_sha256'] == product['sequence_sha256']
        assert availability['length'] == product['protein_length']
    else:
        assert availability is None
    expected = [dict(cds_id=t['cds_id'], ordinal=t['ordinal'], status=t['genomic_status']) for t in targets]
    assert expected == product['original_cds_targets']
    if len(targets) == 1 and targets[0]['joined_status'] == EXACT:
        status = EXACT
    elif targets:
        status = REVIEW
    elif supplement is not None and supplement['exact_derived_translation']:
        status = DERIVED
    else:
        status = MISSING
    return dict(protein_id=product['protein_id'], protein_length=product['protein_length'],
        sequence_sha256=product['sequence_sha256'], selected_representative=selected,
        original_mapping=product['original_mapping'], original_decision=product['original_decision'],
        original_genomic_candidate_count=product['genomic_candidate_count'],
        original_missing_target_reason=product['missing_target_reason'],
        original_target_evidence=targets, separate_derived_evidence=supplement,
        representative_availability=availability,
        availability_scope='selected_representative_baseline' if selected else 'alternative_product_not_assigned_by_this_baseline',
        joined_status=status, biological_codon_eligibility=False, scientific_eligibility=False)


def availability_rows(database, taxon):
    con = sqlite3.connect(Path(database).resolve().as_uri() + '?mode=ro', uri=True)
    result = {}
    for pid, digest, length, afdb, esmfold in con.execute(
            'SELECT protein_id,sequence_sha256,length,has_afdb,has_esmfold FROM proteins WHERE taxon_id=?', (taxon,)):
        assert pid not in result
        state = 'both' if afdb and esmfold else 'afdb_only' if afdb else 'esmfold_only' if esmfold else 'neither'
        result[pid] = dict(protein_id=pid, sequence_sha256=digest, length=length, availability=state)
    con.close()
    return result


def couple_taxon(entry, database, root):
    taxon = entry['taxon_id']
    root = Path(root) / taxon
    root.mkdir(exist_ok=False)
    inherited, supplemental = evidence(entry)
    model_links = availability_rows(database, taxon)
    targets_by_product = defaultdict(list)
    target_counts, product_counts, cross_counts = Counter(), Counter(), Counter()
    target_path = root / 'target_evidence.jsonl.gz'
    target_count = 0
    with gzip.open(target_path, 'wt') as handle:
        for target in json_rows(entry['genomic_targets']):
            assert target['cds_id'] in inherited
            row = target_record(target, inherited[target['cds_id']])
            handle.write(json.dumps(row, separators=(',', ':')) + '\n')
            target_counts[row['joined_status']] += 1
            target_count += 1
            if row['protein_id'] is not None:
                targets_by_product[row['protein_id']].append(row)
    assert target_count == entry['target_records']
    path = root / 'product_evidence.jsonl.gz'
    products = selected = 0
    used_links, used_supplements = set(), set()
    with gzip.open(path, 'wt') as handle:
        for product in json_rows(entry['genomic_products']):
            pid = product['protein_id']
            link = model_links.get(pid) if product['selected_representative'] else None
            supplement = supplemental.get(pid)
            row = product_record(product, targets_by_product.get(pid, []), supplement, link)
            for target in row['original_target_evidence']:
                proof = target['translation_evidence']
                if 'original_protein_sha256' in proof:
                    assert proof['original_protein_sha256'] == product['sequence_sha256']
            handle.write(json.dumps(row, separators=(',', ':')) + '\n')
            products += 1
            selected += int(row['selected_representative'])
            product_counts[row['joined_status']] += 1
            cross_counts[row['joined_status'] + ':' + (link['availability'] if link else 'alternative_not_assigned')] += 1
            if link:
                used_links.add(pid)
            if supplement:
                used_supplements.add(pid)
    assert products == entry['source_products'] and selected == entry['selected_representatives']
    assert used_links == set(model_links) and used_supplements == set(supplemental)
    return dict(taxon_id=taxon, mapping_mode=entry['mapping_mode'], study_role=entry['study_role'],
        manifest_species_name=entry['species_name'], manifest_lineage=entry['lineage'],
        source_products=products, selected_representatives=selected, target_records=target_count,
        target_status_counts=dict(target_counts), product_status_counts=dict(product_counts),
        coding_structure_cross_counts=dict(cross_counts), derived_evidence_products=len(used_supplements),
        artifact_paths=[str(target_path), str(path)])
