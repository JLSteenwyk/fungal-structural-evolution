"""Full source-target translation diagnostics without sequence or code correction."""
from collections import Counter, defaultdict
import gzip
from itertools import zip_longest
import json
from pathlib import Path

from Bio.SeqIO.FastaIO import SimpleFastaParser
from unmodified_cds_translation_v1 import diagnose, digest


def fasta(path):
    path = Path(path)
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as handle:
        for header, sequence in SimpleFastaParser(handle):
            yield header.split()[0], sequence


def json_rows(path):
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            yield json.loads(line)


def emit(handle, row):
    handle.write(json.dumps(row, separators=(',', ':'), allow_nan=False) + '\n')


def audit_taxon(entry, root):
    directory = Path(root) / entry['taxon_id']
    directory.mkdir(exist_ok=False)
    proteins = {}
    for pid, sequence in fasta(entry['proteome']):
        assert pid not in proteins
        proteins[pid] = sequence
    assert len(proteins) == entry['source_products']
    paths = [directory / name for name in ['target_diagnostics.jsonl.gz',
             'changed_original_codons.jsonl.gz', 'product_diagnostics.jsonl.gz']]
    linked = defaultdict(list)
    role_counts, inherited_counts, sites_count = Counter(), Counter(), Counter()
    targets, dna_bases = 0, 0
    originals = fasta(entry['original_cds']) if entry['original_cds'] else iter(())
    with gzip.open(paths[0], 'wt') as out, gzip.open(paths[1], 'wt') as sites:
        for ordinal, pair in enumerate(zip_longest(originals, json_rows(entry['targets'])), 1):
            source, target = pair
            assert source is not None and target is not None, (entry['taxon_id'], ordinal, 'source_count')
            cid, dna = source
            assert target['ordinal'] == ordinal and target['cds_id'] == cid
            assert len(dna) == target['original_target_length']
            pid = target['protein_id']
            assert pid is None or pid in proteins
            protein = proteins[pid] if pid is not None else None
            row, differences = diagnose(target, dna, protein, entry['context'])
            row['original_source_target'] = target
            emit(out, row)
            for change in differences:
                emit(sites, dict(taxon_id=entry['taxon_id'], ordinal=ordinal, cds_id=cid,
                    protein_id=pid, original_target_dna_sha256=target['original_target_dna_sha256'], **change))
                sites_count[change['role']] += 1
            for role, code in row['predeclared_code_roles'].items():
                state = 'unspecified_code_not_tested' if code is None else row['translation_diagnostics'][str(code)]['status']
                role_counts[role + ':' + state] += 1
            old = target['translation_evidence']['inherited_status']
            old_code = row['predeclared_code_roles']['inherited']
            recomputed = 'unspecified_code_not_tested' if old_code is None else row['translation_diagnostics'][str(old_code)]['status']
            inherited_counts[str(old) + ':' + recomputed] += 1
            summary = {k: row[k] for k in ['ordinal', 'cds_id', 'protein_id', 'original_target_dna_sha256',
                'original_dna_length', 'linked_protein_sha256', 'linked_protein_length', 'original_joined_status',
                'predeclared_code_roles', 'translation_diagnostics', 'changed_codon_counts']}
            if pid is not None:
                linked[pid].append(summary)
            targets += 1
            dna_bases += len(dna)
    assert targets == entry['target_records']
    products, representatives = 0, 0
    cross, classifications = Counter(), Counter()
    with gzip.open(paths[2], 'wt') as out:
        for pair in zip_longest(proteins.items(), json_rows(entry['products'])):
            source, original = pair
            assert source is not None and original is not None
            pid, protein = source
            assert original['protein_id'] == pid
            assert original['protein_length'] == len(protein) and original['sequence_sha256'] == digest(protein)
            observations = linked.get(pid, [])
            assert {t['ordinal'] for t in original['original_target_evidence']} == {t['ordinal'] for t in observations}
            row = dict(original_source_product=original, original_target_diagnostics=observations,
                source_target_replaced=False, derived_target_substituted=False, best_code_selected=False,
                protein_modified=False, scientific_eligibility=False, biological_codon_eligibility=False)
            emit(out, row)
            availability = original['representative_availability']
            avail = availability['availability'] if availability is not None else 'alternative_not_assigned'
            for role in ['inherited', 'snapshot_nuclear', 'snapshot_mitochondrial']:
                statuses = set()
                for observation in observations:
                    code = observation['predeclared_code_roles'][role]
                    statuses.add('unspecified_code_not_tested' if code is None else
                                 observation['translation_diagnostics'][str(code)]['status'])
                classification = '|'.join(sorted(statuses)) if statuses else 'no_original_target'
                classifications[role + ':' + classification] += 1
                cross[original['joined_status'] + ':' + avail + ':' + role + ':' + classification] += 1
            products += 1
            representatives += int(original['selected_representative'])
    assert products == entry['source_products'] and representatives == entry['selected_representatives']
    return dict(taxon_id=entry['taxon_id'], source_products=products, selected_representatives=representatives,
        target_records=targets, original_dna_bases=dna_bases, changed_codon_rows=sum(sites_count.values()),
        target_role_status_counts=dict(role_counts), inherited_status_cross_counts=dict(inherited_counts),
        changed_codon_role_counts=dict(sites_count), product_role_status_counts=dict(classifications),
        coding_model_translation_cross_counts=dict(cross), artifact_paths=[str(p) for p in paths])
