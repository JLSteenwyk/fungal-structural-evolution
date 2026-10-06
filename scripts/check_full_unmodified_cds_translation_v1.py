#!/usr/bin/env python3
"""Literal full writer/independent reader integration and semantic corruption checks."""
import argparse
from copy import deepcopy
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path

from full_unmodified_cds_translation_v1 import audit_taxon
from readback_full_unmodified_cds_translation_v1 import check_taxon
from reference_measurement_union_sources import bind


def sha_sequence(sequence):
    return hashlib.sha256(sequence.encode('ascii')).hexdigest()


def write_rows(path, rows):
    with gzip.open(path, 'wt') as handle:
        for row in rows:
            handle.write(json.dumps(row, separators=(',', ':')) + '\n')


def read_rows(path):
    with gzip.open(path, 'rt') as handle:
        return [json.loads(line) for line in handle]


def fixture(root, name, pairs, target_inputs, nuclear, mito, no_original=False):
    source = root / name
    source.mkdir()
    proteome = source / 'protein.faa'
    proteome.write_text(''.join('>' + p + '\n' + s + '\n' for p, s in pairs))
    dna_path = source / 'cds.fna'
    dna_path.write_text(''.join('>' + cid + '\n' + dna + '\n' for cid, pid, dna, code in target_inputs))
    targets = []
    for n, (cid, pid, dna, code) in enumerate(target_inputs, 1):
        targets.append(dict(ordinal=n, cds_id=cid, protein_id=pid,
            original_target_dna_sha256=sha_sequence(dna), original_target_length=len(dna),
            genomic_status='literal_source_review', genomic_candidate_count=1, matching_candidate_indices=[0],
            translation_evidence=dict(inherited_status='literal_saved_test', translation_table=code,
                code_provenance='predeclared_literal_source_test', annotation_flags='retain_literal_flag'),
            joined_status='literal_original_source_disposition', biological_codon_eligibility=False))
    products = []
    for pid, protein in pairs:
        selected = pid != 'p_initiator'
        products.append(dict(protein_id=pid, protein_length=len(protein), sequence_sha256=sha_sequence(protein),
            selected_representative=selected, original_mapping=dict(literal='retain'),
            original_decision=dict(literal='retain_alternative' if not selected else 'retain_representative'),
            original_genomic_candidate_count=1, original_missing_target_reason='literal_no_target' if no_original else None,
            original_target_evidence=[t for t in targets if t['protein_id'] == pid],
            separate_derived_evidence=dict(kind='literal_separate_derived_target', exact_derived_translation=True)
                if pid in ['p_no_original', 'p_derived_only'] else None,
            representative_availability=dict(availability='afdb_only') if selected else None,
            availability_scope='literal_selected' if selected else 'literal_alternative',
            joined_status='literal_original_source_disposition', biological_codon_eligibility=False,
            scientific_eligibility=False))
    target_path, product_path = source / 'targets.jsonl.gz', source / 'products.jsonl.gz'
    write_rows(target_path, targets)
    write_rows(product_path, products)
    return dict(taxon_id=name, proteome=str(proteome), original_cds=None if no_original else str(dna_path),
        targets=str(target_path), products=str(product_path), source_products=len(pairs),
        selected_representatives=sum(p['selected_representative'] for p in products), target_records=len(targets),
        context=dict(nuclear=dict(code=nuclear), mitochondrial=dict(code=mito)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['output', 'receipt']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists() and not args.receipt.exists()
    args.output.mkdir(parents=True)
    code_path = Path('metadata/full_taxonomic_genetic_code_context_20261006_v1.json')
    codebook = json.loads(code_path.read_text())['codebook']
    main_entry = fixture(args.output, 'literal12',
        [('p_exact', 'ML'), ('p_initiator', 'M'), ('p_partial', 'M'), ('p_ambiguity', 'B'),
         ('p_double', 'M'), ('p_no_original', 'W'), ('p_multi', 'M'), ('p_no_code', 'M'), ('p_lower', 'ML')],
        [('exact', 'p_exact', 'ATGCTGTAA', '1'), ('initiator', 'p_initiator', 'GTGTAA', '1'),
         ('partial', 'p_partial', 'ATGT', '1'), ('ambiguity', 'p_ambiguity', 'RAYTAA', '1'),
         ('double', 'p_double', 'ATGTAATAA', '1'), ('multi1', 'p_multi', 'ATGTAA', '1'),
         ('multi2', 'p_multi', 'ATGTAG', '1'), ('unlinked', None, 'ATGNNNTAA', None),
         ('no_code', 'p_no_code', 'ATGTAA', None), ('lowercase', 'p_lower', 'atgctgtaa', '1')], 12, 4)
    entries = [main_entry,
        fixture(args.output, 'literal6', [('p_terminal', 'M')], [('terminal', 'p_terminal', 'ATGTAA', '1')], 6, 0),
        fixture(args.output, 'literal_no_cds', [('p_derived_only', 'M')], [], 1, 1, no_original=True)]
    output = args.output / 'producer'
    output.mkdir()
    reports = [audit_taxon(e, output) for e in entries]
    checked = [check_taxon(e, output, codebook, r) for e, r in zip(entries, reports)]
    assert checked == reports
    assert sum(r['source_products'] for r in reports) == 11
    assert sum(r['target_records'] for r in reports) == 11
    assert sum(r['selected_representatives'] for r in reports) == 10
    assert sum(r['changed_codon_rows'] for r in reports) == 3
    originals = {Path(p): Path(p).read_bytes() for p in reports[0]['artifact_paths']}
    negative = []

    def rejected(name, mutate, report=None):
        for path, contents in originals.items():
            path.write_bytes(contents)
        mutate()
        try:
            check_taxon(main_entry, output, codebook, reports[0] if report is None else report)
        except AssertionError:
            negative.append(name)
        else:
            raise AssertionError('semantic corruption accepted: ' + name)

    target_path, site_path, product_path = map(Path, reports[0]['artifact_paths'])

    def alter(path, index, change):
        rows = read_rows(path)
        change(rows[index])
        write_rows(path, rows)

    rejected('reject_best_code_selection', lambda: alter(target_path, 0, lambda r: r.update(best_code_selected=True)))
    rejected('reject_code_relabeling', lambda: alter(target_path, 0,
        lambda r: r['predeclared_code_roles'].update(snapshot_nuclear=1)))
    rejected('reject_wrong_translation_digest', lambda: alter(target_path, 0,
        lambda r: r['translation_diagnostics']['12'].update(translated_sha256='changed')))
    rejected('reject_original_dna_digest_replacement', lambda: alter(target_path, 0,
        lambda r: r.update(original_target_dna_sha256='changed')))
    rejected('reject_nontriplet_translation', lambda: alter(target_path, 2,
        lambda r: r['translation_diagnostics']['1'].update(status='translated_exact_match_to_normalized_protein')))
    rejected('reject_changed_codon_position', lambda: alter(site_path, 0,
        lambda r: r.update(original_codon_index_1based=3)))
    rejected('reject_changed_codon_biological_admission', lambda: alter(site_path, 0,
        lambda r: r.update(biological_mapping_admitted=True)))
    rejected('reject_removed_alternative_product', lambda: write_rows(product_path,
        [r for r in read_rows(product_path) if r['original_source_product']['protein_id'] != 'p_initiator']))
    rejected('reject_separate_derived_target_substitution', lambda: alter(product_path, 5,
        lambda r: r.update(derived_target_substituted=True)))
    wrong_report = deepcopy(reports[0])
    wrong_report['changed_codon_rows'] += 1
    rejected('reject_taxon_aggregate_disagreement', lambda: None, wrong_report)
    for path, contents in originals.items():
        path.write_bytes(contents)
    assert check_taxon(main_entry, output, codebook, reports[0]) == reports[0]
    pins = {}
    for path in [Path(__file__), code_path, Path('scripts/full_unmodified_cds_translation_v1.py'),
                 Path('scripts/build_full_unmodified_cds_translation_v1.py'),
                 Path('scripts/readback_full_unmodified_cds_translation_v1.py'),
                 Path('scripts/unmodified_cds_translation_v1.py'), Path('scripts/independent_codon_translation_v1.py')]:
        bind(pins, path)
    for path in args.output.rglob('*'):
        if path.is_file():
            bind(pins, path)
    result = dict(status='passed_literal_unmodified_cds_writer_independent_reader_integration',
        checked_utc=datetime.now(timezone.utc).isoformat(), literal_taxa=len(entries), literal_source_products=11,
        literal_target_records=11, literal_representatives=10, literal_changed_codon_rows=3,
        semantic_rejection_controls=negative, semantic_rejection_count=len(negative), source_hashes=pins,
        scientific_eligibility=False, biological_codon_eligibility=False, gpu=False, new_predictions=0,
        scope='Synthetic full writer/independent reader integration covers original target/protein identity, '
              'ambiguity, lowercase, one-terminal-stop-only, no initiation/frame repair, alternatives, multiple '
              'targets, unlinked/no-code targets, no original CDS/separate derived evidence and terminal code '
              'differences outside protein bounds. Ten rehashed semantic corruptions rejected. No real taxon pilot.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
