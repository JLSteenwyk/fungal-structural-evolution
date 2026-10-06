#!/usr/bin/env python3
"""Synthetic whole-source integration for the independent code-context reader."""
import argparse
from datetime import datetime, timezone
import gzip
import io
import json
from pathlib import Path
import tarfile
import csv

from genetic_code_context_v1 import audit_taxon, ancestry, inherited_code, read_taxonomy
from readback_full_genetic_code_context_v1 import source_contexts, check_taxon
from reference_measurement_union_sources import bind


def save_rows(path, rows):
    with gzip.open(path, 'wt') as f:
        for row in rows:
            f.write(json.dumps(row) + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    args = p.parse_args()
    assert not args.receipt.exists()
    args.output.mkdir(exist_ok=False)
    archive = args.output / 'literal_taxdump.tar.gz'
    raw_nodes = []
    for tid, parent, n, nf, m, mf in [(1, 1, 1, 0, 0, 0), (2, 1, 12, 0, 4, 0), (3, 2, 12, 1, 4, 1)]:
        values = [str(tid), str(parent), 'species' if tid != 1 else 'no rank', '', '8', '0',
                  str(n), str(nf), str(m), str(mf), '0', '0', '', '', '', '0', '0', '0']
        raw_nodes.append('\t|\t'.join(values) + '\t|\n')
    data = {'nodes.dmp': ''.join(raw_nodes), 'merged.dmp': '20\t|\t3\t|\n',
            'gencode.dmp': ''.join('\t|\t'.join([str(c), '', 'Literal ' + str(c), 'A' * 64, '-' * 64]) + '\t|\n'
                                   for c in [0, 1, 4, 12])}
    with tarfile.open(archive, 'w:gz') as tar:
        for name, body in data.items():
            blob = body.encode()
            info = tarfile.TarInfo(name)
            info.size = len(blob)
            tar.addfile(info, io.BytesIO(blob))
    identity = args.output / 'literal_identity.tsv'
    label = dict(taxon_id='literal', study_role='ingroup', effective_species_taxid='20',
                 canonical_effective_taxid='3', node_path_taxids='3;2;1', catalogue_taxid='3',
                 frozen_species_name='Literal species', snapshot_scientific_name='Literal species',
                 identity_review_flags='literal_review', biological_species_delimitation='unverified')
    with identity.open('x', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(label), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerow(label)
    independent, codebook = source_contexts(archive, identity)
    nodes, merged, codes = read_taxonomy(archive)
    assert merged[20] == 3
    context = independent['literal']
    assert context['source_parent_path'] == ancestry(3, nodes)
    assert context['nuclear'] == inherited_code(3, nodes, 'nuclear')
    assert context['mitochondrial'] == inherited_code(3, nodes, 'mitochondrial')
    assert codebook == {str(k): v for k, v in codes.items()}
    target = dict(ordinal=1, cds_id='one', protein_id='protein', original_target_dna_sha256='literal-dna',
                  joined_status='original_exact', translation_evidence=dict(inherited_status='exact_translation',
                  translation_table='1', code_provenance='table_1_assumption_no_explicit_code'))
    unlinked = dict(target, ordinal=2, cds_id='unlinked', protein_id=None,
                    translation_evidence=dict(target['translation_evidence'], translation_table=None))
    product = dict(protein_id='protein', sequence_sha256='literal-aa', selected_representative=True,
                   joined_status='original_exact', original_target_evidence=[target],
                   separate_derived_evidence=None, representative_availability=dict(availability='both'))
    alternate = dict(product, protein_id='alternative', selected_representative=False,
                     original_target_evidence=[], representative_availability=None,
                     separate_derived_evidence=dict(kind='separate_strict_annotation_derived_cds',
                     exact_derived_translation=True, shared_original_genome_dependency=True))
    source_targets, source_products = args.output / 'source_targets.jsonl.gz', args.output / 'source_products.jsonl.gz'
    save_rows(source_targets, [target, unlinked])
    save_rows(source_products, [product, alternate])
    entry = dict(taxon_id='literal', targets=str(source_targets), products=str(source_products),
                 target_records=2, source_products=2, selected_representatives=1, context=context)
    (args.output / 'joined').mkdir(exist_ok=False)
    report = audit_taxon(entry, args.output / 'joined')
    proof = check_taxon(entry, report, context)
    assert proof['target_records'] == proof['source_products'] == 2
    assert proof['selected_representatives'] == proof['separate_derived_evidence_products'] == 1
    with gzip.open(report['artifact_paths'][0], 'rt') as f:
        targets = [json.loads(line) for line in f]
    with gzip.open(report['artifact_paths'][1], 'rt') as f:
        products = [json.loads(line) for line in f]
    cases = []
    mutations = [
        ('reject_changed_code_relationship', 'product', 'snapshot_code_classification', 'invented'),
        ('reject_false_code_admission', 'product', 'genetic_code_admission', True),
        ('reject_changed_original_code_provenance', 'target', 'inherited_code_provenance', 'invented')]
    for name, kind, field, value in mutations:
        original = products if kind == 'product' else targets
        altered = json.loads(json.dumps(original))
        altered[0][field] = value
        path = args.output / (name + '.jsonl.gz')
        save_rows(path, altered)
        modified = dict(report, artifact_paths=list(report['artifact_paths']))
        modified['artifact_paths'][1 if kind == 'product' else 0] = str(path)
        try:
            check_taxon(entry, modified, context)
        except AssertionError:
            cases.append(name)
        else:
            raise AssertionError(name)
    pins = {}
    for path in [Path(__file__), Path('scripts/genetic_code_context_v1.py'),
                 Path('scripts/readback_full_genetic_code_context_v1.py'),
                 *[q for q in args.output.rglob('*') if q.is_file()]]:
        bind(pins, path)
    result = dict(status='passed_independent_genetic_code_context_reader_literal_integration',
        checked_utc=datetime.now(timezone.utc).isoformat(), literal_targets=2, literal_products=2,
        unlinked_targets_retained=1, alternatives_retained=1, separate_derived_products_retained=1,
        independent_taxonomy_inheritance=True, rejected_semantic_mutations=cases,
        source_hashes=pins, scientific_eligibility=False, genetic_code_admission=False, gpu=False, new_predictions=0,
        scope='Literal source archive/identity/target/product integration and three semantic mutations '
              'invoked directly without a hash gate or hash-rebinding exercise. No real taxon pilot, '
              'current taxonomy, code/compartment admission, translation or biological inference.')
    with args.receipt.open('x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
