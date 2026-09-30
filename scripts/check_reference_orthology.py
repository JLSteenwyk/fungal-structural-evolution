#!/usr/bin/env python3
"""Synthetic full-adapter checks for reference scope, missingness and native membership."""
import copy
import csv
import json
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path

from run_ortholog_pair_guide_comparison import sha


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + '\n')
    return str(path)


def fixture(root):
    species = root / 'SpeciesIDs.txt'
    species.write_text('0: F1.faa\n1: F2.faa\n2: F3.faa\n')
    sequences = root / 'SequenceIDs.txt'
    sequences.write_text('0_7: a\n1_0: c\n0_2: b\n2_17: d\n2_18: e\n')
    bridge = root / 'bridge.sqlite'
    with sqlite3.connect(bridge) as conn:
        conn.execute('CREATE TABLE structures (taxon_id TEXT, protein_id TEXT, model_id TEXT, version INTEGER)')
        conn.executemany('INSERT INTO structures VALUES (?,?,?,?)', [('F3', 'd', 'M-D', 10), ('F3', 'e', 'M-E', 6)])
    source = root / 'source'
    source.mkdir()
    specs = [(['F3_d'], ['F2_c'], 'provisional_reference_available'),
             (['F3_d', 'F3_e'], ['F3_d', 'F3_e'], 'provisional_reference_available'),
             ([], ['F2_c', 'F3_d'], 'parent_reported_duplication'),
             ([], ['F2_c'], 'no_modeled_nonfocal_sister'),
             ([], [], 'no_modeled_nonfocal_sister')]
    rows = []
    for i, (available, sequence, status) in enumerate(specs):
        assignments = [dict(gene=g, model_id='' if g == 'F2_c' else 'M-' + g[-1].upper(),
                            version='' if g == 'F2_c' else 10 if g == 'F3_d' else 6) for g in sequence]
        chosen = available[0] if available else ''
        rows.append(dict(family='OG0000000', taxon_id='F1', gene_node='n' + str(i), gene_a='F1_a', gene_b='F1_b',
                         status=status, nearest_reference_genes=json.dumps(available), chosen_reference_gene=chosen,
                         reference_model='M-D' if chosen else '', reference_version='10' if chosen else '',
                         sequence_first_nearest_genes=json.dumps(sequence), sequence_first_ties=str(len(sequence)),
                         sequence_first_tied_model_assignments=json.dumps(assignments),
                         sequence_first_chosen_gene=sequence[0] if sequence else '',
                         sequence_first_status='ineligible_parent_context' if i == 2 else 'fixture_choice',
                         parent_reported_duplication=str(int(i == 2)), support='0.8'))
    for guide in ['profile', 'mafft']:
        with (source / (guide + '_sequence_first_references.tsv')).open('w') as f:
            writer = csv.DictWriter(f, list(rows[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader(); writer.writerows(rows)
    tree = root / 'trees.txt'
    tree.write_text('OG0000000: ((F1_a,F1_b),(F2_c,(F3_d,F3_e)));\n')
    mappings = {str(species): sha(species), str(sequences): sha(sequences)}
    tree_proof = Path(write(root / 'tree-proof.json', dict(status='passed_complete_resolved_tree_membership_readback',
                                                         input_hashes=mappings, resolved_tree_file_sha256=sha(tree))))
    review = write(root / 'review-plan.json', dict(guides=[dict(guide=g, trees=str(tree), tree_readback=str(tree_proof))
                                                          for g in ['profile', 'mafft']]))
    reference = write(root / 'reference-plan.json', dict(bridge=str(bridge), review_plan=review))
    source_plan = write(root / 'source-plan.json', dict(output=str(source), reference_plan=reference))
    guides = [dict(guide=g, candidates=5) for g in ['profile', 'mafft']]
    source_receipt = Path(write(source / 'receipt.json', dict(
        status='complete_sequence_first_sister_reference_inventory_pending_readback',
        plan_sha256=sha(source_plan), guides=guides, source_hashes={},
        artifacts={p.name: sha(p) for p in source.iterdir()})))
    proof = Path(write(root / 'source-proof.json', dict(status='passed_full_sequence_first_sister_reference_native_readback',
        producer_receipt_sha256=sha(source_receipt), guides=guides, source_hashes={})))
    completion = write(root / 'source-completion.json', dict(status='complete_verified_expanded_sequence_first_references',
        summary=dict(guides=guides), source_hashes={str(source_receipt): sha(source_receipt), str(proof): sha(proof)}))
    streams = []
    for guide, pairs in [('profile', [(0, 1), (0, 3), (2, 3)]), ('mafft', [(0, 1), (0, 3), (1, 2)])]:
        directory = root / guide
        directory.mkdir()
        stream = directory / 'sorted_pairs.hex'
        stream.write_text(''.join(f'{a:06x}{b:06x}0\n{a:06x}{b:06x}1\n' for a, b in sorted(pairs)))
        snapshot = write(directory / 'snapshot.json', dict(status='complete_separate_small_family_ortholog_supplement',
                                                           input_hashes=mappings))
        identity = write(directory / 'identity-plan.json', dict(output=str(directory / 'identity'), supplement_receipt=snapshot,
                                                               pins={snapshot: sha(snapshot)}))
        id_receipt = Path(write(directory / 'identity' / 'receipt.json', dict(
            status='complete_grouped_ortholog_protein_family_identity_audit', plan_sha256=sha(identity))))
        audit_plan = write(directory / 'audit-plan.json', dict(identity_plan=identity, pins={identity: sha(identity)}))
        audit = write(directory / 'audit.json', dict(status='complete_native_ortholog_pair_multiplicity_audit',
            plan_sha256=sha(audit_plan), identity_receipt_sha256=sha(id_receipt), unique_unordered_pairs=len(pairs),
            duplicate_directed_incidences=0, pairs_missing_reverse=0, pairs_with_repeated_direction=0,
            pairs_with_unequal_multiplicity=0, sorted_pairs_sha256=sha(stream)))
        streams.append(dict(guide=guide, receipt=audit, stream=str(stream), audit_plan=audit_plan,
                            species_ids=str(species), sequence_ids=str(sequences)))
    return dict(source=str(source), source_plan=source_plan, source_readback=str(proof), source_completion=completion,
                bridge=str(bridge), streams=streams, compiler='/usr/bin/g++',
                cpp='scripts/query_ortholog_pair_membership.cpp', fixtures='scripts/check_ortholog_membership_cases.py',
                python=sys.executable, output=str(root / 'output'), pins={}, resources=dict(target_contexts=10))


def main():
    python = sys.executable
    with tempfile.TemporaryDirectory(prefix='reference-orthology-fixture-') as temp:
        root = Path(temp)
        plan = fixture(root)
        plan_path = write(root / 'plan.json', plan)
        producer = [python, 'scripts/query_reference_orthology.py', '--plan', plan_path]
        reader = [python, 'scripts/readback_reference_orthology.py', '--plan', plan_path, '--output']
        for command in [producer, reader + [str(root / 'readback.json')]]:
            run = subprocess.run(command, capture_output=True, text=True)
            if run.returncode:
                raise RuntimeError(run.stderr + run.stdout)
        out = Path(plan['output'])
        receipt = json.loads((out / 'receipt.json').read_text())
        assert receipt['target_contexts'] == 10 and receipt['unique_gene_pair_queries'] == 6
        assert receipt['reference_tie_records'] == 18 and receipt['duplicate_reference_links'] == 36
        # Expected matrices include parent-ineligible, missing models, absence,
        # cross-guide disagreement, nonlexical ties, duplicate physical queries.
        assert receipt['guides']['profile']['present'] == receipt['guides']['mafft']['present'] == 3
        originals = {p.name: p.read_bytes() for p in out.iterdir() if p.is_file()}
        rejected = []
        for label in ['cleared_parent_ineligibility', 'changed_lexical_choice', 'removed_unmodeled_context',
                      'changed_cross_guide_membership', 'deleted_nonlexical_link', 'changed_native_ordinal']:
            for name, data in originals.items():
                (out / name).write_bytes(data)
            if label in ['cleared_parent_ineligibility', 'changed_lexical_choice', 'removed_unmodeled_context']:
                rows = [json.loads(line) for line in (out / 'contexts.jsonl').read_text().splitlines()]
                if label == 'cleared_parent_ineligibility': rows[2]['parent_context_eligible'] = True
                elif label == 'changed_lexical_choice': rows[1]['designs']['sequence_first'][1]['lexical_choice'] = True
                else: rows.pop(3)
                (out / 'contexts.jsonl').write_text(''.join(json.dumps(row) + '\n' for row in rows))
            else:
                with (out / 'reference_membership.tsv').open() as f:
                    read = csv.DictReader(f, delimiter='\t'); fields=read.fieldnames; rows=list(read)
                if label == 'changed_cross_guide_membership': rows[1]['mafft_native_ortholog'] = '1'
                elif label == 'deleted_nonlexical_link': rows.pop(next(i for i,r in enumerate(rows) if r['lexical_choice'] == '0'))
                else: rows[0]['native_pair_key'] = '000002000003'
                with (out / 'reference_membership.tsv').open('w') as f:
                    w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
            altered = copy.deepcopy(receipt)
            altered['artifacts'] = {name: sha(out / name) for name in receipt['artifacts']}
            write(out / 'receipt.json', altered)
            run = subprocess.run(reader + [str(root / (label + '.json'))], capture_output=True, text=True)
            assert run.returncode != 0, label
            rejected.append(label)
        print(json.dumps(dict(status='passed_synthetic_full_reference_orthology_adapter_checks', target_contexts=10,
            unique_queries=6, logical_tie_records=18, duplicate_reference_links=36, rejected_rehashed_exports=rejected,
            scope='Synthetic fixture with stub native-tree/source qualification proofs; exercises actual compiled full '
                  'native stream merge and independent binary searches plus complete context/link reconstruction. '
                  'Not production source qualification, biological orthology or a pilot.'), indent=2))


if __name__ == '__main__':
    main()
