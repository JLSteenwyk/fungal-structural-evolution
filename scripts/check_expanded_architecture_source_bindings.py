#!/usr/bin/env python3
"""Exercise full source binding and rejection around the architecture grid."""
import csv
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value) + '\n')


def table(path, rows):
    with Path(path).open('w') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    policies = ['alignment_evalue', 'alignment_bitscore', 'envelope_evalue', 'envelope_bitscore']
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        annotations = []
        targets = []
        cases = [('shared_conservative_architecture', [1], [1], True),
                 ('neither_annotated', [], [], False),
                 ('one_unannotated', [1], [], True),
                 ('different_ordered_annotations', [1, 2], [2, 1], True),
                 ('shared_but_nonconservative', [1], [1], False)]
        for i, (_, left, right, conservative) in enumerate(cases):
            for suffix, hits in [('a', left), ('b', right)]:
                annotations.append(dict(model_id=f'M{i}{suffix}', version=1, policies={
                    p: [dict(pfam_accession=f'PF{x}.1', pfam_type='Domain',
                             conservative_architecture=int(conservative)) for x in hits] for p in policies}))
            targets.append(dict(guide='profile', family='OG', taxon_id='T0', gene_a=f'G{i}a',
                                gene_b=f'G{i}b', gene_node=f'n{i}', model_a=f'M{i}a', version_a=1,
                                model_b=f'M{i}b', version_b=1))
        for suffix in ['a', 'b']:
            annotations.append(dict(model_id='B' + suffix, version=1, policies={
                p: [dict(pfam_accession='PF1.1', pfam_type='Domain', conservative_architecture=1)]
                for p in policies}))
        annotation_path = root / 'annotations.jsonl'
        annotation_path.write_text(''.join(json.dumps(r) + '\n' for r in annotations))
        targets_path = root / 'targets.tsv'
        table(targets_path, targets)
        background = dict(model_id_a='Ba', version_a=1, model_id_b='Bb', version_b=1,
                          taxon_a='T0', taxon_b='T1', candidate_and_native_ortholog_both=1,
                          both_guides_and_parents_unreported=1, both_parents_unreported=1)
        for guide in ['profile', 'mafft']:
            background.update({guide + '_candidate_and_native_ortholog': 1,
                               guide + '_native_ortholog': 1,
                               guide + '_candidate_status': 'cross_taxon_unreported_candidate',
                               guide + '_family': 'OG', guide + '_sequence_pair_distance': 1})
        background_path = root / 'background.tsv'
        table(background_path, [background])
        settings = ['guide_native_ortholog', 'both_guides_native_ortholog', 'both_guides_unreported_parents']
        support = []
        for target in targets:
            for setting in settings:
                row = {k: target[k] for k in ['guide', 'family', 'taxon_id', 'gene_a', 'gene_b', 'gene_node']}
                row.update(target_same_model=0, target_sequence_distance=1, background_set=setting,
                           distance_rule='multiplicative_distance_range')
                for factor in ['1_25', '1_5', '2_0']:
                    row['within_factor_' + factor] = row['focal_within_factor_' + factor] = 1
                support.append(row)
        support_path = root / 'support.tsv'
        table(support_path, support)
        receipt_path = root / 'source_receipt.json'
        source = dict(artifacts={p.name: sha(p) for p in [annotation_path, targets_path, background_path, support_path]})
        save(receipt_path, source)
        proof_path = root / 'source_proof.json'
        save(proof_path, dict(status='fixture_verified_source', producer_receipt_sha256=sha(receipt_path)))
        output = root / 'output'
        plan_path = root / 'plan.json'
        plan = dict(verified_sources=[dict(receipt=str(receipt_path), readback=str(proof_path),
                                          readback_status='fixture_verified_source')],
                    annotation_tables=[str(annotation_path)], target_links=str(targets_path),
                    background_links=str(background_path), support_table=str(support_path),
                    output=str(output), pins={})
        save(plan_path, plan)
        producer = [sys.executable, 'scripts/assess_architecture_matched_support_v2.py', '--plan']
        subprocess.run(producer + [str(plan_path)], check=True, capture_output=True, text=True)
        records = list(csv.DictReader((output / 'architecture_support.tsv').open(), delimiter='\t'))
        assert len(records) == 60
        for record in records:
            i = int(record['gene_node'][1:])
            assert record['target_architecture_status'] == cases[i][0]
            assert record['within_factor_2_0'] == str(int(i == 0))
        config_path = root / 'config.json'
        config = dict(source_plan=str(plan_path), output=str(root / 'readback.json'), pins={})
        save(config_path, config)
        checker = [sys.executable, 'scripts/readback_architecture_support_v2.py', '--config']
        subprocess.run(checker + [str(config_path)], check=True, capture_output=True, text=True)
        # A future source artifact is absent from static plan pins, but must
        # still fail both stages if it changes relative to its receipt.
        original = background_path.read_bytes()
        background_path.write_bytes(original + b'\n')
        plan['output'] = str(root / 'bad-source-output')
        bad_plan = root / 'bad-plan.json'
        save(bad_plan, plan)
        assert subprocess.run(producer + [str(bad_plan)], capture_output=True).returncode != 0
        config['output'] = str(root / 'bad-source-readback.json')
        save(config_path, config)
        assert subprocess.run(checker + [str(config_path)], capture_output=True).returncode != 0
        background_path.write_bytes(original)
        # Rehashing a false zero-support output must not evade full arithmetic
        # reconstruction merely because its artifact checksum now agrees.
        records[0]['within_factor_2_0'] = '0'
        table(output / 'architecture_support.tsv', records)
        result = json.loads((output / 'receipt.json').read_text())
        result['artifacts']['architecture_support.tsv'] = sha(output / 'architecture_support.tsv')
        save(output / 'receipt.json', result)
        config['output'] = str(root / 'bad-arithmetic-readback.json')
        save(config_path, config)
        assert subprocess.run(checker + [str(config_path)], capture_output=True).returncode != 0
        assert not Path(config['output']).exists()
    print('All five architecture statuses, full 60-row grid, missing-static-pin source changes rejected by both stages, and rehashed false support rejected.')


if __name__ == '__main__':
    main()
