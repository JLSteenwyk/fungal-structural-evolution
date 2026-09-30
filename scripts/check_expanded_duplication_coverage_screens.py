#!/usr/bin/env python3
"""Exercise original-length filtering, excluded orders and full missing/same-model retention."""
import csv
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

from screen_duplication_alignment_reuse import sha


def main():
    with tempfile.TemporaryDirectory() as name:
        root = Path(name)
        queue, summary = root / 'queue', root / 'summary'
        queue.mkdir()
        summary.mkdir()

        def js(path, value):
            path.write_text(json.dumps(value, indent=2) + '\n')

        def table(path, records, compressed=False):
            op = gzip.open(path, 'wt') if compressed else path.open('w')
            with op as handle:
                writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter='\t', lineterminator='\n')
                writer.writeheader()
                writer.writerows(records)

        models = [dict(model_id='A', version=1, length=100), dict(model_id='B', version=1, length=200),
                  dict(model_id='C', version=1, length=100)]
        (queue / 'models.jsonl').write_text(''.join(json.dumps(row) + '\n' for row in models))
        pairs = []
        for a, b in [('A', 'B'), ('A', 'C')]:
            key = hashlib.sha256(json.dumps([(a, 1), (b, 1)], separators=(',', ':')).encode()).hexdigest()
            pairs.append(dict(pair_key=key, model_a=a, version_a=1, model_b=b, version_b=1))
        table(queue / 'model_pairs.tsv', pairs)
        numbers = []
        for ix, pair in enumerate(pairs):
            for mask in ['full', 'plddt70']:
                row = dict(pair_key=pair['pair_key'], mask=mask)
                for order in [0, 1]:
                    excluded = ix == 1 and order == 0
                    n = 100 if mask == 'full' else 50
                    row.update({f'order{order}_status': 'excluded' if excluded else 'aligned',
                                f'order{order}_native_status': 'aligned',
                                f'order{order}_numerical_exclusion_reasons': 'rmsd_discrepancy' if excluded else '',
                                f'order{order}_aligned_length': '' if excluded else str(float(n))})
                numbers.append(row)
        table(summary / 'pair_mask_order_summary.tsv', numbers)
        js(summary / 'receipt.json', dict(artifacts={'pair_mask_order_summary.tsv': sha(summary / 'pair_mask_order_summary.tsv')}))
        proof = root / 'summary_audit.json'
        js(proof, dict(status='passed_full_primary_usable_order_summary_readback',
                       source_receipt_sha256=sha(summary / 'receipt.json'), pair_mask_rows=4))
        events, links = [], []
        for ix, a, b, category in [(0, 'A', 'B', 'both_models'), (1, 'A', 'C', 'both_models'),
                                  (2, 'A', 'A', 'both_models'), (3, 'A', '', 'one_model'),
                                  (4, '', '', 'neither_model')]:
            event = dict(guide='profile', family='family', taxon_id='T1', gene_node='n' + str(ix),
                         gene_a='ga' + str(ix), gene_b='gb' + str(ix), new_model_a=a, new_model_b=b,
                         new_coverage_class=category, new_same_model=int(a == b and bool(a)))
            events.append(event)
            if category == 'both_models':
                key = hashlib.sha256(json.dumps(sorted([(a, 1), (b, 1)]), separators=(',', ':')).encode()).hexdigest()
                links.append({**{k: event[k] for k in ['guide', 'family', 'taxon_id', 'gene_node', 'gene_a', 'gene_b']},
                              'model_a': a, 'model_b': b, 'version_a': '1', 'version_b': '1',
                              'same_model': str(event['new_same_model']), 'pair_key': key,
                              'comparison_status': 'identical_model_no_alignment' if a == b else 'queued_distinct_models'})
        table(queue / 'event_model_pair_links.tsv', links)
        js(queue / 'receipt.json', dict(status='complete_reviewed_duplication_model_pair_queue',
                                       unique_models=3, unique_distinct_model_pairs=2,
                                       counts=dict(event_links=3),
                                       artifacts={p.name: sha(p) for p in queue.iterdir()}))
        event_path = root / 'events.tsv.gz'
        table(event_path, events, compressed=True)
        manifest = root / 'manifest.tsv'
        table(manifest, [dict(taxon_id=t, species_name=t, study_role='ingroup', lineage='lineage') for t in ['T1', 'T2']])
        screens_plan = root / 'screens.json'
        screens = json.loads(Path('metadata/whole_protein_common_fits_plan_20260927.json').read_text())['screens']
        js(screens_plan, dict(screens=screens))
        plan_path = root / 'plan.json'
        out = root / 'output'
        js(plan_path, dict(queue=str(queue), summary=str(summary), summary_audit=str(proof),
                           events=str(event_path), manifest=str(manifest), screens_plan=str(screens_plan),
                           screens=screens, output=str(out), pins={},
                           expected=dict(models=3, pairs=2, pair_mask_rows=4, modeled_events=3, events=5, taxa=2)))
        subprocess.run([sys.executable, 'scripts/screen_expanded_duplication_coverage.py', '--plan', str(plan_path)],
                       check=True, capture_output=True, text=True)
        rows = list(csv.DictReader((out / 'pair_mask_coverage.tsv').open(), delimiter='\t'))
        assert rows[0]['n50_c50_pass'] == '1'
        assert rows[1]['n50_c50_pass'] == '0' and rows[1]['order0_original_coverage_b'] == '0.25'
        assert rows[2]['n30_c50_exclusions'] == 'order0_not_numerically_usable'
        with gzip.open(out / 'event_mask_coverage.tsv.gz', 'rt') as handle:
            records = list(csv.DictReader(handle, delimiter='\t'))
        assert len(records) == 10 and {r['measurement_disposition'] for r in records} == {
            'new_model_pair', 'identical_model_no_alignment', 'one_model', 'no_models'}
        zero = [r for r in csv.DictReader((out / 'taxon_screen_coverage.tsv').open(), delimiter='\t') if r['taxon_id'] == 'T2']
        assert len(zero) == 24 and all(r['events'] == r['passed_events'] == '0' for r in zero)
        command = [sys.executable, 'scripts/readback_expanded_duplication_coverage_screens.py', '--plan', str(plan_path)]
        subprocess.run(command + ['--output', str(root / 'readback.json')], check=True, capture_output=True, text=True)
        records[1]['n50_c50_pass'] = '1'
        records[1]['n50_c50_exclusions'] = ''
        table(out / 'event_mask_coverage.tsv.gz', records, compressed=True)
        receipt = json.loads((out / 'receipt.json').read_text())
        receipt['artifacts']['event_mask_coverage.tsv.gz'] = sha(out / 'event_mask_coverage.tsv.gz')
        js(out / 'receipt.json', receipt)
        bad = subprocess.run(command + ['--output', str(root / 'bad-readback.json')], capture_output=True, text=True)
        assert bad.returncode != 0 and not (root / 'bad-readback.json').exists()
    print('Original-length rather than retained-mask denominator, exact threshold, excluded-order retention, missing/same-model dispositions, full zero-taxon grid and rehashed false-pass rejection passed.')


if __name__ == '__main__':
    main()
