#!/usr/bin/env python3
"""Exercise complete-union, version, overlap and false-reuse rejection fixtures."""
import csv
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

from screen_duplication_alignment_reuse import sha


def save(path, record):
    Path(path).write_text(json.dumps(record) + '\n')


def main():
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        models = [dict(model_id=mid, version=version, sha256='coordinate-' + str(i),
                       sequence_sha256='sequence-' + str(i), length=10 + i)
                  for i, (mid, version) in enumerate([('M', 6), ('M', 10), ('Z', 1), ('Q', 1)])]

        def source(label, model_rows, endpoints):
            folder = root / label
            folder.mkdir()
            (folder / 'models.jsonl').write_text(''.join(json.dumps(m) + '\n' for m in model_rows))
            with (folder / 'model_pairs.tsv').open('w') as handle:
                writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
                writer.writerow(['pair_key', 'model_a', 'version_a', 'model_b', 'version_b'])
                for a, b in endpoints:
                    left, right = sorted([(models[a]['model_id'], models[a]['version']),
                                          (models[b]['model_id'], models[b]['version'])])
                    key = hashlib.sha256(json.dumps([left, right], separators=(',', ':')).encode()).hexdigest()
                    # Intentionally reverse original columns; canonical output must use numeric versions.
                    writer.writerow([key, *right, *left])
            save(folder / 'receipt.json', dict(status='complete_fixture', models=len(model_rows), pairs=len(endpoints),
                                               artifacts={name: sha(folder / name) for name in ['models.jsonl', 'model_pairs.tsv']}))
            return dict(label=label, root=str(folder), status='complete_fixture',
                        model_count_field='models', pair_count_field='pairs')

        original = source('old-original', models, [(0, 1), (0, 2)])
        changed_models = [dict(m, sha256='changed' if m['model_id'] == 'Z' else m['sha256']) for m in models]
        alternate = source('old-alternate', changed_models, [(0, 1)])
        current = source('new-first', changed_models, [(0, 1), (0, 2), (0, 3)])
        overlap = source('new-overlap', changed_models, [(0, 1)])
        plan = dict(old=[original, alternate], new=[current, overlap], output=str(root / 'output'), pins={})
        for spec in plan['old'] + plan['new']:
            for name in ['receipt.json', 'models.jsonl', 'model_pairs.tsv']:
                path = str(Path(spec['root']) / name)
                plan['pins'][path] = sha(path)
        pp = root / 'plan.json'
        save(pp, plan)
        producer = [sys.executable, 'scripts/inventory_alignment_reuse_candidates.py', '--plan', str(pp)]
        reader = [sys.executable, 'scripts/readback_alignment_reuse_candidates.py', '--plan', str(pp)]
        subprocess.run(producer, check=True, capture_output=True)
        proof = root / 'proof.json'
        checked = subprocess.run(reader + ['--output', str(proof)], capture_output=True, text=True)
        assert checked.returncode == 0, checked.stderr
        result = json.loads(proof.read_text())
        assert result['new_union_pairs'] == 3 and result['old_union_pairs'] == 2
        assert set(result['counts'].values()) == {1} and len(result['counts']) == 3
        folder = Path(plan['output'])
        table, receipt = folder / 'pair_reuse_candidates.tsv', folder / 'receipt.json'
        pristine = table.read_bytes(), receipt.read_bytes()
        rows = list(csv.DictReader(table.open(), delimiter='\t'))
        numeric = next(row for row in rows if row['model_a'] == row['model_b'] == 'M')
        assert (numeric['version_a'], numeric['version_b']) == ('6', '10')
        assert len(json.loads(numeric['new_sources'])) == len(json.loads(numeric['matching_old_sources'])) == 2

        def reject(number, edited):
            with table.open('w') as handle:
                writer = csv.DictWriter(handle, rows[0].keys(), delimiter='\t', lineterminator='\n')
                writer.writeheader()
                writer.writerows(edited)
            r = json.loads(receipt.read_text())
            r['artifacts'][table.name] = sha(table)
            save(receipt, r)
            out = root / f'rejected-{number}.json'
            proc = subprocess.run(reader + ['--output', str(out)], capture_output=True)
            assert proc.returncode != 0 and not out.exists(), proc.stdout
            table.write_bytes(pristine[0])
            receipt.write_bytes(pristine[1])

        false = [dict(row) for row in rows]
        false[0]['disposition'] = 'reuse_authorized'
        reject(1, false)
        reject(2, rows[:-1])
        reject(3, rows + [rows[0]])
        false = [dict(row) for row in rows]
        next(row for row in false if row['pair_key'] == numeric['pair_key'])['matching_old_sources'] = '["old-original"]'
        reject(4, false)
        # A model artifact with a refreshed receipt must still fail the frozen plan pin.
        path = Path(current['root']) / 'models.jsonl'
        path.write_text(path.read_text().replace('changed', 'tampered'))
        out = root / 'tampered-source.json'
        proc = subprocess.run(reader + ['--output', str(out)], capture_output=True)
        assert proc.returncode != 0 and not out.exists()
    print('PASS: complete overlapping unions, numeric versions 6/10, changed/new/matching sources; reject false dispositions, missing/duplicate rows, missing source labels and changed frozen sources.')


if __name__ == '__main__':
    main()
