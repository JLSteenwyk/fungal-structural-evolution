#!/usr/bin/env python3
"""Rebuild every serialized whole-protein input and its covariance group labels."""
import hashlib
import json
import subprocess
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import psutil
from readback_whole_protein_input_inventory import reconstruct, KEYS
from screen_duplication_alignment_reuse import sha


def codes(values):
    labels = {label: i for i, label in enumerate(sorted(set(values)))}
    return np.array([labels[label] for label in values], dtype='<i8')


def check_arrays(saved, numeric, identities, backgrounds, families, pattern_rows):
    expected = dict(matrix=numeric, row_identity=np.asarray(identities, dtype='S64'),
                    background=codes(backgrounds), family=codes(families),
                    pattern_rows=np.asarray(pattern_rows, dtype='<i8'))
    assert set(saved) == set(expected)
    for name, value in expected.items():
        assert saved[name].dtype == value.dtype
        np.testing.assert_array_equal(saved[name], value)


def main():
    pp = Path('metadata/whole_protein_materialized_input_readback_plan_20260929.json')
    plan = json.loads(pp.read_text())
    bindings = {str(pp): sha(pp), **plan['pins']}
    launch = json.loads(Path(plan['launch']).read_text())
    while psutil.pid_exists(launch['pid']):
        try:
            p = psutil.Process(launch['pid'])
            if abs(p.create_time() - launch['created']) > .01 or p.status() == psutil.STATUS_ZOMBIE:
                break
            assert p.cmdline() == launch['cmdline']
        except psutil.NoSuchProcess:
            break
        time.sleep(30)
    state = dict(line.split('=', 1) for line in subprocess.check_output(
        ['systemctl', '--user', 'show', launch['unit'], '-p', 'ActiveState', '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
    assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0')
    producer = json.loads(Path(launch['plan']).read_text())
    source_plan = json.loads(Path(producer['inventory_plan']).read_text())
    root = Path(producer['output'])
    receipt = json.loads((root / 'receipt.json').read_text())
    assert receipt['status'] == 'complete_whole_protein_materialized_inputs_pending_readback'
    bindings.update(receipt['source_hashes'])
    bindings[str(root / 'receipt.json')] = sha(root / 'receipt.json')
    bindings.update({str(root / name): h for name, h in receipt['artifacts'].items()})
    def verify():
        for path, digest in bindings.items():
            assert sha(path) == digest, path
    verify()
    inventory = Path(source_plan['output'])
    recipes = {}
    for line in (inventory / 'unique_input_recipes.jsonl').open():
        row = json.loads(line)
        assert row['fit_input_id'] not in recipes
        recipes[row['fit_input_id']] = row
    wanted = defaultdict(list)
    seen_manifest = set()
    for line in (root / 'input_manifest.jsonl').open():
        row = json.loads(line)
        key = row['fit_input_id']
        assert key not in seen_manifest and row['recipe'] == recipes[key]
        assert row['sha256'] == receipt['artifacts'][row['path']]
        seen_manifest.add(key)
        wanted[row['recipe']['partition_index']].append(row)
    assert seen_manifest == set(recipes)
    covroot = Path(source_plan['covariance_root'])
    cov = pd.read_csv(covroot / 'pair_covariance_index.tsv', sep='\t')
    lookup = {}
    for row in cov.to_dict('records'):
        key = row['target_id'], row['background_id']
        assert key not in lookup
        identity = hashlib.sha256(json.dumps([*key, row['family_component'], row['species_pattern_id']], separators=(',', ':')).encode()).hexdigest()
        assert identity == row['row_identity']
        lookup[key] = row
    patterns = pd.read_csv(producer['patterns'], sep='\t')
    assert patterns.species_pattern_id.is_unique and patterns.row_index.is_unique
    pattern_map = dict(zip(patterns.species_pattern_id, patterns.row_index))
    selected = pd.read_csv(source_plan['selections'], sep='\t', usecols=['target_id', 'background_id', 'policy', 'scenario_id'])
    source = Path(source_plan['records_root'])
    seen = set()
    occurrences = byte_count = 0
    for part in json.loads((source / 'partition_manifest.json').read_text()):
        if not wanted.get(part['index']):
            continue
        assert sha(source / part['path']) == part['sha256']
        records = selected.merge(pd.read_parquet(source / part['path']), on=['target_id', 'background_id'], validate='many_to_one')
        groups = records.groupby(KEYS).indices
        for item in wanted[part['index']]:
            recipe = item['recipe']
            frame = records.iloc[groups[tuple(recipe[k] for k in KEYS)]].sort_values('target_id', kind='stable')
            chosen, columns = reconstruct(frame, recipe['variant'])
            spec = recipe['specification']
            assert spec['columns'][:2] == [recipe['outcome'], 'intercept']
            numeric = np.empty((len(chosen), len(spec['columns'])), dtype='<f8')
            numeric[:, 0] = chosen[recipe['outcome']]
            numeric[:, 1] = 1.
            for j, name in enumerate(spec['columns'][2:], 2):
                numeric[:, j] = columns[name]
            assert np.isfinite(numeric).all()
            numeric[numeric == 0] = 0.
            rows = [lookup[pair] for pair in chosen[['target_id', 'background_id']].itertuples(index=False, name=None)]
            assert [r['species_pattern_id'] for r in rows] == list(chosen.species_pattern_id)
            identities = [r['row_identity'] for r in rows]
            assert hashlib.sha256(numeric.tobytes()).hexdigest() == spec['values_sha256']
            assert hashlib.sha256(''.join(identities).encode('ascii')).hexdigest() == spec['ordered_identity_sha256']
            assert len(rows) == spec['records'] == recipe['records']
            assert spec['covariance_receipt_sha256'] == sha(covroot / 'receipt.json')
            assert hashlib.sha256(json.dumps(spec, sort_keys=True, separators=(',', ':')).encode()).hexdigest() == item['fit_input_id']
            mapped = [pattern_map[r['species_pattern_id']] for r in rows]
            assert mapped == [r['species_pattern_row'] for r in rows]
            path = root / item['path']
            with np.load(path, allow_pickle=False) as saved:
                check_arrays(saved, numeric, identities, list(chosen.background_id), [r['family_component'] for r in rows], mapped)
            assert item['fit_input_id'] not in seen
            seen.add(item['fit_input_id'])
            occurrences += len(rows)
            byte_count += path.stat().st_size
        print('Checked materialized whole-protein partition', part['index'], 'inputs', len(seen), '/', len(recipes), flush=True)
    assert seen == set(recipes) and len(seen) == receipt['inputs']
    assert occurrences == receipt['record_occurrences'] and byte_count == receipt['npz_bytes']
    assert receipt['tree_fits'] == len(seen)*5
    verify()
    proof = dict(status='passed_full_whole_protein_materialized_input_readback', inputs=len(seen), record_occurrences=occurrences, npz_bytes=byte_count, source_receipt_sha256=sha(root / 'receipt.json'), checker_sha256=sha(__file__), producer_terminal_state=state, scope='Every numerical matrix, ordered row identity, background/family code and species-pattern index independently reconstructed from audited source records. Full byte/hash bindings checked before and after. No fits, calibrated uncertainty or biological inference.')
    with Path(plan['proof']).open('x') as handle:
        json.dump(proof, handle, indent=2)
        handle.write('\n')
    print(json.dumps(proof), flush=True)


if __name__ == '__main__':
    main()
