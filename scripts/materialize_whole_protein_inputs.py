#!/usr/bin/env python3
"""Materialize every audited unique whole-protein model input without fitting."""
import hashlib
import json
import subprocess
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import psutil
from assess_whole_protein_model_designs_v2 import matrix, KEYS
from inventory_whole_protein_model_inputs import signature
from screen_duplication_alignment_reuse import sha


def arrays_for(chosen, numeric, patterns):
    values = np.ascontiguousarray(numeric, dtype='<f8')
    assert np.isfinite(values).all() and np.all(values[:, 1] == 1)
    values[values == 0] = 0.
    assert chosen.groupby('background_id').family_component.nunique().eq(1).all()
    mapped = chosen.species_pattern_id.map(patterns)
    assert mapped.notna().all()
    np.testing.assert_array_equal(mapped.to_numpy(), chosen.species_pattern_row.to_numpy())
    return dict(matrix=values, row_identity=np.asarray(chosen.row_identity, dtype='S64'),
                background=pd.factorize(chosen.background_id, sort=True)[0].astype('<i8'),
                family=pd.factorize(chosen.family_component, sort=True)[0].astype('<i8'),
                pattern_rows=mapped.to_numpy(dtype='<i8'))


def main():
    pp = Path('metadata/whole_protein_materialized_inputs_plan_20260929.json')
    plan = json.loads(pp.read_text())
    bindings = {str(pp): sha(pp), **plan['pins']}
    launch = json.loads(Path(plan['audit_launch']).read_text())
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
    source_plan = json.loads(Path(plan['inventory_plan']).read_text())
    source = Path(source_plan['output'])
    receipt = json.loads((source / 'receipt.json').read_text())
    proof = json.loads(Path(plan['audit_proof']).read_text())
    assert proof['status'] == 'passed_full_whole_protein_input_inventory_readback'
    assert proof['source_receipt_sha256'] == sha(source / 'receipt.json')
    bindings.update(receipt['source_hashes'])
    bindings.update({str(source / n): h for n, h in receipt['artifacts'].items()})
    bindings[str(source / 'receipt.json')] = sha(source / 'receipt.json')
    bindings[plan['audit_proof']] = sha(plan['audit_proof'])
    def verify():
        for path, digest in bindings.items():
            assert sha(path) == digest, path
    verify()
    recipes = defaultdict(list)
    expected = set()
    for line in (source / 'unique_input_recipes.jsonl').open():
        row = json.loads(line)
        assert row['fit_input_id'] not in expected
        expected.add(row['fit_input_id'])
        recipes[row['partition_index']].append(row)
    assert len(expected) == receipt['unique_inputs']
    covariance = Path(source_plan['covariance_root'])
    cov = pd.read_csv(covariance / 'pair_covariance_index.tsv', sep='\t', usecols=[
        'target_id', 'background_id', 'row_identity', 'family_component', 'species_pattern_id', 'species_pattern_row']).rename(columns={'species_pattern_id': 'covariance_species_pattern_id'})
    covariance_sha = sha(covariance / 'receipt.json')
    patterns = pd.read_csv(plan['patterns'], sep='\t').set_index('species_pattern_id').row_index
    assert patterns.index.is_unique and patterns.is_unique
    selected = pd.read_csv(source_plan['selections'], sep='\t', usecols=['target_id', 'background_id', 'policy', 'scenario_id'])
    assert len(selected) == 2786912
    records_root = Path(source_plan['records_root'])
    out = Path(plan['output'])
    out.mkdir(exist_ok=False)
    seen, artifacts = set(), {}
    occurrences = byte_count = 0
    with (out / 'input_manifest.jsonl').open('w') as manifest:
        for part in json.loads((records_root / 'partition_manifest.json').read_text()):
            wanted = recipes.get(part['index'], [])
            if not wanted:
                continue
            assert sha(records_root / part['path']) == part['sha256']
            values = pd.read_parquet(records_root / part['path'])
            records = selected.merge(values, on=['target_id', 'background_id'], validate='many_to_one').merge(cov, on=['target_id', 'background_id'], validate='many_to_one').sort_values('target_id', kind='stable')
            assert records.species_pattern_id.equals(records.covariance_species_pattern_id)
            groups = records.groupby(KEYS).indices
            for recipe in wanted:
                frame = records.iloc[groups[tuple(recipe[k] for k in KEYS)]]
                chosen, x, names = matrix(frame, recipe['variant'])
                assert chosen.target_id.is_monotonic_increasing
                spec = recipe['specification']
                columns = spec['columns']
                assert columns[:2] == [recipe['outcome'], 'intercept']
                numeric = np.column_stack([chosen[recipe['outcome']].to_numpy(), np.ones(len(chosen)), x[:, [names.index(c) for c in columns[2:]]]])
                identifier, _, reconstructed = signature(chosen, numeric, columns, recipe['outcome'], covariance_sha)
                assert identifier == recipe['fit_input_id'] and reconstructed == spec and identifier not in seen
                arrays = arrays_for(chosen, numeric, patterns)
                assert hashlib.sha256(arrays['matrix'].tobytes()).hexdigest() == spec['values_sha256']
                assert hashlib.sha256(arrays['row_identity'].tobytes()).hexdigest() == spec['ordered_identity_sha256']
                folder = out / identifier[:2]
                folder.mkdir(exist_ok=True)
                path = folder / (identifier + '.npz')
                np.savez(path, **arrays)
                with np.load(path, allow_pickle=False) as saved:
                    assert set(saved.files) == set(arrays)
                    for key, value in arrays.items():
                        np.testing.assert_array_equal(saved[key], value)
                name = str(path.relative_to(out))
                digest = sha(path)
                artifacts[name] = digest
                manifest.write(json.dumps(dict(fit_input_id=identifier, recipe=recipe, path=name, sha256=digest), separators=(',', ':')) + '\n')
                seen.add(identifier)
                occurrences += len(chosen)
                byte_count += path.stat().st_size
            manifest.flush()
            print('Materialized whole-protein partition', part['index'], 'inputs', len(seen), '/', len(expected), flush=True)
    assert seen == expected
    verify()
    artifacts['input_manifest.jsonl'] = sha(out / 'input_manifest.jsonl')
    result = dict(status='complete_whole_protein_materialized_inputs_pending_readback', inputs=len(seen), tree_fits=len(seen)*5, record_occurrences=occurrences, npz_bytes=byte_count, source_hashes=bindings, script_sha256=sha(__file__), artifacts=artifacts, scope='All audited unique numerical inputs and covariance group indices serialized; same-process array roundtrip checked. Independent reconstruction required before fitting. No model fits or biological inference.')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'artifacts']}), flush=True)


if __name__ == '__main__':
    main()
