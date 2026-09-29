#!/usr/bin/env python3
"""Inventory exact whole-protein likelihood inputs after independent design audit."""
import argparse
import hashlib
import json
import subprocess
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import psutil
from assess_whole_protein_model_designs_v2 import matrix, KEYS, SETTING, VARIANTS
from screen_duplication_alignment_reuse import sha

OUTCOMES = ['rmsd_difference', 'mean_endpoint_tm_divergence_difference']


def signature(frame, numeric, columns, outcome, covariance_sha):
    values = np.ascontiguousarray(numeric, dtype='<f8')
    assert np.isfinite(values).all() and values.shape == (len(frame), len(columns))
    values[values == 0] = 0.0
    identities = np.asarray(frame.row_identity, dtype='S64')
    specification = dict(records=len(frame), columns=columns, outcome=outcome,
                         values_sha256=hashlib.sha256(values.tobytes()).hexdigest(),
                         ordered_identity_sha256=hashlib.sha256(identities.tobytes()).hexdigest(),
                         covariance_receipt_sha256=covariance_sha,
                         model='equal-record intercept model; residual/background/family-component/species working covariance; five trees separately')
    encoded = json.dumps(specification, sort_keys=True, separators=(',', ':'))
    return hashlib.sha256(encoded.encode()).hexdigest(), encoded, specification


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--plan', type=Path, required=True)
    args = ap.parse_args()
    plan = json.loads(args.plan.read_text())
    bindings = {str(args.plan): sha(args.plan), **plan['pins']}
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
    state = dict(x.split('=', 1) for x in subprocess.check_output(
        ['systemctl', '--user', 'show', launch['unit'], '-p', 'ActiveState', '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
    assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0')
    design = Path(plan['design_root'])
    proof = json.loads(Path(plan['audit_proof']).read_text())
    assert proof['status'] == 'passed_full_whole_protein_model_design_readback'
    assert proof['source_receipt_sha256'] == sha(design / 'receipt.json')
    dr = json.loads((design / 'receipt.json').read_text())
    bindings.update(dr['source_hashes'])
    for name, digest in dr['artifacts'].items():
        bindings[str(design / name)] = digest
    for path in [plan['audit_proof'], str(design / 'receipt.json')]:
        bindings[path] = sha(path)
    covariance = Path(plan['covariance_root'])
    cr = json.loads((covariance / 'receipt.json').read_text())
    covariance_sha = sha(covariance / 'receipt.json')
    assert cr['status'] == 'complete_whole_protein_covariance_index_with_full_serialized_readback'
    bindings.update(cr['source_hashes'])
    for name, digest in cr['artifacts'].items():
        bindings[str(covariance / name)] = digest
    def verify():
        for path, digest in bindings.items():
            assert sha(path) == digest, path
    verify()
    cov = pd.read_csv(covariance / 'pair_covariance_index.tsv', sep='\t', usecols=['target_id', 'background_id', 'row_identity'])
    assert len(cov) == 52675 and not cov.duplicated(['target_id', 'background_id']).any()
    selected = pd.read_csv(plan['selections'], sep='\t', usecols=['target_id', 'background_id', 'policy', 'scenario_id'])
    assert len(selected) == 2786912 and not selected.duplicated(['target_id', 'policy', 'scenario_id']).any()
    source = Path(plan['records_root'])
    parts = json.loads((source / 'partition_manifest.json').read_text())
    out = Path(plan['output'])
    out.mkdir(exist_ok=False)
    recipes, encoded_signatures, counts = {}, {}, Counter()
    total = 0
    with (design / 'designs.jsonl').open() as designs, (out / 'setting_input_map.jsonl').open('w') as mapping:
        for part in parts:
            values = pd.read_parquet(source / part['path'])
            assert sha(source / part['path']) == part['sha256']
            records = selected.merge(values, on=['target_id', 'background_id'], validate='many_to_one')
            records = records.merge(cov, on=['target_id', 'background_id'], validate='many_to_one').sort_values('target_id', kind='stable')
            groups = records.groupby(KEYS).indices
            seen = set()
            for _ in range(432 * 5):
                row = json.loads(next(designs))
                key = tuple(row[k] for k in KEYS)
                assert all(row[k] == part[k] for k in SETTING) and (*key, row['variant']) not in seen
                seen.add((*key, row['variant']))
                assert row['variant'] in VARIANTS
                frame = records.iloc[groups.get(key, [])]
                assert len(frame) == row['screen_retained_records']
                chosen, x, names = matrix(frame, row['variant'])
                assert len(chosen) == row['records'] and names == row['columns']
                for outcome in OUTCOMES:
                    entry = {k: row[k] for k in KEYS + SETTING}
                    entry.update(variant=row['variant'], outcome=outcome, design_status=row['status'], records=len(chosen), fit_input_id=None, partition_index=part['index'])
                    if row['status'] == 'full_rank_with_positive_residual_df':
                        assert chosen.target_id.is_monotonic_increasing
                        active = [names.index(n) for n in row['active_columns']]
                        numeric = np.column_stack([chosen[outcome].to_numpy(), np.ones(len(chosen)), x[:, active]])
                        columns = [outcome, 'intercept'] + row['active_columns']
                        fit_id, encoded, spec = signature(chosen, numeric, columns, outcome, covariance_sha)
                        if fit_id in recipes:
                            assert encoded_signatures[fit_id] == encoded
                        else:
                            encoded_signatures[fit_id] = encoded
                            recipes[fit_id] = dict(entry, fit_input_id=fit_id, specification=spec)
                        entry['fit_input_id'] = fit_id
                    counts[row['status']] += 1
                    mapping.write(json.dumps(entry, separators=(',', ':')) + '\n')
                    total += 1
            assert len(seen) == 2160
            print('Inventoried whole-protein setting', part['index'] + 1, '/96', flush=True)
        assert designs.read() == '' and total == 414720
    with (out / 'unique_input_recipes.jsonl').open('w') as handle:
        for key in sorted(recipes):
            handle.write(json.dumps(recipes[key], sort_keys=True, separators=(',', ':')) + '\n')
    verify()
    receipt = dict(status='complete_whole_protein_input_inventory_pending_independent_readback',
                   source_hashes=bindings, script_sha256=sha(__file__), settings=total,
                   outcomes=OUTCOMES, design_status_counts=dict(counts), unique_inputs=len(recipes),
                   tree_alternatives=5, unique_tree_fits=len(recipes) * 5,
                   artifacts={p.name: sha(p) for p in out.iterdir()},
                   scope='All design/outcome settings retained, with explicit non-estimable statuses. Only exact ordered identities and numeric input bytes share recipes; signed zero normalized. Five trees remain separate. Inventory only: no fits, runtime estimate, calibrated uncertainty or biological inference. Full independent readback required before reuse.')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: v for k, v in receipt.items() if k != 'source_hashes'}), flush=True)


if __name__ == '__main__':
    main()
