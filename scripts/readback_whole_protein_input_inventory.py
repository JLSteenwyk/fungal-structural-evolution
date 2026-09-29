#!/usr/bin/env python3
"""Reconstruct every inventory entry without importing its matrix/hash helpers."""
import hashlib
import json
import subprocess
import time
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
import psutil
from screen_duplication_alignment_reuse import sha

KEYS = ['guide', 'policy', 'scenario_id']
SETTING = ['mask', 'cohort', 'screen', 'target_order', 'background_order']
OUTCOMES = ['rmsd_difference', 'mean_endpoint_tm_divergence_difference']
NUISANCE = ['original_coverage_difference', 'log_aligned_length_ratio', 'confidence_fraction_difference']


def reconstruct(frame, variant):
    if variant == 'positive_log_gene_distance_linear':
        frame = frame.loc[frame.target_sequence_distance.gt(0) & frame.background_sequence_distance.gt(0)]
    target = frame.target_sequence_distance.to_numpy()
    control = frame.background_sequence_distance.to_numpy()
    columns = {}
    if variant.startswith('gene_distance_'):
        degree = {'gene_distance_linear': 1, 'gene_distance_quadratic': 2, 'gene_distance_cubic': 3}[variant]
        for power in range(1, degree + 1):
            columns['gene_distance_power_' + str(power) + '_difference'] = target ** power - control ** power
    elif variant == 'positive_log_gene_distance_linear':
        columns['log_gene_distance_difference'] = np.log(target) - np.log(control)
    else:
        assert variant == 'alignment_identity_linear'
        columns['alignment_identity_difference'] = frame.identity_difference.to_numpy()
    for name in NUISANCE:
        columns[name] = frame[name].to_numpy()
    return frame, columns


def main():
    lp = Path('metadata/whole_protein_model_input_inventory_launch_20260928.json')
    launch = json.loads(lp.read_text())
    bindings = {str(lp): sha(lp), launch['plan']: launch['plan_sha256']}
    while psutil.pid_exists(launch['pid']):
        try:
            p = psutil.Process(launch['pid'])
            if abs(p.create_time() - launch['created']) > .01 or p.status() == psutil.STATUS_ZOMBIE:
                break
            assert p.cmdline() == launch['cmdline']
        except psutil.NoSuchProcess:
            break
        time.sleep(30)
    state = dict(x.split('=', 1) for x in subprocess.check_output(['systemctl', '--user', 'show', launch['unit'], '-p', 'ActiveState', '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
    assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0')
    plan = json.loads(Path(launch['plan']).read_text())
    root = Path(plan['output'])
    receipt = json.loads((root / 'receipt.json').read_text())
    assert receipt['status'] == 'complete_whole_protein_input_inventory_pending_independent_readback'
    bindings.update(receipt['source_hashes'])
    bindings[str(root / 'receipt.json')] = sha(root / 'receipt.json')
    bindings.update({str(root / name): digest for name, digest in receipt['artifacts'].items()})
    for path, digest in bindings.items():
        assert sha(path) == digest, path
    recipes = {}
    for line in (root / 'unique_input_recipes.jsonl').open():
        row = json.loads(line)
        assert row['fit_input_id'] not in recipes
        recipes[row['fit_input_id']] = row
    cov = pd.read_csv(Path(plan['covariance_root']) / 'pair_covariance_index.tsv', sep='\t')
    identities = {}
    for row in cov.to_dict('records'):
        identity = hashlib.sha256(json.dumps([row[k] for k in ['target_id', 'background_id', 'family_component', 'species_pattern_id']], separators=(',', ':')).encode()).hexdigest()
        assert identity == row['row_identity']
        key = (row['target_id'], row['background_id'])
        assert key not in identities
        identities[key] = identity
    assert len(identities) == 52675
    covariance_sha = sha(Path(plan['covariance_root']) / 'receipt.json')
    selected = pd.read_csv(plan['selections'], sep='\t', usecols=['target_id', 'background_id', 'policy', 'scenario_id'])
    assert len(selected) == 2786912
    source = Path(plan['records_root'])
    parts = json.loads((source / 'partition_manifest.json').read_text())
    counts, used = Counter(), set()
    total = 0
    with (root / 'setting_input_map.jsonl').open() as mapping, (Path(plan['design_root']) / 'designs.jsonl').open() as designs:
        for part in parts:
            values = pd.read_parquet(source / part['path'])
            assert sha(source / part['path']) == part['sha256']
            records = selected.merge(values, on=['target_id', 'background_id'], validate='many_to_one')
            groups = records.groupby(KEYS).indices
            for _ in range(2160):
                design = json.loads(next(designs))
                assert all(design[k] == part[k] for k in SETTING)
                key = tuple(design[k] for k in KEYS)
                frame = records.iloc[groups.get(key, [])].sort_values('target_id', kind='stable')
                assert len(frame) == design['screen_retained_records']
                chosen, columns = reconstruct(frame, design['variant'])
                assert len(chosen) == design['records'] and list(columns) == design['columns']
                for outcome in OUTCOMES:
                    actual = json.loads(next(mapping))
                    expected = {k: design[k] for k in KEYS + SETTING}
                    expected.update(variant=design['variant'], outcome=outcome, design_status=design['status'], records=len(chosen), fit_input_id=None, partition_index=part['index'])
                    if design['status'] == 'full_rank_with_positive_residual_df':
                        names = [outcome, 'intercept'] + design['active_columns']
                        array = np.empty((len(chosen), len(names)), dtype='<f8')
                        array[:, 0] = chosen[outcome]
                        array[:, 1] = 1.0
                        for i, name in enumerate(names[2:], 2):
                            array[:, i] = columns[name]
                        assert np.isfinite(array).all()
                        array[array == 0] = 0.0
                        identity_bytes = b''.join(identities[pair].encode('ascii') for pair in chosen[['target_id', 'background_id']].itertuples(index=False, name=None))
                        spec = dict(records=len(chosen), columns=names, outcome=outcome,
                                    values_sha256=hashlib.sha256(array.tobytes(order='C')).hexdigest(),
                                    ordered_identity_sha256=hashlib.sha256(identity_bytes).hexdigest(),
                                    covariance_receipt_sha256=covariance_sha,
                                    model='equal-record intercept model; residual/background/family-component/species working covariance; five trees separately')
                        fit_id = hashlib.sha256(json.dumps(spec, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
                        expected['fit_input_id'] = fit_id
                        assert recipes[fit_id]['specification'] == spec
                        if fit_id not in used:
                            assert recipes[fit_id] == dict(expected, specification=spec)
                        used.add(fit_id)
                    assert actual == expected
                    counts[design['status']] += 1
                    total += 1
            print('Verified whole-protein inventory setting', part['index'] + 1, '/96', flush=True)
        assert mapping.read() == designs.read() == ''
    assert total == receipt['settings'] == 414720 and dict(counts) == receipt['design_status_counts']
    assert used == set(recipes) and len(used) == receipt['unique_inputs']
    assert receipt['tree_alternatives'] == 5 and receipt['unique_tree_fits'] == len(used) * 5
    for path, digest in bindings.items():
        assert sha(path) == digest, path
    result = dict(status='passed_full_whole_protein_input_inventory_readback', source_receipt_sha256=sha(root / 'receipt.json'), checker_sha256=sha(__file__), producer_terminal_state=state, settings=total, unique_inputs=len(used), scope='All settings and every unique representative recipe reconstructed independently; exact bytes and row identities checked. Polynomial/log operation order intentionally retained for byte equality; alternate algebra was checked by the preceding design audit. No fitted-effect or calibrated inference claim.')
    with Path('metadata/whole_protein_input_inventory_completed_readback_20260928.json').open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
