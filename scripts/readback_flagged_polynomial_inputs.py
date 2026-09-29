#!/usr/bin/env python3
"""Check all reconstructed flagged inputs against audited recipes and source identities."""
import csv
import hashlib
import json
import subprocess
import time
from pathlib import Path
import numpy as np
import pandas as pd
import psutil
from screen_duplication_alignment_reuse import sha


def codes(values):
    lookup = {value: i for i, value in enumerate(sorted(set(values)))}
    return np.asarray([lookup[v] for v in values], dtype=np.int64)


def main():
    lp = Path('metadata/flagged_polynomial_input_preparation_launch_20260928.json')
    launch = json.loads(lp.read_text())
    lh = sha(lp)
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
    root = Path('results/model_validation/flagged-polynomial-inputs-20260928-v1')
    r = json.loads((root / 'receipt.json').read_text())
    rh = sha(root / 'receipt.json')
    assert r['status'] == 'complete_flagged_polynomial_input_reconstruction_pending_independent_readback'
    bindings = dict(r['source_hashes'])
    bindings.update({str(root / n): h for n, h in r['artifacts'].items()})
    for path, digest in bindings.items():
        assert sha(path) == digest, path
    plan = json.loads(Path('metadata/full_polynomial_ml_plan_20260927.json').read_text())
    recipes = {}
    for folder in ['inventory', 'nonlinear_inventory']:
        for line in (Path(plan[folder]) / 'unique_fit_recipes.jsonl').open():
            row = json.loads(line)
            assert row['fit_input_id'] not in recipes
            recipes[row['fit_input_id']] = dict(row, polynomial_degree=row.get('polynomial_degree', 1))
    components = {row['node_id']: row['family_component'] for row in csv.DictReader(Path(plan['nodes']).open(), delimiter='\t') if row['role'] == 'target'}
    patterns = {(row['target_id'], row['background_id']): row['species_pattern_id'] for row in csv.DictReader(Path(plan['pairs']).open(), delimiter='\t')}
    pattern_rows = {row['species_pattern_id']: int(row['row_index']) for row in csv.DictReader((Path(plan['factors']) / 'patterns.tsv').open(), delimiter='\t')}
    selected = pd.read_csv(plan['selections'], sep='\t', usecols=['target_id', 'background_id', 'policy', 'scenario_id', 'domain_config_id'])
    guides = {row['node_id']: row['guide'] for row in csv.DictReader(Path(plan['nodes']).open(), delimiter='\t') if row['role'] == 'target'}
    selected['guide'] = selected.target_id.map(guides)
    assert selected.guide.notna().all() and len(selected) == 2786912
    parts = {p['index']: p for p in json.loads((Path(plan['summaries']) / 'partition_manifest.json').read_text())}
    manifest = json.loads((root / 'input_manifest.json').read_text())
    checked, occurrences = set(), 0
    for index in sorted({row['recipe']['partition_index'] for row in manifest}):
        part = parts[index]
        path = Path(plan['summaries']) / part['path']
        assert sha(path) == part['sha256']
        valid = set(pd.read_parquet(path, columns=['domain_config_id']).domain_config_id)
        members = selected.loc[selected.domain_config_id.isin(valid)]
        groups = members.groupby(['guide', 'policy', 'scenario_id']).indices
        for entry in manifest:
            recipe = entry['recipe']
            if recipe['partition_index'] != index:
                continue
            identifier = recipe['fit_input_id']
            assert identifier not in checked and recipe == recipes[identifier]
            frame = members.iloc[groups[tuple(recipe[k] for k in ['guide', 'policy', 'scenario_id'])]].sort_values('target_id')
            assert len(frame) == recipe['records'] and not frame.target_id.duplicated().any()
            targets = frame.target_id.tolist()
            backgrounds = frame.background_id.tolist()
            families = [components[t] for t in targets]
            species = [patterns[(t, b)] for t, b in zip(targets, backgrounds)]
            identity_bytes = b''.join(hashlib.sha256(json.dumps([t, b, f, s], separators=(',', ':')).encode()).hexdigest().encode('ascii') for t, b, f, s in zip(targets, backgrounds, families, species))
            assert hashlib.sha256(identity_bytes).hexdigest() == recipe['ordered_identity_sha256']
            path = root / entry['path']
            assert sha(path) == entry['sha256']
            with np.load(path) as saved:
                assert set(saved.files) == {'matrix', 'row_identity', 'background', 'family', 'pattern_rows'}
                assert saved['matrix'].shape == (len(frame), recipe['polynomial_degree'] + 4)
                assert np.isfinite(saved['matrix']).all()
                assert hashlib.sha256(saved['matrix'].tobytes()).hexdigest() == recipe['values_sha256']
                assert saved['row_identity'].tobytes() == identity_bytes
                np.testing.assert_array_equal(saved['background'], codes(backgrounds))
                np.testing.assert_array_equal(saved['family'], codes(families))
                np.testing.assert_array_equal(saved['pattern_rows'], [pattern_rows[s] for s in species])
            checked.add(identifier)
            occurrences += len(frame)
        print('Verified flagged input partition', index, 'inputs', len(checked), '/449', flush=True)
    flags = pd.read_csv(root / 'flagged_tree_fits.tsv', sep='\t')
    original_flags = pd.read_csv('results/model_validation/polynomial-optimization-flag-census-20260928-v1/flagged_fits.tsv', sep='\t')
    pd.testing.assert_frame_equal(flags, original_flags)
    assert set(flags.fit_input_id) == checked and len(flags) == r['tree_fits'] == 601
    assert len(checked) == r['inputs'] == 449 and occurrences == r['record_occurrences'] == 1081570
    for path, digest in bindings.items():
        assert sha(path) == digest, path
    assert sha(root / 'receipt.json') == rh and sha(lp) == lh
    result = dict(status='passed_full_flagged_polynomial_input_readback', source_receipt_sha256=rh, checker_sha256=sha(__file__), producer_terminal_state=state, inputs=len(checked), tree_fits=len(flags), record_occurrences=occurrences, scope='All original numerical and identity recipe hashes matched; all memberships, background/family codes and species-factor rows independently reconstructed. All flagged fit links retained. No production refinement or flag resolution.')
    with Path('metadata/flagged_polynomial_inputs_completed_readback_20260928.json').open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
