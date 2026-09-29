#!/usr/bin/env python3
"""Restartable complete whole-protein ordinary-ML grid; requires audited inputs."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import fcntl
import hashlib
import json
import multiprocessing
from pathlib import Path
import subprocess
import traceback

import numpy as np
from threadpoolctl import threadpool_limits
from fit_whole_protein_ml import fit_input
from screen_duplication_alignment_reuse import sha

FACTORS = {}


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def initialize(root):
    global FACTORS
    threadpool_limits(limits=1)
    FACTORS = {}
    for path in sorted(Path(root).glob('*.npz')):
        with np.load(path, allow_pickle=False) as saved:
            FACTORS[path.stem] = saved['factor']
    assert len(FACTORS) == 5


def fit_case(task):
    item, input_root, output_root, plan_hash = task
    identifier = item['fit_input_id']
    path = Path(input_root) / item['path']
    assert sha(path) == item['sha256']
    with np.load(path, allow_pickle=False) as saved:
        arrays = {name: saved[name] for name in saved.files}
    spec = item['recipe']['specification']
    assert digest(spec) == identifier
    assert hashlib.sha256(arrays['matrix'].tobytes()).hexdigest() == spec['values_sha256']
    assert hashlib.sha256(arrays['row_identity'].tobytes()).hexdigest() == spec['ordered_identity_sha256']
    assert len(arrays['matrix']) == spec['records']
    directory = Path(output_root) / 'fits' / identifier[:2] / identifier
    directory.mkdir(parents=True, exist_ok=True)
    outcomes = []
    for tree, full_factor in FACTORS.items():
        result_path = directory / (tree + '.json')
        identity = dict(plan_sha256=plan_hash, fit_input_id=identifier, tree=tree,
                        input_sha256=item['sha256'], specification=spec)
        if result_path.exists():
            saved = json.loads(result_path.read_text())
            assert all(saved[k] == v for k, v in identity.items())
            assert saved['payload_sha256'] == digest(saved['payload'])
        else:
            try:
                rows = arrays['pattern_rows']
                assert rows.dtype.kind in 'iu' and np.all(rows >= 0) and np.all(rows < len(full_factor))
                payload = fit_input(arrays['background'], arrays['family'], full_factor[rows], arrays['matrix'], spec['columns'])
            except Exception as error:
                payload = dict(status='fit_error_requires_review', error_type=type(error).__name__, error=str(error), traceback=traceback.format_exc())
            saved = dict(identity, payload=payload, payload_sha256=digest(payload))
            temporary = result_path.with_suffix('.json.tmp')
            temporary.write_text(json.dumps(saved, indent=2, allow_nan=False) + '\n')
            temporary.replace(result_path)
        outcomes.append(dict(fit_input_id=identifier, tree=tree, status=saved['payload']['status'], path=str(result_path), sha256=sha(result_path)))
    return outcomes


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    ph = sha(args.plan)
    bindings = {str(args.plan): ph, **plan['pins']}
    # Resource assessment is a pinned, concrete record, written after workload review.
    assessment = json.loads(Path(plan['resource_assessment']).read_text())
    assert plan['resource_assessment'] in plan['pins']
    assert assessment['status'] == 'whole_protein_fit_resources_assessed'
    assert assessment['workers'] == plan['workers'] and plan['workers'] > 0
    assert assessment['memory_gib'] > 0 and assessment['output_gib'] > 0
    assert assessment['runtime_scenarios'] and assessment['uncovered_workload_handling']
    root = Path(plan['inputs'])
    receipt = json.loads((root / 'receipt.json').read_text())
    proof = json.loads(Path(plan['input_audit']).read_text())
    assert proof['status'] == 'passed_full_whole_protein_materialized_input_readback'
    assert proof['source_receipt_sha256'] == sha(root / 'receipt.json')
    assert receipt['inputs'] == 75070 and receipt['tree_fits'] == 375350
    workload = json.loads(Path(plan['workload']).read_text())
    checked = json.loads(Path(plan['workload_audit']).read_text())
    assert checked['status'] == 'passed_whole_protein_workload_arithmetic_readback'
    assert checked['source_result_sha256'] == sha(plan['workload'])
    assert workload['unique_inputs'] == receipt['inputs'] and workload['tree_fits'] == receipt['tree_fits']
    for name in ['input_audit', 'workload', 'workload_audit', 'factor_audit']:
        assert plan[name] in plan['pins']
    audit_launch = json.loads(Path(plan['input_audit_launch']).read_text())
    state = dict(line.split('=', 1) for line in subprocess.check_output(
        ['systemctl', '--user', 'show', audit_launch['unit'], '-p', 'ActiveState', '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
    assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0')
    bindings.update(receipt['source_hashes'])
    bindings[str(root / 'receipt.json')] = sha(root / 'receipt.json')
    bindings.update({str(root / name): h for name, h in receipt['artifacts'].items()})
    factor_root = Path(plan['factors'])
    factor_receipt = json.loads((factor_root / 'receipt.json').read_text())
    factor_proof = json.loads(Path(plan['factor_audit']).read_text())
    assert factor_proof['source_receipt_sha256'] == sha(factor_root / 'receipt.json')
    assert factor_proof['status'].startswith('passed_full_')
    bindings[str(factor_root / 'receipt.json')] = sha(factor_root / 'receipt.json')
    bindings.update({str(factor_root / name): h for name, h in factor_receipt['artifacts'].items()})
    def verify():
        for path, expected in bindings.items():
            assert sha(path) == expected, path
    verify()
    output = Path(plan['output'])
    output.mkdir(parents=True, exist_ok=True)
    lock = (output / 'run.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    for filename, content in [('run_plan.json', args.plan.read_text()), ('source_bindings.json', json.dumps(bindings, sort_keys=True, indent=2)+'\n')]:
        path = output / filename
        if path.exists():
            assert path.read_text() == content
        else:
            path.write_text(content)
    pending, submitted = set(), set()
    counts = Counter()
    completed = 0
    with (output / 'fit_manifest.jsonl').open('w') as manifest, ProcessPoolExecutor(
            max_workers=plan['workers'], mp_context=multiprocessing.get_context('spawn'),
            initializer=initialize, initargs=(str(factor_root),)) as pool:
        def collect():
            nonlocal pending, completed
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                for row in future.result():
                    counts[row['status']] += 1
                    completed += 1
                    manifest.write(json.dumps(row, sort_keys=True)+'\n')
                manifest.flush()
                print('Completed whole-protein ML dispositions', completed, '/', receipt['tree_fits'], dict(counts), flush=True)
            errors = counts['fit_error_requires_review']
            if errors >= 10 and errors / completed > .1:
                raise RuntimeError('Numerical errors exceed ten percent after ten errors; outputs preserved')
        for line in (root / 'input_manifest.jsonl').open():
            item = json.loads(line)
            assert item['fit_input_id'] not in submitted
            submitted.add(item['fit_input_id'])
            pending.add(pool.submit(fit_case, (item, str(root), str(output), ph)))
            if len(pending) >= 2 * plan['workers']:
                collect()
        while pending:
            collect()
    assert len(submitted) == receipt['inputs'] and completed == receipt['tree_fits']
    verify()
    result = dict(status='complete_whole_protein_ml_dispositions_pending_full_audit', plan_sha256=ph,
                  unique_inputs=len(submitted), tree_fit_dispositions=completed, status_counts=dict(counts),
                  artifacts={name:sha(output/name) for name in ['run_plan.json', 'source_bindings.json', 'fit_manifest.jsonl']},
                  scope='All whole-protein exact inputs across five trees attempted; errors and review flags retained. Full output audit, optimization review, matched-observation model comparisons, uncertainty calibration and biological interpretation remain required.')
    (output / 'receipt.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
