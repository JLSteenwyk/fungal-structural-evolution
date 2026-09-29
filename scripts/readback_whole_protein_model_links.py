#!/usr/bin/env python3
"""Validate every whole-protein comparison link against serialized observations."""
import hashlib
import json
import subprocess
import time
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import psutil
from screen_duplication_alignment_reuse import sha

GROUP = ['guide', 'policy', 'scenario_id', 'mask', 'cohort', 'screen', 'target_order', 'background_order', 'outcome']
VARIANTS = ['gene_distance_linear', 'gene_distance_quadratic', 'gene_distance_cubic', 'positive_log_gene_distance_linear', 'alignment_identity_linear']


def check_relation(actual, left, right):
    same_rows = left['records'] == right['records'] and left['ordered_identity_sha256'] == right['ordered_identity_sha256']
    if not same_rows:
        assert actual == 'different_observations_no_direct_comparison'
        return
    for name in ['outcome', 'covariance_receipt_sha256', 'response_sha256']:
        assert left[name] == right[name], name
    a, b = left['design_column_hashes'], right['design_column_hashes']
    common = a.keys() & b.keys()
    assert all(a[k] == b[k] for k in common)
    left_only, right_only = a.keys()-b.keys(), b.keys()-a.keys()
    expected = ('identical_named_design' if not left_only and not right_only else
                'left_named_columns_nested_in_right' if not left_only else
                'right_named_columns_nested_in_left' if not right_only else
                'same_observations_no_named_column_nesting')
    assert actual == expected


def main():
    pp = Path('metadata/whole_protein_model_links_readback_plan_20260929.json')
    plan = json.loads(pp.read_text())
    bindings = {str(pp):sha(pp), **plan['pins']}
    launch = json.loads(Path(plan['launch']).read_text())
    while psutil.pid_exists(launch['pid']):
        try:
            p = psutil.Process(launch['pid'])
            if abs(p.create_time()-launch['created']) > .01 or p.status() == psutil.STATUS_ZOMBIE:
                break
            assert p.cmdline() == launch['cmdline']
        except psutil.NoSuchProcess:
            break
        time.sleep(30)
    state = dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert state == dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
    producer = json.loads(Path(launch['plan']).read_text())
    root = Path(producer['output'])
    receipt = json.loads((root/'receipt.json').read_text())
    assert receipt['status'] == 'complete_whole_protein_comparison_input_links_pending_readback'
    bindings.update(receipt['source_hashes'])
    bindings.update({str(root/n):h for n,h in receipt['artifacts'].items()})
    bindings[str(root/'receipt.json')] = sha(root/'receipt.json')
    def verify():
        for path,h in bindings.items():
            assert sha(path) == h,path
    verify()
    signatures = {}
    for line in (root/'input_column_hashes.jsonl').open():
        row = json.loads(line)
        key = row.pop('fit_input_id')
        assert key not in signatures
        signatures[key] = row
    inputs = Path(producer['inputs'])
    seen = set()
    for line in (inputs/'input_manifest.jsonl').open():
        item = json.loads(line)
        key = item['fit_input_id']
        assert key not in seen
        seen.add(key)
        spec = item['recipe']['specification']
        with np.load(inputs/item['path'],allow_pickle=False) as arrays:
            numeric = arrays['matrix']
            assert numeric.shape == (spec['records'],len(spec['columns']))
            assert hashlib.sha256(numeric.tobytes()).hexdigest() == spec['values_sha256']
            assert hashlib.sha256(arrays['row_identity'].tobytes()).hexdigest() == spec['ordered_identity_sha256']
            hashes = {name:hashlib.sha256(numeric[:,i].astype('<f8').tobytes()).hexdigest() for i,name in enumerate(spec['columns'])}
        expected = {k:spec[k] for k in ['records','outcome','ordered_identity_sha256','covariance_receipt_sha256']}
        expected['response_sha256'] = hashes.pop(spec['columns'][0])
        expected['design_column_hashes'] = hashes
        assert expected == signatures[key]
        if len(seen)%5000 == 0:
            print('Checked model-link column signatures',len(seen),flush=True)
    assert seen == set(signatures) and len(seen) == receipt['inputs']
    groups = defaultdict(dict)
    for line in (Path(producer['inventory'])/'setting_input_map.jsonl').open():
        row = json.loads(line)
        group = tuple(row[k] for k in GROUP)
        assert row['variant'] not in groups[group]
        groups[group][row['variant']] = row['fit_input_id']
    assert len(groups) == receipt['setting_outcome_groups'] == 82944
    for models in groups.values():
        assert set(models) == set(VARIANTS)
    visited = set()
    counts = Counter()
    for line in (root/'comparison_input_map.jsonl').open():
        row = json.loads(line)
        group = tuple(row[k] for k in GROUP)
        left, right = row['left_variant'],row['right_variant']
        assert VARIANTS.index(left) < VARIANTS.index(right)
        key = (*group,left,right)
        assert key not in visited
        visited.add(key)
        assert row['left_input'] == groups[group][left] and row['right_input'] == groups[group][right]
        check_relation(row['relation'],signatures[row['left_input']],signatures[row['right_input']])
        counts[row['relation']] += 1
    assert len(visited) == len(groups)*10 == receipt['comparisons'] == 829440
    assert dict(counts) == receipt['relation_counts']
    verify()
    proof = dict(status='passed_full_whole_protein_comparison_input_link_readback',inputs=len(seen),comparisons=len(visited),relation_counts=dict(counts),source_receipt_sha256=sha(root/'receipt.json'),checker_sha256=sha(__file__),producer_terminal_state=state,scope='Every serialized column signature and within-setting pair independently checked, including response/common-column equality and observation exclusions. No fitted model preference, general column-space equivalence, calibrated test or biological inference.')
    with Path(plan['proof']).open('x') as handle:
        json.dump(proof,handle,indent=2)
        handle.write('\n')
    print(json.dumps(proof),flush=True)


if __name__=='__main__':
    main()
