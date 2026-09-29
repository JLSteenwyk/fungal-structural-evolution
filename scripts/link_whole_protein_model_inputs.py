#!/usr/bin/env python3
"""Link model inputs only after exact response, row and shared-column checks."""
import hashlib
import itertools
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


def relation(a, b):
    if a['records'] != b['records'] or a['ordered_identity_sha256'] != b['ordered_identity_sha256']:
        return 'different_observations_no_direct_comparison'
    assert a['outcome'] == b['outcome'] and a['covariance_receipt_sha256'] == b['covariance_receipt_sha256']
    assert a['response_sha256'] == b['response_sha256']
    left, right = a['design_column_hashes'], b['design_column_hashes']
    for column in set(left) & set(right):
        assert left[column] == right[column], column
    if set(left) == set(right):
        return 'identical_named_design'
    if set(left) < set(right):
        return 'left_named_columns_nested_in_right'
    if set(right) < set(left):
        return 'right_named_columns_nested_in_left'
    return 'same_observations_no_named_column_nesting'


def main():
    pp = Path('metadata/whole_protein_model_links_plan_20260929.json')
    plan = json.loads(pp.read_text())
    bindings = {str(pp): sha(pp), **plan['pins']}
    launch = json.loads(Path(plan['audit_launch']).read_text())
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
    root = Path(plan['inputs'])
    receipt = json.loads((root/'receipt.json').read_text())
    proof = json.loads(Path(plan['audit_proof']).read_text())
    assert proof['status'] == 'passed_full_whole_protein_materialized_input_readback'
    assert proof['source_receipt_sha256'] == sha(root/'receipt.json')
    bindings.update(receipt['source_hashes'])
    bindings.update({str(root/n):h for n,h in receipt['artifacts'].items()})
    bindings[str(root/'receipt.json')] = sha(root/'receipt.json')
    bindings[plan['audit_proof']] = sha(plan['audit_proof'])
    def verify():
        for path, h in bindings.items():
            assert sha(path) == h, path
    verify()
    out = Path(plan['output'])
    out.mkdir(exist_ok=False)
    signatures = {}
    with (out/'input_column_hashes.jsonl').open('w') as handle:
        for line in (root/'input_manifest.jsonl').open():
            item = json.loads(line)
            key = item['fit_input_id']
            assert key not in signatures
            spec = item['recipe']['specification']
            with np.load(root/item['path'], allow_pickle=False) as a:
                matrix = a['matrix']
                assert hashlib.sha256(matrix.tobytes()).hexdigest() == spec['values_sha256']
                assert hashlib.sha256(a['row_identity'].tobytes()).hexdigest() == spec['ordered_identity_sha256']
                hashes = [hashlib.sha256(np.ascontiguousarray(matrix[:,i], dtype='<f8').tobytes()).hexdigest() for i in range(matrix.shape[1])]
            sig = {k:spec[k] for k in ['records','outcome','ordered_identity_sha256','covariance_receipt_sha256']}
            sig.update(response_sha256=hashes[0], design_column_hashes=dict(zip(spec['columns'][1:],hashes[1:])))
            signatures[key] = sig
            handle.write(json.dumps(dict(fit_input_id=key,**sig))+'\n')
            if len(signatures)%5000 == 0:
                print('Hashed whole-protein model columns',len(signatures),flush=True)
    assert len(signatures) == receipt['inputs']
    inventory = Path(plan['inventory'])
    ir = json.loads((inventory/'receipt.json').read_text())
    assert sha(inventory/'setting_input_map.jsonl') == ir['artifacts']['setting_input_map.jsonl']
    groups = defaultdict(dict)
    for line in (inventory/'setting_input_map.jsonl').open():
        row = json.loads(line)
        key = tuple(row[k] for k in GROUP)
        assert row['variant'] not in groups[key] and row['fit_input_id'] in signatures
        groups[key][row['variant']] = row['fit_input_id']
    assert len(groups) == 82944
    counts = Counter()
    with (out/'comparison_input_map.jsonl').open('w') as handle:
        for key, inputs in sorted(groups.items()):
            assert set(inputs) == set(VARIANTS)
            for left,right in itertools.combinations(VARIANTS,2):
                status = relation(signatures[inputs[left]],signatures[inputs[right]])
                counts[status] += 1
                row = dict(zip(GROUP,key))
                row.update(left_variant=left,right_variant=right,left_input=inputs[left],right_input=inputs[right],relation=status)
                handle.write(json.dumps(row,separators=(',',':'))+'\n')
    assert sum(counts.values()) == 829440
    verify()
    result = dict(status='complete_whole_protein_comparison_input_links_pending_readback',inputs=len(signatures),setting_outcome_groups=len(groups),comparisons=sum(counts.values()),relation_counts=dict(counts),source_hashes=bindings,script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in out.iterdir()},scope='All ten within-setting predictor pairs for both outcomes. Same-observation links require identical response and shared named columns. Named-column nesting is sufficient but not exhaustive column-space equivalence. Different observations excluded from direct likelihood comparison. No fitted preference, significance, calibration or biological effect; full readback pending.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'}),flush=True)


if __name__=='__main__':
    main()
