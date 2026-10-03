#!/usr/bin/env python3
"""Check whole-grid recovery selection contracts without launching native jobs."""
import argparse
import copy
from datetime import datetime, timezone
import json
from pathlib import Path
from audit_baliphy_memory_recovery_20261002_v2 import build_overlay, changed_configuration, CHECKED, launch_plan_path
from run_ortholog_pair_guide_comparison import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    source_path=Path('metadata/baliphy_independent_chain_plan_20260927.json')
    retry_plan_path=Path('metadata/baliphy_memory_recovery_plan_20261002.json')
    source=json.loads(source_path.read_text());retry_plan=json.loads(retry_plan_path.read_text())
    original_path=Path(source['output'])/'receipt.json';retry_path=Path(retry_plan['output'])/'receipt.json'
    jobs=json.loads(Path(source['jobs']).read_text());original=json.loads(original_path.read_text())['chains']
    retry=json.loads(retry_path.read_text())['recovery_chains'];failed=retry_plan['failed_chain_ids']
    rows,groups,summary=build_overlay(jobs,original,retry,failed)
    assert summary==dict(full_native_chains=1620,full_quartets=405,original_checked_chains=1617,
        original_failed_chains=3,recovery_attempts=3,recovered_integrity_checked_chains=1,
        unresolved_failed_chains=2,selected_checked_chains=1618,complete_quartets=403,unresolved_quartets=2)
    original_index={r['chain_id']:r for r in original};retry_index={r['chain_id']:r for r in retry}
    for row in rows:
        cid=row['chain']['chain_id'];assert row['original_disposition']==original_index[cid]
        assert row['scientific_eligibility'] is False
        if cid not in retry_index:assert row['selected_disposition']==original_index[cid] and row['recovery_disposition'] is None
        else:
            assert row['selected_disposition']['receipt']==retry_index[cid]['new_attempt']
            assert row['recovery_disposition']==retry_index[cid]
    successful=next(i for i,r in enumerate(retry) if r['status']==CHECKED)
    bad=next(i for i,r in enumerate(retry) if r['status']=='failed')
    mutations=[]
    def reject(label, mutate):
        args=copy.deepcopy([jobs,original,retry,failed]);mutate(*args)
        try:build_overlay(*args)
        except (AssertionError,KeyError):mutations.append(label)
        else:raise AssertionError('Invalid full-grid overlay accepted: '+label)
    reject('missing_original_chain',lambda j,o,r,f:o.pop())
    reject('duplicate_original_chain',lambda j,o,r,f:o.append(copy.deepcopy(o[0])))
    reject('missing_recovery_disposition',lambda j,o,r,f:r.pop())
    reject('duplicate_recovery_disposition',lambda j,o,r,f:r.append(copy.deepcopy(r[0])))
    reject('incompatible_original_attempt_lineage',lambda j,o,r,f:r[successful].update(original_failed_attempt='foreign/receipt.json'))
    reject('changed_original_failed_receipt_digest',lambda j,o,r,f:r[successful].update(original_failed_attempt_sha256='0'*64))
    reject('failed_native_attempt_promoted',lambda j,o,r,f:r[bad].update(status=CHECKED))
    reject('failed_attempt_given_checked_audit',lambda j,o,r,f:r[bad].update(sample_audit='foreign/audit.json'))
    reject('successful_attempt_missing_audit',lambda j,o,r,f:r[successful].pop('sample_audit'))
    reject('invented_partial_recovery_status',lambda j,o,r,f:r[successful].update(status='partial_samples_usable'))
    reject('changed_quartet_seed',lambda j,o,r,f:j[0]['chain'].update(seed=j[1]['chain']['seed']))
    reject('changed_native_seed',lambda j,o,r,f:j[0]['config']['command'].__setitem__(6,'1'))
    reject('changed_horizon',lambda j,o,r,f:j[0]['config']['command'].__setitem__(10,'500'))
    reject('changed_quartet_identity',lambda j,o,r,f:j[0]['config'].update(model_input_identity='foreign'))
    reject('missing_original_alignment_alias',lambda j,o,r,f:j[0]['chain']['original_configuration_ids'].pop())
    reject('changed_native_program_within_quartet',lambda j,o,r,f:j[0]['chain'].update(program='foreign/model.hs'))
    changed=[]
    for job in jobs:
        if job['chain']['chain_id'] not in failed:continue
        config=changed_configuration(job['config'],48*2**30)
        actual=copy.deepcopy(config);actual['command'][1]='--as=12884901888'
        assert actual==job['config'] and config['seed']==job['chain']['seed']
        changed.append(job['chain']['chain_id'])
    try:changed_configuration(jobs[0]['config'],96*2**30)
    except AssertionError:mutations.append('unapproved_address_space_configuration')
    else:raise AssertionError('Unexpected memory configuration accepted')
    launch_formats=[]
    for launch_path in [*retry_plan['original_launches'], 'metadata/baliphy_memory_recovery_launch_20261002.json']:
        launch=json.loads(Path(launch_path).read_text());bound=launch_plan_path(launch)
        assert sha(bound)==launch['plan_sha256']
        launch_formats.append(dict(launch=launch_path,command_plan=bound,explicit_plan_present='plan' in launch))
        broken=copy.deepcopy(launch);broken['plan']='foreign/plan.json'
        try:launch_plan_path(broken)
        except AssertionError:pass
        else:raise AssertionError('Conflicting command/explicit plan accepted')
    mutations.append('conflicting_explicit_command_plan')
    paths=[Path(__file__),Path('scripts/audit_baliphy_memory_recovery_20261002_v2.py'),
        source_path,retry_plan_path,Path(source['jobs']),original_path,retry_path]
    result=dict(status='passed_full_baliphy_recovery_overlay_contracts',checked_utc=datetime.now(timezone.utc).isoformat(),
        full_original_grid=summary,unchanged_original_checked_dispositions=1617,
        resource_only_configuration_changes=changed,malformed_overlays_rejected=mutations,
        original_launch_formats_checked=launch_formats,source_hashes={str(p):sha(p) for p in paths},
        scope='Full original1620-chain/405-quartet identity and whole-attempt overlay contract; no native job launch, output integrity parsing, diagnostic metric acceptance, or posterior qualification.')
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
