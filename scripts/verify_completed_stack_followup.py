#!/usr/bin/env python3
"""Fresh complete hash and original terminal-handle check of the closed followup."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
from collections import Counter

from ancestral_chain_attempt import sha
from baliphy_joint_sampler_gates_v3 import closed_stage
from record_project_runtime_checkpoint_v4 import fingerprint,journal_terminal
from reference_measurement_union_sources import bind,verify


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();assert not a.output.exists()
    plan,completion,pins=closed_stage('metadata/baliphy_stack_followup_plan_20261003_v1.json',
                                    'complete_verified_full_role_stack_followup_v1')
    assert completion['followup_roles']==24
    assert completion['selected_checked_sampler_attempts']==1596 and completion['selected_unsuccessful_sampler_attempts']==24
    verify(pins)
    rows=list((Path(plan['output'])/'attempts').rglob('receipt.json'));assert len(rows)==24
    counts=Counter(json.loads(q.read_text())['status'] for q in rows);assert counts=={'timeout':24}
    handles=[];inventory=Path(plan['launch_inventory']);bind(pins,inventory)
    for lp in json.loads(inventory.read_text())['launches']:
        launch=json.loads(Path(lp).read_text());bind(pins,lp)
        assert fingerprint(launch) is None and sha(launch['plan'])==launch['plan_sha256']
        observed=journal_terminal(launch)
        assert observed['status']=='verified_original_terminal_success_with_bound_completed_artifacts'
        handles.append(observed)
    bind(pins,Path(__file__));verify(pins)
    result=dict(status='fresh_verified_complete_original_stack_followup_hashes_and_terminal_handles',
        checked_utc=datetime.now(timezone.utc).isoformat(),original_completion_sha256=sha(plan['completion']),
        full_bindings_checked=len(pins),source_hashes=pins,original_terminal_handles=handles,
        original_successful_roles=1596,original_failed_roles_retained=24,followup_roles=24,
        followup_native_status_counts=dict(counts),selected_complete_quartets=399,selected_unresolved_quartets=6,
        scientific_eligibility=False,posterior_qualified=False,gpu=False,original_jobs_restarted=False,
        scope='Freshly rehashes the complete closed original followup archive and every bound file; '
              'checks all three exact original terminal controllers and24native timeout receipts. '
              'Inherited full output/resource/readback qualification is not reimplemented here. '
              'All24fresh attempts timed out; no additional successful sampler, adequate posterior, '
              'crash repair or completed biological aim is inferred. Original outputs unchanged.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','original_terminal_handles']},indent=2))


if __name__=='__main__':main()
