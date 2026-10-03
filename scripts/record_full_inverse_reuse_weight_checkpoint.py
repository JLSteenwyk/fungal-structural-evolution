#!/usr/bin/env python3
"""Observe original control-catalog handles and checkpoint counts without replay."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from record_process_covariance_pipeline_checkpoint import observe
from reference_measurement_union_sources import verify


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    assert not args.output.exists()
    pp=Path('metadata/full_inverse_reuse_weight_plan_20261003_v1.json');plan=json.loads(pp.read_text())
    verify(plan['pins'])
    ip=Path('metadata/full_inverse_reuse_weight_launches_20261003_v1.json');inv=json.loads(ip.read_text())
    assert inv['source_plan']==str(pp) and inv['source_plan_sha256']==sha(pp)
    expected={'cpu.max':'200000 100000','memory.max':str(16*2**30),'memory.swap.max':'0'}
    handles=[observe(p,expected) for p in inv['launches']]
    root=Path(plan['output']);reports=list((root/'cohorts').glob('*.json'))
    gates={}
    for name,path in [('producer',root/'receipt.json'),('reader',root/'readback.json'),('completion',Path(plan['completion']))]:
        gates[name]=dict(path=str(path),present=path.exists())
        if path.exists():
            value=json.loads(path.read_text());assert value['scientific_eligibility'] is False
            gates[name].update(status=value['status'],sha256=sha(path),
                cohorts=value['cohorts'],case_row_occurrences=value['case_row_occurrences'],
                case_control_occurrences=value['case_control_occurrences'])
            if name=='completion':
                assert sha(value['full_hash_archive'])==value['full_hash_archive_sha256']
                gates[name].update(bound_source_hashes=value['bound_source_hashes'],
                    exact_process_journals_checked=value['exact_process_journals_checked'])
    result=dict(status='verified_original_full_inverse_reuse_control_runtime',
        checked_utc=datetime.now(timezone.utc).isoformat(),plan_sha256=sha(pp),original_handles=handles,
        frozen_pins_checked=len(plan['pins']),cohort_checkpoint_files=len(reports),expected=plan['expected'],
        completion_gates=gates,weighted_fitting_launched=False,gpu=False,all_eight_aims_incomplete=True,
        scope='Original PID/create/command or invocation journals, actual live cgroup caps and available checkpoint counts. Full stage closure archive digest if present; no replay of full saved controls or broader upstream artifacts. No calibrated precision, accepted fits, posterior adequacy, restart or final peak claim.')
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='original_handles'},indent=2))


if __name__=='__main__':main()
