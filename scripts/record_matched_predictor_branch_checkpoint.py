#!/usr/bin/env python3
"""Check both original matched-predictor workflows without native restarts."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from matched_predictor_branch_inputs import verify
from record_project_runtime_checkpoint_v4 import fingerprint,journal_terminal,live_record


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    handles=[];closures=[];plans=[]
    for phase in ['inputs','fits']:
        path=Path('metadata/matched_predictor_branch_'+phase+'_plan_20261003_v1.json')
        plan=json.loads(path.read_text());verify(plan['pins']);plans.append(dict(plan=str(path),sha256=sha(path),pins_checked=len(plan['pins'])))
        inventory=json.loads(Path(plan['launch_inventory']).read_text());assert inventory['source_plan_sha256']==sha(path)
        for index,lp in enumerate(inventory['launches']):
            record=json.loads(Path(lp).read_text());record['launch']=lp
            assert sha(record['plan'])==record['plan_sha256']
            process=fingerprint(record)
            if process is None:observed=journal_terminal(record)
            else:
                observed=live_record(record,process)
                group=next(x[3:] for x in (Path('/proc')/str(record['pid'])/'cgroup').read_text().splitlines() if x.startswith('0::'))
                cg=Path('/sys/fs/cgroup')/group.lstrip('/')
                limits={key:(cg/key).read_text().strip() for key in ['cpu.max','memory.max','memory.swap.max']}
                assert limits==record['actual_cgroup_limits'];observed['current_cgroup_limits']=limits
            handles.append(observed)
        cp=Path(plan['completion'])
        if cp.exists():
            c=json.loads(cp.read_text());assert sha(c['full_hash_archive'])==c['full_hash_archive_sha256']
            assert c['exact_process_journals_checked']==2 and c['scientific_eligibility'] is False
            closures.append(dict(phase=phase,completion=str(cp),completion_sha256=sha(cp),summary=c))
    result=dict(status='verified_original_matched_predictor_branch_execution_checkpoint',checked_utc=datetime.now(timezone.utc).isoformat(),
        source_plans=plans,original_handles=handles,completed_closures=closures,native_restarts=0,gpu=False,scientific_eligibility=False,
        scope='Exact original six controller identities/journals, frozen pins and live limits if present; compact full archive hashes checked for completed phases. No new native execution or full native replay/optimization, model/framework acceptance or project completion.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],original_handles=len(handles),completed_phases=[c['phase'] for c in closures],checked_utc=result['checked_utc'])),flush=True)


if __name__=='__main__':main()
