#!/usr/bin/env python3
"""Wait for completed native execution, then produce separate small-family deltas."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
from assess_small_family_output_exposure import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',required=True,type=Path)
    a=p.parse_args();plan=json.loads(a.plan.read_text());pins={str(a.plan):sha(a.plan),**plan['pins']}
    def verify():
        for path,digest in pins.items():
            if sha(path)!=digest:raise ValueError('Changed pin: '+path)
    verify();out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    started=time.time();completed=[]
    def state(status,**details):
        temp=out/'state.tmp'
        temp.write_text(json.dumps(dict(status=status,completed=completed,elapsed_seconds=time.time()-started,**details),indent=2)+'\n')
        temp.replace(out/'state.json')
    try:
        state('waiting_for_exact_native_controller')
        proc=Path('/proc')/str(plan['predecessor_pid'])/'stat'
        while proc.exists():
            try:fields=proc.read_text().rsplit(')',1)[1].split()
            except FileNotFoundError:break
            if fields[19]!=str(plan['predecessor_start_ticks']) or fields[0]=='Z':break
            time.sleep(20)
        verify()
        execution=json.loads(Path(plan['execution_plan']).read_text())
        base=Path(execution['output']);completion=base/'state.json'
        native=json.loads(completion.read_text())
        if native['status']!='both_native_executions_complete_pending_full_output_readback':
            raise ValueError('Both native executions must have completed successfully')
        for guide in ['profile','mafft']:
            if shutil.disk_usage(out).free<plan['resources']['minimum_free_disk_gib']*2**30:
                raise ValueError('Insufficient free disk')
            verify();state('scanning_completed_native_orthologs',guide=guide)
            target=out/guide
            command=[sys.executable,'scripts/supplement_small_family_orthologs.py',
                '--source',str(base/guide/'Source/WorkingDirectory'),
                '--results',str(base/guide/('Results_expanded_'+guide)),
                '--completion',str(completion),'--execution-plan',plan['execution_plan'],
                '--output',str(target)]
            with (out/(guide+'.log')).open('w') as log:
                subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
            r=json.loads((target/'receipt.json').read_text())
            if r['status']!='complete_separate_small_family_ortholog_supplement':raise ValueError('Incomplete supplement')
            completed.append(dict(guide=guide,receipt_sha256=sha(target/'receipt.json'),
                                  supplemental_directed_pairs=r['supplemental_directed_pairs']))
        verify();state('complete_separate_supplements_pending_full_reconciliation_readback')
    except Exception as exc:
        state('failed_requires_review',error=repr(exc));raise


if __name__=='__main__':main()
