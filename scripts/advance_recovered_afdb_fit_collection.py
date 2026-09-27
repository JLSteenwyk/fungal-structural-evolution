#!/usr/bin/env python3
"""Wait for the original fit audit, combine selected sources, and read back every table field."""
import argparse,json,subprocess,sys,time
from pathlib import Path
import psutil
from run_ortholog_pair_guide_comparison import sha

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
plan=json.loads(a.plan.read_text());ph=sha(a.plan)
def verify():
    if sha(a.plan)!=ph:raise ValueError('Changed collection plan')
    for path,h in plan['pins'].items():
        if sha(path)!=h:raise ValueError('Changed source: '+path)
verify();dep=plan['producer']
while psutil.pid_exists(dep['pid']):
    try:
        proc=psutil.Process(dep['pid'])
        if abs(proc.create_time()-dep['created'])>.01 or proc.status()==psutil.STATUS_ZOMBIE:break
        if proc.cmdline()!=dep['cmdline']:raise ValueError('Changed recovery identity')
    except psutil.NoSuchProcess:break
    time.sleep(30)
verify();recovery=json.loads(Path(plan['recovery_receipt']).read_text())
if recovery['status']!='complete_paired_fit_checkpoint_recovery' or recovery['plan_sha256']!=sha(plan['recovery_plan']):raise ValueError('Incomplete fit recovery')
control=Path(plan['control']);control.mkdir(parents=True,exist_ok=False)
args=plan['combine_arguments'];command=[sys.executable,'scripts/combine_audited_paired_fit_sources.py']
for key,value in args.items():command.extend(['--'+key,value])
commands=[('combine',command),('readback',[sys.executable,'scripts/readback_combined_paired_fit_sources.py','--collection',args['output'],'--inventory',args['inventory'],'--output',str(control/'table_readback.json')])]
for stage,command in commands:
    verify();(control/'state.json').write_text(json.dumps(dict(stage=stage,command=command))+'\n')
    with (control/(stage+'.log')).open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
verify();result=json.loads((control/'table_readback.json').read_text())
if result['status']!='passed_complete_combined_paired_fit_source_readback' or result['markers']!=125:raise ValueError('Incomplete source collection readback')
receipt=dict(status='complete_recovered_afdb_fit_collection_and_readback',plan_sha256=ph,markers=125,collection_receipt_sha256=sha(Path(args['output'])/'receipt.json'),readback_sha256=sha(control/'table_readback.json'),scope='Source-preserving audited point-estimate collection and independent full-table readback. Branch uncertainty and biological acceleration tests remain separate.')
(control/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
