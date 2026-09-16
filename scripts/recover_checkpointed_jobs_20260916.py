#!/usr/bin/env python3
"""Recover explicitly reviewed restartable jobs after the September 16 reboot."""
import hashlib,json,os,subprocess,sys,time
from pathlib import Path
import psutil
base=Path(__file__).resolve().parents[1];os.chdir(base)
out=base/'results/recovery-20260916';launch=out/'launch.json'
if launch.exists():raise FileExistsError('Inspect existing recovery launch before another attempt')
old=json.loads((base/'metadata/project_pause_20260915_manifest.json').read_text())
by_pid={r['pid']:r for r in old['roots']}
for r in old['roots']:
    try:
        if psutil.Process(r['pid']).create_time()==r['created']:raise RuntimeError('An original process is still alive')
    except psutil.NoSuchProcess:pass
cool=json.loads((base/'metadata/project_resume_20260916_cooling_config.json').read_text())
os.sched_setaffinity(0,cool['cpus'])
apps=subprocess.check_output(['nvidia-smi','--query-compute-apps=pid','--format=csv,noheader'],text=True)
if apps.strip():raise RuntimeError('GPU occupied; inspect before recovery')
ids=[2006275,3150233,3765266,3392402,3799913,3799915,18695,287662]
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
record={'status':'recovering_reviewed_checkpointable_jobs','boot_id':Path('/proc/sys/kernel/random/boot_id').read_text().strip(),'started_at':time.time(),'jobs':[],'pending_manual_review':['supported PMSF checkpoint recovery','OrthoFinder assignment continuation','replacement completion controllers']}
launch.write_text(json.dumps(record,indent=2)+'\n')
for old_pid in ids:
    command=by_pid[old_pid]['command'];jobenv=dict(env)
    if old_pid==2006275:jobenv['CUDA_VISIBLE_DEVICES']='GPU-56e78ad3-b4d7-54f6-38fa-e95729009960'
    if command[0]=='python':command[0]=sys.executable
    log=out/(str(old_pid)+'.log')
    with log.open('w') as f:p=subprocess.Popen(command,stdout=f,stderr=subprocess.STDOUT,env=jobenv,cwd=base)
    proc=psutil.Process(p.pid)
    record['jobs'].append({'old_pid':old_pid,'pid':p.pid,'created':proc.create_time(),'command':command,'log':str(log.relative_to(base))})
    launch.write_text(json.dumps(record,indent=2)+'\n')
cool['roots']=record['jobs'];cool['gpu']=next({'pid':r['pid'],'created':r['created']} for r in record['jobs'] if r['old_pid']==2006275)
cool['purpose']='Recovered existing jobs after reboot; same 24-core and 30-second run/rest limits with saved remaining cooling duration'
cp=base/'metadata/recovery_20260916_cooling_config.json';cp.write_text(json.dumps(cool,indent=2)+'\n')
record['status']='launched_reviewed_jobs_with_cooling';record['cooling_config_sha256']=hashlib.sha256(cp.read_bytes()).hexdigest();launch.write_text(json.dumps(record,indent=2)+'\n')
os.execv(sys.executable,[sys.executable,str(base/'scripts/temporary_run_resource_limit.py'),'--config',str(cp),'--state',str(out/'cooling_state.json')])
