#!/usr/bin/env python3
"""Start the frozen ecology queue after a successful predecessor frees its GPU."""
import argparse,fcntl,json,os,shutil,subprocess,time
from pathlib import Path
from audit_busco_gene_copies import ROOT,sha
from assess_pae_sensitivity import checked_receipt


def process_identity(pid):
    try:
        stat=(Path('/proc')/str(pid)/'stat').read_text()
    except FileNotFoundError:return None
    fields=stat[stat.rfind(')')+2:].split()
    return None if fields[0]=='Z' else fields[19]


def completed_chunk(row,config_sha):
    return (row.get('status')=='production_chunk_finished' and row.get('config_sha256')==config_sha
            and row.get('interrupted') is False and row.get('remaining_eligible')==0 and row.get('oom_deferred')==0)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,required=True);a=p.parse_args();c=json.loads(a.config.read_text())
    control=ROOT/c['control_output'];control.mkdir(parents=True,exist_ok=True)
    lock=(control/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    if (control/'launch.json').exists():raise FileExistsError('A launch already occurred; inspect its outcome before recovery')
    if sha(Path(__file__))!=c['controller_sha256']:raise ValueError('Changed controller')
    config_sha=sha(a.config)
    print('Waiting for the pinned predecessor process to finish successfully',flush=True)
    while process_identity(c['predecessor_pid'])==c['predecessor_start_ticks']:time.sleep(20)
    pred=ROOT/c['predecessor_output'];pc=pred/'config.json';chunk=pred/'last_chunk.json'
    if sha(pc)!=c['predecessor_config_sha256'] or not completed_chunk(json.loads(chunk.read_text()),sha(pc)):raise ValueError('Predecessor did not complete cleanly; inspect before scheduling')
    inputs=ROOT/c['inputs'];checked_receipt(inputs)
    for name,digest in c['pinned_files'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Changed pinned dependency: '+name)
    r=json.loads((inputs/'receipt.json').read_text())
    if r['prediction_candidates']!=c['limit']:raise ValueError('Queue size differs')
    if (ROOT/c['prediction_output']).exists():raise FileExistsError('Use a fresh prediction output')
    if shutil.disk_usage(ROOT).free<c['minimum_free_disk_bytes']:raise ValueError('Insufficient planned output headroom')
    while True:
        apps=subprocess.run(['nvidia-smi','--query-compute-apps=gpu_uuid,pid','--format=csv,noheader'],capture_output=True,text=True,check=True).stdout
        occupied={line.split(',')[0].strip() for line in apps.splitlines() if line.strip()}
        if c['gpu_uuid'] not in occupied:break
        time.sleep(20)
    # Recheck pins immediately before launch, after any device wait.
    for name,digest in c['pinned_files'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Dependency changed during device wait')
    command=[c['python'],str(ROOT/'scripts/run_marker_predictions.py'),'--inputs',str(inputs),'--checkpoint',str(ROOT/c['checkpoint']),'--output',str(ROOT/c['prediction_output']),'--max-length','512','--limit',str(c['limit'])]
    launch={'status':'launching_ecology_predictions','controller_config_sha256':config_sha,'predecessor_chunk_sha256':sha(chunk),'gpu_process_observation':apps,'command':command,'gpu_uuid':c['gpu_uuid'],'unix_time':time.time()}
    (control/'launch.json').write_text(json.dumps(launch,indent=2)+'\n');print('Launching ecology predictions',flush=True)
    env=os.environ.copy();env['CUDA_VISIBLE_DEVICES']=c['gpu_uuid']
    with (control/'prediction.log').open('w') as f:result=subprocess.run(command,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
    out=ROOT/c['prediction_output'];final=out/'last_chunk.json'
    success=result.returncode==0 and final.exists() and completed_chunk(json.loads(final.read_text()),sha(out/'config.json'))
    receipt={'status':'complete_ecology_prediction_chunk' if success else 'ecology_prediction_requires_review','exit_code':result.returncode,'controller_config_sha256':config_sha,'launch_sha256':sha(control/'launch.json'),'prediction_chunk_sha256':sha(final) if final.exists() else None,'interpretation':'Controller completion verifies clean chunk termination, not an independent model/PAE artifact audit.'}
    (control/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
    if not success:raise SystemExit(1)


if __name__=='__main__':main()
