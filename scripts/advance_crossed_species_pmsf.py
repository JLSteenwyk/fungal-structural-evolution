#!/usr/bin/env python3
"""Run remaining crossed PMSF combinations serially after a pinned predecessor."""
import argparse,fcntl,hashlib,json,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def identity(pid):
    try:
        p=Path('/proc')/str(pid);text=(p/'stat').read_text();fields=text[text.rfind(')')+2:].split();cmd=[x for x in (p/'cmdline').read_bytes().decode().split('\0') if x]
        return {'start_ticks':fields[19],'state':fields[0],'command':cmd}
    except FileNotFoundError:return None


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,required=True);a=p.parse_args()
    c=json.loads(a.config.read_text());ch=sha(a.config)
    def verify():
        if sha(a.config)!=ch:raise ValueError('Controller configuration changed')
        for name,digest in c['pinned_files'].items():
            if sha(ROOT/name)!=digest:raise ValueError('Pinned file changed: '+name)
    verify();out=ROOT/c['controller_output'];out.mkdir(parents=True,exist_ok=True)
    lock=(out/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    if (out/'receipt.json').exists():raise FileExistsError('Controller already complete')
    while True:
        current=identity(c['producer_pid'])
        if current is None or current['start_ticks']!=c['producer_start_ticks'] or current['state']=='Z':break
        if current['command']!=c['producer_command']:raise ValueError('Producer identity changed')
        print('waiting_for_first_pmsf',c['producer_pid'],flush=True);time.sleep(30)
    verify();pr=ROOT/c['producer_output']/'receipt.json';r=json.loads(pr.read_text())
    if r['status']!='complete_pmsf_execution_pending_full_audit' or r['returncode']!=0 or r['config_sha256']!=sha(ROOT/c['producer_output']/'config.json'):raise ValueError('First PMSF did not complete successfully')
    for name,digest in r['artifacts'].items():
        if sha(ROOT/c['producer_output']/name)!=digest:raise ValueError('Completed predecessor artifact changed')
    stages=[]
    for job in c['remaining_runs']:
        verify();checkpoint=out/(job['label']+'.json');receipt=ROOT/job['output']/'receipt.json'
        if checkpoint.exists():
            done=json.loads(checkpoint.read_text())
            if done['config_sha256']!=ch or done['receipt_sha256']!=sha(receipt):raise ValueError('Completed stage changed')
        else:
            if (ROOT/job['output']).exists():raise FileExistsError('Uncheckpointed stage requires review')
            command=[sys.executable,str(ROOT/'scripts/run_species_pmsf.py'),'--matrix',job['matrix'],'--guide-audit',job['guide_audit'],'--resources',c['resources'],'--output',job['output']]
            print('starting',job['label'],flush=True)
            with (out/(job['label']+'.log')).open('w') as log:subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
            result=json.loads(receipt.read_text())
            if result['status']!='complete_pmsf_execution_pending_full_audit' or result['returncode']!=0:raise ValueError('PMSF stage did not complete')
            done={'label':job['label'],'config_sha256':ch,'command':command,'receipt_path':str(receipt.relative_to(ROOT)),'receipt_sha256':sha(receipt)}
            checkpoint.write_text(json.dumps(done,indent=2)+'\n')
        stages.append(done)
    verify();result={'status':'complete_remaining_crossed_pmsf_execution_pending_full_audits','config_sha256':ch,'first_run_receipt_sha256':sha(pr),'remaining_runs':stages,'interpretation':'Serial execution of the remaining three baseline alignment/guide combinations. Every run still requires full model/profile/support audit; no final species-tree or adequacy claim.'}
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
