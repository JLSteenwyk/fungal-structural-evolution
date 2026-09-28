"""Wait for terminal successful catalog readback, then compare full catalogs."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time
import psutil
from readback_whole_proteome_catalog import sha


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args()
    p=json.loads(a.plan.read_text());digest=sha(a.plan)
    def verify():
        assert sha(a.plan)==digest
        for name,h in p['pins'].items():assert sha(name)==h,name
    verify();launch=json.loads(Path(p['audit_launch']).read_text())
    while True:
        state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',launch['unit'],
            '-p','MainPID','-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
        if state['ActiveState'] in ['inactive','failed']:
            assert state['ActiveState']=='inactive' and state['Result']=='success' and state['ExecMainStatus']=='0',state
            break
        assert state['ActiveState']=='active' and int(state['MainPID'])==launch['pid']
        try:
            proc=psutil.Process(launch['pid'])
            assert proc.create_time()==launch['created'] and proc.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:
            time.sleep(1);continue
        print('Waiting for refreshed catalog readback',flush=True);time.sleep(30)
    verify()
    subprocess.run([sys.executable,'scripts/compare_whole_proteome_catalogs.py','--plan',p['comparison_plan']],check=True)
    verify()

if __name__=='__main__':main()
