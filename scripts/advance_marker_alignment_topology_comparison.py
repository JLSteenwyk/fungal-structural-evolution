#!/usr/bin/env python3
"""Run the full marker topology comparison after complete MAFFT source audits."""
import argparse,hashlib,json,subprocess,sys,time
from pathlib import Path
import psutil


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args()
    plan=json.loads(a.plan.read_text());ph=sha(a.plan);root=Path(plan['controller_output']);root.mkdir(parents=True,exist_ok=False)
    def verify():
        assert sha(a.plan)==ph
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify();identity=plan['producer']
    while True:
        try:
            p=psutil.Process(identity['pid'])
            if p.create_time()!=identity['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==identity['cmdline']
        except psutil.NoSuchProcess:break
        print('waiting_for_full_mafft_marker_audit',identity['pid'],flush=True);time.sleep(30)
    verify();source=Path(plan['handoff_receipt']);r=json.loads(source.read_text())
    assert r['status']=='complete_full_mafft_marker_tree_audit_handoff' and r['markers']==125
    assert r['plan_sha256']==sha(plan['producer_plan'])
    source_hash=sha(source)
    for i,command in enumerate(plan['commands']):
        verify()
        with (root/f'stage_{i}.log').open('x') as f:subprocess.run([sys.executable,*command],stdout=f,stderr=subprocess.STDOUT,check=True)
    verify();audit=json.loads(Path(plan['readback']).read_text());comparison=Path(plan['comparison'])
    assert audit['status']=='passed_full_marker_topology_comparison_independent_pruning_readback' and audit['markers']==125
    assert audit['source_receipt_sha256']==sha(comparison/'receipt.json')
    assert sha(source)==source_hash
    result=dict(status='complete_full_marker_alignment_topology_comparison_and_readback',markers=125,changed_topologies=audit['changed_topologies'],comparison_receipt_sha256=sha(comparison/'receipt.json'),readback_sha256=sha(plan['readback']),source_handoff_sha256=source_hash,plan_sha256=ph,scope='All 125 markers compared and independently checked on shared taxa. Pruning is not refitting and differences do not identify a preferred alignment or biological cause.')
    (root/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main()
