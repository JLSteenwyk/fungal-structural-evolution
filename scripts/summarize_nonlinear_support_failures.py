"""Summarize all certificate failures after full independent matrix readback."""
import argparse
from collections import Counter
import json
import math
from pathlib import Path
import subprocess
import time
import psutil
from screen_duplication_domain_alignment_coverage import sha


def reasons(r):
    if r['solver_status'] != 0:
        return ['solver_failure']
    flags = []
    w = r['support_weights']
    if any(v < 0 for v in w): flags.append('negative_weight')
    if abs(math.fsum(w)-1) > 1e-8: flags.append('weight_sum')
    if abs(r['primal_distance']-r['solver_objective']) > 1e-8: flags.append('primal_objective')
    if r['direction_l1_norm'] > 1+1e-8: flags.append('direction_norm')
    if r['separating_lower_bound'] > r['primal_distance']+1e-8: flags.append('inconsistent_bounds')
    assert bool(flags) == (not r['certificate_valid'])
    if r['classification'] == 'unresolved_near_boundary': flags.append('near_boundary')
    return flags


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--script-sha256',required=True)
    p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--launch',type=Path,required=True)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();assert sha(__file__)==a.script_sha256
    lh=sha(a.launch); launch=json.loads(a.launch.read_text())
    ph=sha(a.plan);plan=json.loads(a.plan.read_text());assert launch['plan_sha256']==ph
    for path,h in plan['pins'].items():assert sha(path)==h,path
    while True:
        try:
            process=psutil.Process(launch['pid'])
            if process.create_time()!=launch['created'] or process.status()==psutil.STATUS_ZOMBIE: break
            assert process.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess: break
        time.sleep(30)
    state=subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True)
    assert dict(l.split('=',1) for l in state.splitlines())==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
    assert sha(a.launch)==lh
    receipt_path=a.source/'receipt.json'; rh=sha(receipt_path)
    receipt=json.loads(receipt_path.read_text())
    assert receipt['status']=='passed_full_nonlinear_joint_support_certificate_readback'
    for name,h in receipt['artifacts'].items(): assert sha(a.source/name)==h,name
    rows=[json.loads(l) for l in (a.source/'joint_support.jsonl').read_text().splitlines()]
    assert len(rows)==len({r['fit_input_id'] for r in rows})==receipt['unique_inputs']==57616
    classification=Counter(r['classification'] for r in rows)
    assert dict(classification)==receipt['classification_counts']
    counts=Counter(); details=[];degrees=Counter()
    for r in rows:
        flags=reasons(r); counts.update(flags)
        degrees[str(r['polynomial_degree'])+':'+r['classification']]+=1
        if not flags: continue
        w=r.get('support_weights',[])
        details.append(dict(polynomial_degree=r['polynomial_degree'],fit_input_id=r['fit_input_id'],partition_index=r['partition_index'],
            classification=r['classification'],reasons=flags,negative_weight_count=sum(v<0 for v in w),
            negative_weight_mass=math.fsum(-v for v in w if v<0),
            weight_sum=math.fsum(w) if w else None,primal_distance=r.get('primal_distance')))
    a.output.mkdir(parents=True,exist_ok=False)
    dest=a.output/'unresolved_inputs.jsonl'
    dest.write_text(''.join(json.dumps(r,allow_nan=False)+'\n' for r in details))
    assert [json.loads(l) for l in dest.read_text().splitlines()]==details
    assert sha(receipt_path)==rh
    for name,h in receipt['artifacts'].items(): assert sha(a.source/name)==h,name
    assert sha(__file__)==a.script_sha256
    assert sha(a.plan)==ph
    for path,h in plan['pins'].items():assert sha(path)==h,path
    summary=dict(status='complete_full_nonlinear_joint_support_failure_summary',plan_sha256=ph,degree_classification_counts=dict(degrees),unique_inputs=len(rows),
        unresolved_inputs=len(details),classification_counts=dict(classification),failure_reason_counts=dict(counts),
        maximum_negative_mass=max([r['negative_weight_mass'] for r in details],default=0),
        source_receipt_sha256=rh,launch_sha256=lh,script_sha256=sha(__file__),
        artifacts={dest.name:sha(dest)},scope='All independently checked certificates retained. Failure flags describe numerical certificates, not biological effects. No weights clipped, certificates repaired, inputs dropped, or unresolved results promoted.')
    (a.output/'receipt.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary),flush=True)

if __name__=='__main__': main()
