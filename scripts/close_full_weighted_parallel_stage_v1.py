#!/usr/bin/env python3
"""Close full parallel weighted qualification, including all worker checkpoints."""
import argparse
import json
import subprocess
import sys
from pathlib import Path
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    plan=json.loads(a.plan.read_text());bindings=dict(plan['pins']);bind(bindings,a.plan)
    rp,ap=map(Path,[plan['producer_receipt'],plan['independent_readback']])
    producer,reader=[json.loads(p.read_text()) for p in [rp,ap]]
    assert producer['status']==plan['producer_status'] and reader['status']==plan['reader_status']
    assert reader['producer_receipt_sha256']==sha(rp)
    assert producer['plan_sha256']==reader['plan_sha256']==sha(plan['source_plan'])
    summary={k:reader[k] for k in plan['summary_fields']}
    assert all(producer[k]==v for k,v in summary.items())
    for path,r in [(rp,producer),(ap,reader)]:
        assert r['scientific_eligibility'] is False;bind(bindings,path)
        for q,d in r['source_hashes'].items():bind(bindings,q,d)
    for name,d in producer['artifacts'].items():bind(bindings,rp.parent/name,d)
    checkpoints=reader['reader_checkpoint_artifacts']
    assert len(checkpoints)==producer['cohorts']==4340
    for name,d in checkpoints.items():bind(bindings,ap.parent/name,d)
    verify(bindings)
    archive=rp.parent/'completion_archive.json';inner=rp.parent/'completion_closure_plan.json'
    spec=dict(output=str(archive),completed_status=plan['completed_status']+'_archive',
        evidence=dict(producer=dict(path=str(rp),expected=dict(status=producer['status'],**summary)),
                      reader=dict(path=str(ap),expected=dict(status=reader['status'],**summary))),
        links=[dict(**{'from':'reader'},field='producer_receipt_sha256',to=str(rp))],
        launches=plan['launches'],pins=bindings,summary=summary,scope=plan['scope'])
    with inner.open('x') as f:f.write(json.dumps(spec,indent=2)+'\n')
    subprocess.run([sys.executable,'scripts/record_completed_process_handoffs_v2.py','--plan',str(inner)],check=True)
    proof=json.loads(archive.read_text());assert len(proof['services'])==2
    result=dict(status=plan['completed_status'],**summary,full_hash_archive=str(archive),
        full_hash_archive_sha256=sha(archive),bound_source_hashes=len(proof['source_hashes']),exact_process_journals_checked=2,
        producer_receipt=str(rp),producer_receipt_sha256=sha(rp),independent_readback=str(ap),independent_readback_sha256=sha(ap),
        completion_plan_sha256=sha(a.plan),scientific_eligibility=False,scope=plan['scope'])
    with Path(plan['output']).open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['scope']}),flush=True)


if __name__=='__main__':main()
