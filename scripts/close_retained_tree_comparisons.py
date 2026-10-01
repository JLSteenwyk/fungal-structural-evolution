#!/usr/bin/env python3
"""Close complete retained baseline comparisons after full readback/two journals."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
from retained_tree_comparison_sources import load,SUMMARY_FIELDS
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);args=parser.parse_args()
    plan=json.loads(args.plan.read_text());config=json.loads(Path(plan['source_plan']).read_text());groups,universe,bindings=load(config,plan['source_plan'])
    bind(bindings,args.plan)
    for path,digest in plan['pins'].items():bind(bindings,path,digest)
    root=Path(config['output']);rp=root/'receipt.json';ap=Path(plan['readback']);r,a=[json.loads(p.read_text()) for p in [rp,ap]]
    assert r['status']=='complete_full_retained_baseline_tree_comparisons_pending_readback' and a['status']=='passed_full_retained_baseline_tree_comparison_independent_readback'
    assert r['plan_sha256']==a['plan_sha256']==sha(plan['source_plan']) and a['producer_receipt_sha256']==sha(rp)
    summary={k:r[k] for k in SUMMARY_FIELDS};assert all(a[k]==v for k,v in summary.items()) and r['policy_summaries']==a['policy_summaries']
    assert summary['policies']==4 and summary['tree_views']==32 and summary['comparison_rows']==64
    for p in [rp,ap]:bind(bindings,p)
    for record in [r,a]:
        assert record['scientific_eligibility'] is False
        for path,digest in record['source_hashes'].items():bind(bindings,path,digest)
    for name,digest in r['artifacts'].items():bind(bindings,root/name,digest)
    verify(bindings);archive=root/'completion_archive.json';inner=root/'completion_closure_plan.json'
    closure=dict(output=str(archive),completed_status='complete_verified_retained_tree_comparison_archive',
        evidence=dict(producer=dict(path=str(rp),expected=dict(status=r['status'],**summary)),reader=dict(path=str(ap),expected=dict(status=a['status'],**summary))),
        links=[dict(**{'from':'reader'},field='producer_receipt_sha256',to=str(rp))],launches=plan['launches'],pins=bindings,summary=summary,scope=plan['scope'])
    with inner.open('x') as handle:handle.write(json.dumps(closure,indent=2)+'\n')
    subprocess.run([sys.executable,'scripts/record_completed_order_marker_handoffs.py','--plan',str(inner)],check=True)
    closed=json.loads(archive.read_text());assert len(closed['services'])==2
    result=dict(status='complete_verified_full_retained_baseline_tree_comparisons',**summary,policy_summaries=r['policy_summaries'],
        full_hash_archive=str(archive),full_hash_archive_sha256=sha(archive),bound_source_hashes=len(closed['source_hashes']),exact_process_journals_checked=2,
        producer_receipt=str(rp),producer_receipt_sha256=sha(rp),independent_readback=str(ap),independent_readback_sha256=sha(ap),
        completion_plan_sha256=sha(args.plan),scientific_eligibility=False,scope=plan['scope'])
    with Path(plan['output']).open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
