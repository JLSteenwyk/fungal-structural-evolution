#!/usr/bin/env python3
"""Close complete coalescent/reference comparisons and both original journals."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

from coalescent_tree_comparison_sources import load,SUMMARY_FIELDS
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);args=parser.parse_args()
    plan=json.loads(args.plan.read_text());source=json.loads(Path(plan['source_plan']).read_text())
    groups,universe,bindings=load(source,plan['source_plan']);bind(bindings,args.plan)
    for p,d in plan['pins'].items():bind(bindings,p,d)
    root=Path(source['output']);rp=root/'receipt.json';ap=root/'readback.json'
    producer,reader=[json.loads(p.read_text()) for p in [rp,ap]]
    assert producer['status']=='complete_coalescent_reference_tree_comparisons_pending_independent_readback'
    assert reader['status']=='passed_coalescent_reference_full_raw_tree_comparison_readback'
    assert producer['plan_sha256']==reader['plan_sha256']==sha(plan['source_plan'])
    assert reader['producer_receipt_sha256']==sha(rp)
    summary={k:reader[k] for k in SUMMARY_FIELDS};assert all(producer[k]==v for k,v in summary.items())
    if source['mode']=='full_batch':
        assert summary['cohorts']==5 and summary['tree_views']==summary['role_boundary_rows']==70 and summary['comparison_rows']==315
    else:
        assert source['mode']=='complete_named_case_preflight'
        assert summary['cohorts']==1 and summary['tree_views']==summary['role_boundary_rows']==9 and summary['comparison_rows']==8
    for p,r in [(rp,producer),(ap,reader)]:
        bind(bindings,p);assert r['scientific_eligibility'] is False
        for path,digest in r['source_hashes'].items():bind(bindings,path,digest)
    for name,digest in producer['artifacts'].items():bind(bindings,root/name,digest)
    verify(bindings)
    archive=root/'completion_archive.json';inner=root/'completion_closure_plan.json'
    status=('complete_verified_full_coalescent_reference_comparison_archive' if source['mode']=='full_batch' else
            'complete_verified_named_coalescent_reference_comparison_archive')
    closure=dict(output=str(archive),completed_status=status,
        evidence=dict(producer=dict(path=str(rp),expected=dict(status=producer['status'],**summary)),
                      reader=dict(path=str(ap),expected=dict(status=reader['status'],**summary))),
        links=[dict(**{'from':'reader'},field='producer_receipt_sha256',to=str(rp))],
        launches=plan['launches'],pins=bindings,summary=summary,scope=plan['scope'])
    with inner.open('x') as handle:handle.write(json.dumps(closure,indent=2)+'\n')
    subprocess.run([sys.executable,'scripts/record_completed_process_handoffs_v2.py','--plan',str(inner)],check=True)
    proof=json.loads(archive.read_text());assert len(proof['services'])==2
    result=dict(status=('complete_verified_full_coalescent_reference_tree_comparisons' if source['mode']=='full_batch' else
                        'complete_verified_named_coalescent_reference_tree_comparisons'),**summary,
        full_hash_archive=str(archive),full_hash_archive_sha256=sha(archive),bound_source_hashes=len(proof['source_hashes']),
        exact_process_journals_checked=2,producer_receipt=str(rp),producer_receipt_sha256=sha(rp),
        independent_readback=str(ap),independent_readback_sha256=sha(ap),completion_plan_sha256=sha(args.plan),
        scientific_eligibility=False,scope=plan['scope'])
    with Path(plan['output']).open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if isinstance(v,(str,int,bool))}),flush=True)


if __name__=='__main__':main()
