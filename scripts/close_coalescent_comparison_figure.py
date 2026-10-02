#!/usr/bin/env python3
"""Close the full 315-cell figure readback and original two process journals."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);args=parser.parse_args()
    plan=json.loads(args.plan.read_text());source=json.loads(Path(plan['source_plan']).read_text())
    root=Path(source['output']);bindings=dict(plan['pins']);bind(bindings,args.plan);bind(bindings,plan['source_plan'])
    rp,ap=root/'receipt.json',root/'readback.json'
    producer,reader=[json.loads(p.read_text()) for p in [rp,ap]]
    assert producer['status']=='complete_full_coalescent_reference_comparison_figure_pending_cell_readback'
    assert reader['status']=='passed_full_315_coalescent_comparison_figure_cell_and_pdf_readback'
    assert producer['plan_sha256']==reader['plan_sha256']==sha(plan['source_plan'])
    assert reader['producer_receipt_sha256']==sha(rp)
    summary=dict(comparisons=315,cross_reference_cells=240,internal_candidate_cells=75,
                 panels=10,pdf_pages=2,tree_views=70,cohorts=5)
    assert all(producer[k]==reader[k]==value for k,value in summary.items())
    for p,r in [(rp,producer),(ap,reader)]:
        bind(bindings,p);assert r['scientific_eligibility'] is False
        for path,digest in r['source_hashes'].items():bind(bindings,path,digest)
    for name,digest in producer['artifacts'].items():bind(bindings,root/name,digest)
    verify(bindings)
    archive=root/'completion_archive.json';inner=root/'completion_closure_plan.json'
    closure=dict(output=str(archive),completed_status='complete_verified_full_coalescent_comparison_figure_archive',
        evidence=dict(producer=dict(path=str(rp),expected=dict(status=producer['status'],**summary)),
                      reader=dict(path=str(ap),expected=dict(status=reader['status'],**summary))),
        links=[dict(**{'from':'reader'},field='producer_receipt_sha256',to=str(rp))],
        launches=plan['launches'],pins=bindings,summary=summary,scope=plan['scope'])
    with inner.open('x') as f:f.write(json.dumps(closure,indent=2)+'\n')
    subprocess.run([sys.executable,'scripts/record_completed_process_handoffs_v2.py','--plan',str(inner)],check=True)
    proof=json.loads(archive.read_text());assert len(proof['services'])==2
    result=dict(status='complete_verified_full_315_coalescent_comparison_figure_pending_visual_inspection',**summary,
        full_hash_archive=str(archive),full_hash_archive_sha256=sha(archive),bound_source_hashes=len(proof['source_hashes']),
        exact_process_journals_checked=2,producer_receipt=str(rp),producer_receipt_sha256=sha(rp),
        independent_readback=str(ap),independent_readback_sha256=sha(ap),completion_plan_sha256=sha(args.plan),
        comparison_completion=source['completion'],comparison_completion_sha256=sha(source['completion']),
        visual_inspection_required=True,scientific_eligibility=False,scope=plan['scope'])
    with Path(plan['output']).open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if isinstance(v,(str,int,bool))}),flush=True)


if __name__=='__main__':main()
