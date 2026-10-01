#!/usr/bin/env python3
"""Close complete matched-taxon baseline projections with two original journals."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
from retained_taxon_projection_sources import load
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);args=parser.parse_args()
    plan=json.loads(args.plan.read_text());config=json.loads(Path(plan['source_plan']).read_text())
    sources,universe,policies,bindings=load(config,plan['source_plan']);bind(bindings,args.plan)
    for path,digest in plan['pins'].items():bind(bindings,path,digest)
    root=Path(config['output']);rp=root/'receipt.json';ap=Path(plan['readback']);r,a=[json.loads(p.read_text()) for p in [rp,ap]]
    assert r['status']=='complete_full_retained_taxon_baseline_projections_pending_independent_readback'
    assert a['status']=='passed_full_retained_taxon_baseline_projection_independent_readback'
    assert r['plan_sha256']==a['plan_sha256']==sha(plan['source_plan']) and a['producer_receipt_sha256']==sha(rp)
    fields=['original_baseline_runs','original_raw_bootstrap_trees','policy_baseline_combinations',
            'projected_bootstrap_tree_states','tree_views','edge_rows','bootstrap_split_rows',
            'role_boundary_rows','boundary_views_with_role_split']
    summary={k:r[k] for k in fields};assert all(a[k]==v for k,v in summary.items()) and r['summaries']==a['summaries']
    assert summary['original_baseline_runs']==4 and summary['original_raw_bootstrap_trees']==4000 and summary['policy_baseline_combinations']==16
    assert summary['projected_bootstrap_tree_states']==16000 and summary['tree_views']==summary['role_boundary_rows']==32 and summary['edge_rows']==33088
    for p in [rp,ap]:bind(bindings,p)
    for record in [r,a]:
        assert record['scientific_eligibility'] is False
        for path,digest in record['source_hashes'].items():bind(bindings,path,digest)
    for name,digest in r['artifacts'].items():bind(bindings,root/name,digest)
    verify(bindings)
    archive=root/'completion_archive.json';inner=root/'completion_closure_plan.json'
    closure=dict(output=str(archive),completed_status='complete_verified_retained_taxon_projection_archive',
        evidence=dict(producer=dict(path=str(rp),expected=dict(status=r['status'],**summary)),
                      reader=dict(path=str(ap),expected=dict(status=a['status'],**summary))),
        links=[dict(**{'from':'reader'},field='producer_receipt_sha256',to=str(rp))],
        launches=plan['launches'],pins=bindings,summary=summary,scope=plan['scope'])
    with inner.open('x') as handle:handle.write(json.dumps(closure,indent=2)+'\n')
    subprocess.run([sys.executable,'scripts/record_completed_order_marker_handoffs.py','--plan',str(inner)],check=True)
    closed=json.loads(archive.read_text());assert len(closed['services'])==2
    result=dict(status='complete_verified_full_retained_taxon_baseline_projections',**summary,full_hash_archive=str(archive),
        full_hash_archive_sha256=sha(archive),bound_source_hashes=len(closed['source_hashes']),exact_process_journals_checked=2,
        producer_receipt=str(rp),producer_receipt_sha256=sha(rp),independent_readback=str(ap),
        independent_readback_sha256=sha(ap),completion_plan_sha256=sha(args.plan),scientific_eligibility=False,scope=plan['scope'])
    with Path(plan['output']).open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
