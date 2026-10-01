#!/usr/bin/env python3
"""Close full post-screen balance only after complete independent and journal proof."""
import argparse
import json
import subprocess
import sys
from pathlib import Path
from full_screened_balance_sources import load, SUMMARY_FIELDS
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);args=parser.parse_args();plan=json.loads(args.plan.read_text())
    config=json.loads(Path(plan['source_plan']).read_text());source,bindings=load(config,plan['source_plan']);bind(bindings,args.plan)
    for path,digest in plan['pins'].items():bind(bindings,path,digest)
    root=Path(config['output']);rp,ap=root/'receipt.json',Path(plan['readback']);r,a=[json.loads(p.read_text()) for p in [rp,ap]]
    assert r['status']=='complete_full_screened_balance_pending_independent_readback' and a['status']=='passed_full_screened_balance_independent_readback'
    assert r['plan_sha256']==a['plan_sha256']==sha(plan['source_plan']) and a['producer_receipt_sha256']==sha(rp) and r['matching_receipt_sha256']==sha(source['prior_receipt'])
    summary={k:r[k] for k in SUMMARY_FIELDS};assert all(a[k]==v for k,v in summary.items())
    for key in ['target_nodes','background_nodes','selected_records','scenarios']:assert summary[key]==config['expected'][key]
    assert summary['strata']==len(config['guides'])*len(config['policies'])*config['expected']['scenarios']
    assert summary['coverage_rows']==summary['strata']*3*len(config['screens']) and summary['balance_rows']==summary['coverage_rows']*len(config['features'])
    assert summary['retained_selection_screen_cells']==sum(int(row['joint_pass_matched_records']) for row in source['attrition'].values())
    assert 0<=summary['reuse_rows']<=summary['retained_selection_screen_cells']<=summary['selected_records']*3*len(config['screens'])
    for path in [rp,ap]:bind(bindings,path)
    for record in [r,a]:
        assert record['scientific_eligibility'] is False
        for path,digest in record['source_hashes'].items():bind(bindings,path,digest)
    verify(bindings);archive,inner=root/'completion_archive.json',root/'completion_closure_plan.json'
    closure=dict(output=str(archive),completed_status='complete_verified_full_screened_balance_archive',evidence=dict(producer=dict(path=str(rp),expected=dict(status=r['status'],**summary)),reader=dict(path=str(ap),expected=dict(status=a['status'],**summary))),links=[dict(**{'from':'reader'},field='producer_receipt_sha256',to=str(rp))],launches=plan['launches'],pins=bindings,summary=summary,scope=plan['scope'])
    with inner.open('x') as handle:handle.write(json.dumps(closure,indent=2)+'\n')
    subprocess.run([sys.executable,'scripts/record_completed_order_marker_handoffs.py','--plan',str(inner)],check=True);closed=json.loads(archive.read_text());assert len(closed['services'])==2
    result=dict(status='complete_verified_full_fixed_matching_post_screen_balance_and_reuse',**summary,full_hash_archive=str(archive),full_hash_archive_sha256=sha(archive),bound_source_hashes=len(closed['source_hashes']),exact_process_journals_checked=2,producer_receipt=str(rp),producer_receipt_sha256=sha(rp),independent_readback=str(ap),independent_readback_sha256=sha(ap),source_plan=plan['source_plan'],source_plan_sha256=sha(plan['source_plan']),completion_plan_sha256=sha(args.plan),scientific_eligibility=False,scope=plan['scope'])
    with Path(plan['output']).open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
