#!/usr/bin/env python3
"""Close all 30 native coalescent quartet/numeric audits with actual journals."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    args=parser.parse_args();plan=json.loads(args.plan.read_text())
    source_plan=json.loads(Path(plan['source_plan']).read_text())
    assert source_plan['mode']=='full_batch'
    native_plan=json.loads(Path(source_plan['native_plan']).read_text())
    root=Path(source_plan['output'])
    np=Path(native_plan['output'])/'receipt.json';ap=root/'receipt.json'
    native,reader=[json.loads(p.read_text()) for p in [np,ap]]
    assert native['status']=='complete_full_native_coalescent_species_sensitivities_pending_independent_scoring'
    assert reader['status']=='passed_full_30_coalescent_global_local_quartet_numeric_readback'
    assert reader['plan_sha256']==sha(plan['source_plan'])
    assert reader['native_plan_sha256']==native['plan_sha256']==sha(source_plan['native_plan'])
    assert reader['global_numerator_independently_verified'] is True
    assert reader['native_search_optimality_verified'] is False
    assert reader['scientific_eligibility'] is native['scientific_eligibility'] is False
    expected=dict(cases=30,markers_per_case=125,internal_branches=15510,branch_gene_states=1938750,
                  global_numerator_independently_verified=True,native_search_optimality_verified=False)
    assert all(reader[k]==v for k,v in expected.items())
    assert native['cases']==len(native['runs'])==30
    assert [r['case']['case'] for r in native['runs']]==[r['case'] for r in reader['summaries']]
    bindings=dict(plan['pins']);bind(bindings,args.plan);bind(bindings,plan['source_plan'])
    for p,r in [(np,native),(ap,reader)]:
        bind(bindings,p)
        for path,digest in r['source_hashes'].items():bind(bindings,path,digest)
    assert reader['source_hashes'][str(np)]==sha(np)
    for row in reader['summaries']:
        assert row['global_numerator_independently_verified'] is True
        assert row['native_search_optimality_verified'] is False
        bind(bindings,row['local_readback_receipt'],row['local_readback_receipt_sha256'])
        bind(bindings,row['global_readback'],row['global_readback_sha256'])
    verify(bindings)
    archive=root/'completion_archive.json';inner=root/'completion_closure_plan.json'
    closure=dict(output=str(archive),completed_status='complete_verified_full_native_coalescent_quartet_numeric_archive',
        evidence=dict(native=dict(path=str(np),expected=dict(status=native['status'],cases=30,markers_per_case=125)),
                      reader=dict(path=str(ap),expected=dict(status=reader['status'],**expected))),
        links=[dict(**{'from':'reader'},field='native_plan_sha256',to=source_plan['native_plan'])],
        launches=plan['launches'],pins=bindings,summary=expected,scope=plan['scope'])
    with inner.open('x') as f:f.write(json.dumps(closure,indent=2)+'\n')
    subprocess.run([sys.executable,'scripts/record_completed_process_handoffs_v2.py','--plan',str(inner)],check=True)
    proof=json.loads(archive.read_text());assert len(proof['services'])==2
    result=dict(status='complete_verified_full_native_coalescent_quartet_numerics',**expected,
        full_hash_archive=str(archive),full_hash_archive_sha256=sha(archive),
        bound_source_hashes=len(proof['source_hashes']),exact_process_journals_checked=2,
        native_receipt=str(np),native_receipt_sha256=sha(np),
        independent_readback=str(ap),independent_readback_sha256=sha(ap),
        completion_plan_sha256=sha(args.plan),scientific_eligibility=False,scope=plan['scope'])
    with Path(plan['output']).open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if isinstance(v,(str,int,bool))}),flush=True)


if __name__=='__main__':main()
