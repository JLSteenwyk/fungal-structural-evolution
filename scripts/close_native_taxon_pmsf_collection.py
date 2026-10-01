#!/usr/bin/env python3
"""Close all 16 native sensitivities after independent readback and two journals."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);args=parser.parse_args()
    plan=json.loads(args.plan.read_text());bindings=dict(plan['pins']);bind(bindings,args.plan)
    cp=Path(plan['collection_plan']);config=json.loads(cp.read_text());np=Path(config['native_plan']);native=json.loads(np.read_text())
    root=Path(native['output']);rp=root/'receipt.json';ap=Path(config['output'])/'receipt.json'
    r,a=[json.loads(p.read_text()) for p in [rp,ap]]
    assert r['status']=='complete_native_taxon_pmsf_batch_pending_full_independent_collection_readback' and r['plan_sha256']==sha(np)
    assert a['status']=='passed_full_native_taxon_pmsf_independent_collection_pending_journal_closure' and a['plan_sha256']==sha(cp) and a['native_batch_receipt_sha256']==sha(rp)
    assert len(r['runs'])==len(a['runs'])==len(native['jobs'])==a['run_count']==16
    assert a['tree_views']==a['role_boundary_rows']==32 and a['raw_bootstrap_trees']==16000 and a['site_profiles']==902216
    assert r['scientific_eligibility'] is a['scientific_eligibility'] is False
    for job,raw,readback in zip(native['jobs'],r['runs'],a['runs']):
        assert all(raw[k]==readback[k]==v for k,v in job.items())
        folder=root/job['label'];c=json.loads((folder/'config.json').read_text());command=c['command']
        for flag,value in [('-m','LG+C20+F+G4'),('-st','AA'),('-T',str(native['resources']['threads'])),
                           ('--mem',native['resources']['iqtree_mixture_memory_limit']),('--seed',str(native['seed'])),
                           ('--alrt','1000'),('-B','1000')]:
            assert command.count(flag)==1 and command[command.index(flag)+1]==value
        assert command[0]==native['executable'] and '--bnni' in command and '--boot-trees' in command
        assert Path(command[command.index('--prefix')+1]).resolve()==(folder/'pmsf').resolve()
        native_launches=list(folder.glob('pmsf.native_launch*.json'));assert native_launches
        matching=[]
        for p in native_launches:
            captured=json.loads(p.read_text());bind(bindings,p)
            if captured['command']==captured['cmdline']==command and captured['plan_sha256']==sha(np):
                assert isinstance(captured['pid'],int) and captured['created']>0;matching.append(captured)
        assert matching,'No captured native child command for '+job['label']
        checkpoint=root/(job['label']+'.completed.json');assert json.loads(checkpoint.read_text())==raw;bind(bindings,checkpoint)
        assert raw['run_receipt_sha256']==readback['source_receipt_sha256']==sha(folder/'receipt.json')
        assert raw['audit_receipt_sha256']==readback['first_audit_receipt_sha256']==sha(folder/'audit'/'receipt.json')
        assert raw['taxa']==readback['taxa'] and raw['columns']==readback['sites'] and readback['bootstrap_trees']==1000
    for p in [cp,np,rp,ap]:bind(bindings,p)
    for record in [r,a]:
        for path,digest in record['source_hashes'].items():bind(bindings,path,digest)
    for name,digest in a['artifacts'].items():bind(bindings,ap.parent/name,digest)
    verify(bindings)
    summary={k:a[k] for k in ['run_count','tree_views','raw_bootstrap_trees','site_profiles','internal_support_rows',
                            'role_boundary_rows','boundary_views_with_role_split']}
    archive=ap.parent/'completion_archive.json';inner=ap.parent/'completion_closure_plan.json'
    closure=dict(output=str(archive),completed_status='complete_verified_native_taxon_pmsf_collection_archive',
        evidence=dict(native=dict(path=str(rp),expected=dict(status=r['status'],plan_sha256=sha(np))),
                      reader=dict(path=str(ap),expected=dict(status=a['status'],**summary))),
        links=[dict(**{'from':'reader'},field='native_batch_receipt_sha256',to=str(rp))],
        launches=plan['launches'],pins=bindings,summary=summary,scope=plan['scope'])
    with inner.open('x') as handle:handle.write(json.dumps(closure,indent=2)+'\n')
    subprocess.run([sys.executable,'scripts/record_completed_order_marker_handoffs.py','--plan',str(inner)],check=True)
    closed=json.loads(archive.read_text());assert len(closed['services'])==2
    result=dict(status='complete_verified_full_native_taxon_pmsf_collection',**summary,
        full_hash_archive=str(archive),full_hash_archive_sha256=sha(archive),bound_source_hashes=len(closed['source_hashes']),
        exact_process_journals_checked=2,native_batch_receipt=str(rp),native_batch_receipt_sha256=sha(rp),
        independent_readback=str(ap),independent_readback_sha256=sha(ap),completion_plan_sha256=sha(args.plan),
        scientific_eligibility=False,scope=plan['scope'])
    with Path(plan['output']).open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
