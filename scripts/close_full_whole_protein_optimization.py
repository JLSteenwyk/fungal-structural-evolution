#!/usr/bin/env python3
"""Bind every original fit/audit and all follow-ups before four-original-journal closure."""
import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess
import sys
import shutil
from run_whole_protein_ml import digest
from whole_protein_flag_followup import full_scope,FLAG,ERROR
from reference_measurement_union_sources import bind,verify
from screen_duplication_alignment_reuse import sha


def index(path):
    records={}
    with Path(path).open() as f:
        for line in f:
            r=json.loads(line);key=r['fit_input_id'],r['tree'];assert key not in records;records[key]=r
    return records


def original_proofs(production,audited,recipes,inputs,audit_root,fit_hash,audit_hash,bindings):
    """Every original five-tree audit shard must match both complete manifests."""
    trees={k[1] for k in production};seen=set()
    for identifier,recipe in recipes.items():
        assert digest(recipe['recipe']['specification'])==identifier
        bind(bindings,Path(inputs)/recipe['path'],recipe['sha256'])
        path=Path(audit_root)/'inputs'/identifier[:2]/(identifier+'.json');bind(bindings,path)
        proof=json.loads(path.read_text());keys={(identifier,t) for t in trees}
        assert keys<=production.keys() and keys<=audited.keys()
        expected=dict(fit_input_id=identifier,input_sha256=recipe['sha256'],
            fit_hashes={production[k]['path']:production[k]['sha256'] for k in keys},
            fit_plan_sha256=fit_hash,audit_plan_sha256=audit_hash)
        assert proof['identity']==expected and proof['results_sha256']==digest(proof['results'])
        results={(r['fit_input_id'],r['tree']):r for r in proof['results']}
        assert len(results)==len(proof['results'])==len(trees) and set(results)==keys
        for key in keys:
            assert results[key]==audited[key]
            assert results[key]['numerical_fit_verified'] is (production[key]['status']!=ERROR)
            bind(bindings,production[key]['path'],production[key]['sha256'])
        seen.update(keys)
    assert seen==set(production)==set(audited)
    return len(recipes),len(seen)


def followup_proofs(production,flags,recipes,root,checked_root,plan_hash,bindings):
    scope=index(Path(root)/'scope_dispositions.jsonl');assert set(scope)==set(production)
    for key,row in production.items():
        expected='optimization_flag_followup_required' if row['status']==FLAG else 'original_fit_error_retained_requires_review' if row['status']==ERROR else 'original_numerical_pass_retained'
        assert scope[key]==dict(**row,followup_scope_status=expected,scientific_eligibility=False)
    manifest=index(Path(root)/'followup_manifest.jsonl');checked=index(Path(checked_root)/'readback_manifest.jsonl')
    expected={(r['fit_input_id'],r['tree']) for r in flags};assert set(manifest)==set(checked)==expected
    counts=Counter()
    for key,row in manifest.items():
        path=Path(root)/row['path'];bind(bindings,path,row['sha256']);saved=json.loads(path.read_text());original=production[key]
        assert saved['identity']==dict(fit_input_id=key[0],tree=key[1],original_fit=original['path'],original_fit_sha256=original['sha256'],input_sha256=recipes[key[0]]['sha256'],plan_sha256=plan_hash)
        assert saved['original_status']==FLAG and saved['scientific_eligibility'] is False
        assert saved['result_sha256']==digest(saved['result'])
        status=saved['result']['status'];assert status==row['status']==checked[key]['status']
        assert status in ['numerical_followup_passed_pending_independent_readback','numerical_followup_requires_review','numerical_followup_error_requires_review']
        assert checked[key]['source_sha256']==row['sha256']
        assert checked[key]['numerical_result_verified'] is (status!='numerical_followup_error_requires_review')
        counts[status]+=1
    return len(manifest),dict(counts)


def run(path):
    plan=json.loads(Path(path).read_text());bindings=dict(plan['pins']);bind(bindings,path)
    assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
    fp,ap,pp,vp=[json.loads(Path(plan[k]).read_text()) for k in ['fit_plan','audit_plan','followup_plan','readback_plan']]
    assert pp['fit_plan']==plan['fit_plan'] and pp['audit_plan']==plan['audit_plan'] and vp['source_plan']==plan['followup_plan']
    roots=[Path(r['output']) for r in [fp,ap,pp,vp]];receipts=[]
    for root in roots:
        p=root/'receipt.json';bind(bindings,p);r=json.loads(p.read_text());receipts.append(r)
        for key in ['source_hashes']:
            for p,h in r.get(key,{}).items():bind(bindings,p,h)
        for n,h in r['artifacts'].items():bind(bindings,root/n,h)
    fit,audit,producer,reader=receipts
    assert fit['status']=='complete_whole_protein_ml_dispositions_pending_full_audit' and fit['plan_sha256']==sha(plan['fit_plan'])
    assert audit['status']=='complete_full_whole_protein_ml_output_audit' and audit['audit_plan_sha256']==sha(plan['audit_plan'])
    assert audit['source_receipt_sha256']==sha(roots[0]/'receipt.json')
    assert fit['unique_inputs']==75070 and fit['tree_fit_dispositions']==audit['tree_fit_dispositions']==375350
    assert fit['status_counts']==audit['status_counts']
    assert producer['status']=='complete_full_grid_optimization_followup_pending_independent_readback' and producer['plan_sha256']==sha(plan['followup_plan'])
    assert reader['status']=='passed_full_grid_optimization_followup_readback_with_review_statuses_retained'
    assert reader['source_receipt_sha256']==sha(roots[2]/'receipt.json')
    assert producer['scientific_eligibility'] is reader['scientific_eligibility'] is False
    for root in [roots[0],roots[2]]:
        for p,h in json.loads((root/'source_bindings.json').read_text()).items():bind(bindings,p,h)
    ip=Path(fp['inputs']);ir=json.loads((ip/'receipt.json').read_text());bind(bindings,ip/'receipt.json')
    assert ir['inputs']==75070 and ir['tree_fits']==375350
    for p,h in ir['source_hashes'].items():bind(bindings,p,h)
    for n,h in ir['artifacts'].items():bind(bindings,ip/n,h)
    recipes={}
    for line in (ip/'input_manifest.jsonl').open():
        r=json.loads(line);identifier=r['fit_input_id'];assert identifier not in recipes;recipes[identifier]=r
    assert len(recipes)==75070
    trees={p.stem for p in Path(fp['factors']).glob('*.npz')};assert len(trees)==5
    production=index(roots[0]/'fit_manifest.jsonl');audited=index(roots[1]/'audit_manifest.jsonl')
    production,flags,errors=full_scope(production.values(),audited.values(),trees,set(recipes))
    assert dict(Counter(r['status'] for r in production.values()))==fit['status_counts']
    inputs,total=original_proofs(production,audited,recipes,fp['inputs'],roots[1],sha(plan['fit_plan']),sha(plan['audit_plan']),bindings)
    nf,counts=followup_proofs(production,flags,recipes,roots[2],roots[3],sha(plan['followup_plan']),bindings)
    assert total==producer['full_dispositions']==reader['full_dispositions']==375350
    assert nf==len(flags)==producer['flagged_fits']==reader['flagged_fits']
    assert len(errors)==producer['original_fit_errors_retained']==reader['original_fit_errors_retained']
    assert counts==producer['followup_status_counts']==reader['counts']
    verify(bindings)
    summary=dict(unique_inputs=inputs,full_dispositions=total,original_fit_status_counts=fit['status_counts'],flagged_fits=nf,
        original_fit_errors_retained=len(errors),followup_status_counts=counts,candidate_likelihoods_replayed=reader['candidate_likelihoods_replayed'],
        original_candidate_likelihoods_checked=audit['candidate_likelihoods_checked'],maximum_objective_error=reader['maximum_objective_error'])
    out=Path(plan['archive_directory']);out.mkdir(exist_ok=True,parents=True);inner=out/'completion_plan.json';archive=out/'completion_archive.json'
    labels=['fit','audit','followup','readback'];evidence={label:dict(path=str(root/'receipt.json'),expected=dict(status=r['status'])) for label,root,r in zip(labels,roots,receipts)}
    spec=dict(output=str(archive),completed_status='complete_verified_full_whole_protein_optimization_archive',evidence=evidence,
        links=[dict(**{'from':'audit'},field='source_receipt_sha256',to=str(roots[0]/'receipt.json')),dict(**{'from':'readback'},field='source_receipt_sha256',to=str(roots[2]/'receipt.json'))],
        launches=plan['launches'],pins=bindings,summary=summary,scope=plan['scope'])
    with inner.open('x') as f:f.write(json.dumps(spec,indent=2)+'\n')
    subprocess.run([sys.executable,'scripts/record_completed_process_handoffs_v2.py','--plan',str(inner)],check=True)
    proof=json.loads(archive.read_text());assert len(proof['services'])==4
    assert sum(p.stat().st_size for p in out.iterdir() if p.is_file())<=plan['resources']['maximum_output_gib']*2**30
    result=dict(status='complete_verified_full_whole_protein_optimization',**summary,full_hash_archive=str(archive),full_hash_archive_sha256=sha(archive),
        bound_source_hashes=len(proof['source_hashes']),exact_process_journals_checked=4,producer_receipt=str(roots[2]/'receipt.json'),producer_receipt_sha256=sha(roots[2]/'receipt.json'),
        independent_readback=str(roots[3]/'receipt.json'),independent_readback_sha256=sha(roots[3]/'receipt.json'),completion_plan_sha256=sha(path),scientific_eligibility=False,scope=plan['scope'])
    with Path(plan['output']).open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);run(p.parse_args().plan)
