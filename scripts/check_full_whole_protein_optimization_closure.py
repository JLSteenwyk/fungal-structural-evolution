#!/usr/bin/env python3
"""Contracts for complete original audit coverage and truthful full-follow-up provenance."""
import argparse
import copy
import json
from pathlib import Path
import tempfile
from close_full_whole_protein_optimization import original_proofs,followup_proofs
from whole_protein_flag_followup import full_scope,PASS,FLAG,ERROR
from run_whole_protein_ml import digest
from screen_duplication_alignment_reuse import sha


def put(path,value):
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value)+'\n')


def rejected(function):
    try:function()
    except (AssertionError,KeyError,ValueError):return
    raise AssertionError('False closure evidence accepted')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    with tempfile.TemporaryDirectory() as temp:
        root=Path(temp);inputs=root/'inputs';audit=root/'audit';follow=root/'follow';checked=root/'checked'
        recipes={};production={};audited={};trees={str(i) for i in range(5)}
        for i in range(3):
            spec=dict(columns=['y','intercept'],records=8,fixture=i);identifier=digest(spec);path=inputs/(identifier+'.npz');put(path,spec)
            recipes[identifier]=dict(path=path.name,sha256=sha(path),recipe=dict(specification=spec));results=[]
            for tree in sorted(trees):
                status=FLAG if i==0 and tree in ['0','1','2'] else ERROR if i==1 and tree=='0' else PASS
                f=root/'fits'/(identifier+'-'+tree+'.json');put(f,dict(status=status));key=(identifier,tree)
                production[key]=dict(fit_input_id=identifier,tree=tree,path=str(f),sha256=sha(f),status=status)
                audited[key]=dict(fit_input_id=identifier,tree=tree,status=status,numerical_fit_verified=status!=ERROR,candidates_checked=0 if status==ERROR else 22);results.append(audited[key])
            identity=dict(fit_input_id=identifier,input_sha256=sha(path),fit_hashes={production[identifier,t]['path']:production[identifier,t]['sha256'] for t in trees},fit_plan_sha256='fit',audit_plan_sha256='audit')
            put(audit/'inputs'/identifier[:2]/(identifier+'.json'),dict(identity=identity,results=results,results_sha256=digest(results)))
        production,flags,errors=full_scope(production.values(),audited.values(),trees,set(recipes));assert len(production)==15 and len(flags)==3 and len(errors)==1
        bindings={};assert original_proofs(production,audited,recipes,inputs,audit,'fit','audit',bindings)==(3,15)
        identifier=next(iter(recipes));proofpath=audit/'inputs'/identifier[:2]/(identifier+'.json');baseline=json.loads(proofpath.read_text())
        def mutate_original(mutator):
            bad=copy.deepcopy(baseline);mutator(bad);bad['results_sha256']=digest(bad['results']);put(proofpath,bad)
            rejected(lambda:original_proofs(production,audited,recipes,inputs,audit,'fit','audit',{}));put(proofpath,baseline)
        mutate_original(lambda r:r['results'].pop())
        mutate_original(lambda r:r['results'].append(copy.deepcopy(r['results'][0])))
        mutate_original(lambda r:r['identity'].update(input_sha256='changed'))
        mutate_original(lambda r:r['results'][0].update(numerical_fit_verified=False))
        mutate_original(lambda r:r['identity'].update(fit_hashes={}))
        scope=[]
        for r in production.values():
            status='optimization_flag_followup_required' if r['status']==FLAG else 'original_fit_error_retained_requires_review' if r['status']==ERROR else 'original_numerical_pass_retained'
            scope.append(dict(**r,followup_scope_status=status,scientific_eligibility=False))
        follow.mkdir();checked.mkdir();(follow/'scope_dispositions.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in scope))
        statuses=['numerical_followup_passed_pending_independent_readback','numerical_followup_requires_review','numerical_followup_error_requires_review'];manifest=[];readbacks=[]
        for r,status in zip(flags,statuses):
            identity=dict(fit_input_id=r['fit_input_id'],tree=r['tree'],original_fit=r['path'],original_fit_sha256=r['sha256'],input_sha256=recipes[r['fit_input_id']]['sha256'],plan_sha256='follow')
            path=follow/(r['tree']+'.json');result=dict(status=status);put(path,dict(identity=identity,original_status=FLAG,scientific_eligibility=False,result=result,result_sha256=digest(result)))
            manifest.append(dict(fit_input_id=r['fit_input_id'],tree=r['tree'],path=path.name,sha256=sha(path),status=status))
            readbacks.append(dict(fit_input_id=r['fit_input_id'],tree=r['tree'],source_sha256=sha(path),status=status,numerical_result_verified=status!=statuses[2]))
        def save(rows,path):path.write_text(''.join(json.dumps(r)+'\n' for r in rows))
        save(manifest,follow/'followup_manifest.jsonl');save(readbacks,checked/'readback_manifest.jsonl')
        check=lambda:followup_proofs(production,flags,recipes,follow,checked,'follow',{})
        assert check()==(3,{s:1 for s in statuses})
        for change in [lambda r:r.pop(),lambda r:r.append(copy.deepcopy(r[0])),lambda r:r[0].update(source_sha256='changed'),lambda r:r[2].update(numerical_result_verified=True)]:
            bad=copy.deepcopy(readbacks);change(bad);save(bad,checked/'readback_manifest.jsonl');rejected(check);save(readbacks,checked/'readback_manifest.jsonl')
        bad=copy.deepcopy(scope);bad.pop();save(bad,follow/'scope_dispositions.jsonl');rejected(check);save(scope,follow/'scope_dispositions.jsonl')
        path=follow/manifest[0]['path'];baseline=json.loads(path.read_text());bad=copy.deepcopy(baseline);bad['identity']['original_fit_sha256']='changed';put(path,bad)
        changed=copy.deepcopy(manifest);changed[0]['sha256']=sha(path);save(changed,follow/'followup_manifest.jsonl')
        changedproof=copy.deepcopy(readbacks);changedproof[0]['source_sha256']=sha(path);save(changedproof,checked/'readback_manifest.jsonl');rejected(check)
    result=dict(status='passed_full_whole_protein_optimization_closure_contracts',original_inputs=3,original_fit_audit_states=15,flagged_cases=3,original_fit_errors_retained=1,
        passed_review_and_error_followups_preserved=True,rejected_false_evidence_cases=11,source_hashes={str(Path(__file__).relative_to(Path.cwd())):sha(__file__),'scripts/close_full_whole_protein_optimization.py':sha('scripts/close_full_whole_protein_optimization.py')},scientific_eligibility=False,scope='Software provenance/census contracts; production full375350-grid and original journals are not stubbed or claimed complete.')
    args.output.parent.mkdir(exist_ok=True,parents=True)
    with args.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
