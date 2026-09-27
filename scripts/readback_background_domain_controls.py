#!/usr/bin/env python3
"""Independently read back every background domain annotation and pair control."""
import argparse,csv,hashlib,json,sqlite3,time
from collections import Counter,defaultdict
from pathlib import Path
import psutil
from run_ortholog_pair_guide_comparison import sha
from readback_duplication_domain_controls import POLICIES,pair_summary

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        if sha(a.plan)!=ph:raise ValueError('Changed readback plan')
        for path,h in plan['pins'].items():
            if sha(path)!=h:raise ValueError('Changed source: '+path)
    verify();dep=plan.get('producer')
    while dep is not None:
        try:
            proc=psutil.Process(dep['pid'])
            if proc.create_time()!=dep['created'] or proc.status()==psutil.STATUS_ZOMBIE:break
            if proc.cmdline()!=dep['cmdline']:raise ValueError('Producer identity changed')
        except psutil.NoSuchProcess:break
        time.sleep(30)
    verify();sourceplan=json.loads(Path(plan['source_plan']).read_text());folder=Path(sourceplan['output']);receipt=json.loads((folder/'receipt.json').read_text())
    if receipt['status']!='complete_background_domain_control_inventory_pending_readback' or receipt['plan_sha256']!=sha(plan['source_plan']):raise ValueError('Unbound producer')
    for name,h in receipt['artifacts'].items():
        if sha(folder/name)!=h:raise ValueError('Changed output artifact')
    models={};pairs={};scopes=defaultdict(set)
    for label,root in [('background',Path(sourceplan['inventory']))]:
        r=json.loads((root/'receipt.json').read_text())
        for name in ['models.jsonl','model_pairs.tsv']:
            if sha(root/name)!=r['artifacts'][name]:raise ValueError('Changed queue artifact')
        with (root/'models.jsonl').open() as f:
            for line in f:
                m=json.loads(line);key=m['model_id'],m['version']
                if key in models and models[key]!=m:raise ValueError('Inconsistent shared model')
                models[key]=m
        with (root/'model_pairs.tsv').open() as f:
            for row in csv.DictReader(f,delimiter='\t'):
                endpoints=tuple(sorted([(row['model_a'],int(row['version_a'])),(row['model_b'],int(row['version_b']))]))
                key=hashlib.sha256(json.dumps(endpoints,separators=(',',':')).encode()).hexdigest()
                if key!=row['pair_key'] or label in scopes[key] or endpoints[0]==endpoints[1]:raise ValueError('Invalid pair membership')
                pairs[key]=endpoints;scopes[key].add(label)
    db=sqlite3.connect('file:'+str((Path(sourceplan['registry'])/'structure_domains.sqlite').resolve())+'?mode=ro',uri=True);db.row_factory=sqlite3.Row
    seen=set();annotations={};counts=Counter()
    with (folder/'model_annotations.jsonl').open() as f:
        for line in f:
            row=json.loads(line);key=row['model_id'],row['version']
            if key not in models or key in seen:raise ValueError('Unexpected/repeated model')
            m=models[key];native=db.execute('SELECT * FROM models WHERE sequence_sha256=?',(m['sequence_sha256'],)).fetchone()
            if native is None or (native['model_id'],native['version'],native['length'],native['path'])!=(key[0],key[1],m['length'],m['path']):raise ValueError('Registry model provenance differs')
            segments={h['hit_id']:dict(h) for h in db.execute('SELECT * FROM segments WHERE model_key=?',(native['model_key'],))}
            expected={policy:[] for policy in POLICIES}
            for hit in db.execute('SELECT * FROM policy_hits WHERE model_key=?',(native['model_key'],)):
                h=dict(segments[hit['hit_id']]);h.pop('model_key');h.update(conservative_architecture=hit['conservative_architecture'],domain_interval_candidate=hit['domain_interval_candidate'])
                expected[hit['policy']].append(h)
            for hits in expected.values():hits.sort(key=lambda h:(h['alignment_start'],h['alignment_end'],h['hit_id']))
            entire=dict(model_id=key[0],version=key[1],sequence_sha256=m['sequence_sha256'],raw_hits=native['raw_hits'],policies_agree=native['policies_agree'],policies=expected)
            if row!=entire:raise ValueError('Exported model annotation differs')
            annotations[key]=expected;seen.add(key);counts['models']+=1;counts['models_without_hits']+=not native['raw_hits'];counts['policy_disagreement_models']+=not native['policies_agree']
            if len(seen)%10000==0:print('Checked models',len(seen),flush=True)
    db.close()
    if seen!=set(models) or dict(counts)!=receipt['models']:raise ValueError('Incomplete model universe')
    architecture_counts=Counter();domain_counts=Counter();total=0
    with (folder/'pair_architecture_controls.tsv').open() as pf,(folder/'single_copy_domain_pairs.tsv').open() as df:
        pairrows=iter(csv.DictReader(pf,delimiter='\t'));domainrows=iter(csv.DictReader(df,delimiter='\t'))
        for key,(left,right) in sorted(pairs.items()):
            scope='+'.join(sorted(scopes[key]))
            for policy in POLICIES:
                values,matches=pair_summary(annotations[left][policy],annotations[right][policy]);expected=dict(pair_key=key,scope=scope,policy=policy,**values)
                actual=next(pairrows,None)
                if actual!={k:str(v) for k,v in expected.items()}:raise ValueError('Pair control record differs')
                total+=1;architecture_counts[policy+':'+values['annotation_class']]+=1;domain_counts[policy]+=len(matches)
                for accession,x,y in matches:
                    expected=dict(pair_key=key,scope=scope,policy=policy,pfam_accession=accession,model_left=left[0],version_left=left[1],hit_left=x['hit_id'],alignment_start_left=x['alignment_start'],alignment_end_left=x['alignment_end'],model_right=right[0],version_right=right[1],hit_right=y['hit_id'],alignment_start_right=y['alignment_start'],alignment_end_right=y['alignment_end'])
                    if next(domainrows,None)!={k:str(v) for k,v in expected.items()}:raise ValueError('Single-copy domain match differs')
        if next(pairrows,None) is not None or next(domainrows,None) is not None:raise ValueError('Extra exported record')
    if total!=receipt['pair_policy_rows'] or len(pairs)!=receipt['unique_model_pairs'] or dict(architecture_counts)!=receipt['architecture_counts'] or dict(domain_counts)!=receipt['single_copy_domain_matches_by_policy']:raise ValueError('Summary counts differ')
    if dict(Counter('+'.join(sorted(s)) for s in scopes.values()))!=receipt['pair_scope_counts']:raise ValueError('Scope counts differ')
    verify()
    for name,h in receipt['artifacts'].items():
        if sha(folder/name)!=h:raise ValueError('Output changed during readback')
    result=dict(status='passed_full_background_domain_control_readback',plan_sha256=ph,producer_receipt_sha256=sha(folder/'receipt.json'),models=len(models),unique_model_pairs=len(pairs),pair_policy_rows=total,single_copy_domain_matches_by_policy=dict(domain_counts),scope='Every model export independently rejoined from registry segments and policy membership; full verified background pair universe, all architecture categories/counts and every single-copy domain match independently reconstructed. Same frozen Pfam source, not independent annotation searches or structural domain validation. No biological gain/loss/rearrangement or duplication-effect inference.')
    with Path(plan['output']).open('x') as f:json.dump(result,f,indent=2);f.write('\n')


if __name__=='__main__':main()
