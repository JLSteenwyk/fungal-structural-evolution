#!/usr/bin/env python3
"""Join verified background models and pairs to the same four domain policies as targets."""
import argparse,csv,hashlib,json,sqlite3
from collections import Counter,defaultdict
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha
from inventory_duplication_domain_controls import POLICIES,architecture_class,single_copy_matches

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        if sha(a.plan)!=ph:raise ValueError('Changed plan')
        for path,h in plan['pins'].items():
            if sha(path)!=h:raise ValueError('Changed source: '+path)
    verify();registry=Path(plan['registry']);rr=json.loads((registry/'receipt.json').read_text());audit=json.loads(Path(plan['registry_readback']).read_text());dbpath=registry/'structure_domains.sqlite'
    if audit['status']!='passed_full_structure_domain_registry_readback' or audit['producer_receipt_sha256']!=sha(registry/'receipt.json') or audit['database_sha256']!=rr['database_sha256'] or sha(dbpath)!=rr['database_sha256']:raise ValueError('Unbound domain registry')
    inventory=Path(plan['inventory']);r=json.loads((inventory/'receipt.json').read_text());proof=json.loads(Path(plan['inventory_readback']).read_text())
    if r['status']!='complete_background_measurement_inventory_pending_readback' or proof['status']!='passed_full_background_measurement_inventory_readback' or proof['producer_receipt_sha256']!=sha(inventory/'receipt.json'):raise ValueError('Unverified background inventory')
    for name,digest in r['artifacts'].items():
        if sha(inventory/name)!=digest:raise ValueError('Changed background artifact')
    models={}
    with (inventory/'models.jsonl').open() as f:
        for line in f:
            m=json.loads(line);key=m['model_id'],m['version']
            if key in models:raise ValueError('Duplicate background model')
            models[key]=m
    if len(models)!=r['all_models'] or len(models)!=proof['all_models']:raise ValueError('Model count differs')
    pairs={};scope=defaultdict(set)
    with (inventory/'model_pairs.tsv').open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            ends=tuple(sorted([(row['model_a'],int(row['version_a'])),(row['model_b'],int(row['version_b']))]))
            key=hashlib.sha256(json.dumps(ends,separators=(',',':')).encode()).hexdigest()
            if key!=row['pair_key'] or key in pairs or ends[0]==ends[1] or any(e not in models for e in ends):raise ValueError('Invalid background pair')
            pairs[key]=ends;scope[key].add('background')
    if len(pairs)!=proof['distinct_eligible_model_pairs']:raise ValueError('Pair scope differs')
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    db=sqlite3.connect('file:'+str(dbpath.resolve())+'?mode=ro',uri=True);db.row_factory=sqlite3.Row
    cache={};model_counts=Counter()
    with (out/'model_annotations.jsonl').open('w') as handle:
        for n,(key,m) in enumerate(sorted(models.items()),1):
            row=db.execute('SELECT * FROM models WHERE sequence_sha256=?',(m['sequence_sha256'],)).fetchone()
            if row is None or (row['model_id'],row['version'],row['length'],row['path'])!=(m['model_id'],m['version'],m['length'],m['path']):raise ValueError('Model-registry provenance differs')
            annotations={policy:[] for policy in POLICIES}
            query='SELECT p.policy,p.conservative_architecture,p.domain_interval_candidate,s.* FROM policy_hits p JOIN segments s USING(model_key,hit_id) WHERE p.model_key=? ORDER BY s.alignment_start,s.alignment_end,s.hit_id'
            for hit in db.execute(query,(row['model_key'],)):
                h=dict(hit);policy=h.pop('policy');h.pop('model_key')
                if policy not in annotations:raise ValueError('Unknown annotation policy')
                annotations[policy].append(h)
            if bool(row['raw_hits'])!=any(annotations.values()):raise ValueError('Missing retained annotations')
            cache[key]=annotations
            model_counts['models']+=1;model_counts['models_without_hits']+=not row['raw_hits'];model_counts['policy_disagreement_models']+=not row['policies_agree']
            handle.write(json.dumps(dict(model_id=key[0],version=key[1],sequence_sha256=m['sequence_sha256'],raw_hits=row['raw_hits'],policies_agree=row['policies_agree'],policies=annotations),separators=(',',':'))+'\n')
            if n%10000==0:print('Annotated models',n,'/',len(models),flush=True)
    counts=Counter();matched=Counter()
    with (out/'pair_architecture_controls.tsv').open('w') as handle,(out/'single_copy_domain_pairs.tsv').open('w') as domainfile:
        fields=['pair_key','scope','policy','annotation_class','left_hits','right_hits','left_conservative','right_conservative','both_conservative','left_domain_hits','right_domain_hits','single_copy_candidate_domain_matches']
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
        dw=csv.writer(domainfile,delimiter='\t',lineterminator='\n');dw.writerow(['pair_key','scope','policy','pfam_accession','model_left','version_left','hit_left','alignment_start_left','alignment_end_left','model_right','version_right','hit_right','alignment_start_right','alignment_end_right'])
        for pair,(left,right) in sorted(pairs.items()):
            label='+'.join(sorted(scope[pair]))
            for policy in POLICIES:
                x,y=cache[left][policy],cache[right][policy]
                signatures=[[(h['pfam_accession'],h['pfam_type']) for h in hits] for hits in [x,y]]
                category=architecture_class(*signatures)
                conservative=[int(bool(hits) and all(h['conservative_architecture'] for h in hits)) for hits in [x,y]]
                matches=single_copy_matches(x,y)
                writer.writerow(dict(pair_key=pair,scope=label,policy=policy,annotation_class=category,left_hits=len(x),right_hits=len(y),left_conservative=conservative[0],right_conservative=conservative[1],both_conservative=int(all(conservative)),left_domain_hits=sum(h['pfam_type']=='Domain' for h in x),right_domain_hits=sum(h['pfam_type']=='Domain' for h in y),single_copy_candidate_domain_matches=len(matches)))
                counts[policy+':'+category]+=1;matched[policy]+=len(matches)
                for accession,l,r in matches:dw.writerow([pair,label,policy,accession,*left,l['hit_id'],l['alignment_start'],l['alignment_end'],*right,r['hit_id'],r['alignment_start'],r['alignment_end']])
    db.close();verify()
    result=dict(status='complete_background_domain_control_inventory_pending_readback',plan_sha256=ph,models=dict(model_counts),unique_model_pairs=len(pairs),pair_policy_rows=len(pairs)*len(POLICIES),pair_scope_counts=dict(Counter('+'.join(sorted(s)) for s in scope.values())),architecture_counts=dict(counts),single_copy_domain_matches_by_policy=dict(matched),artifacts={p.name:sha(p) for p in out.iterdir()},scope='All verified qualified distinct background model pairs and all candidate models joined by exact model/version/sequence to the audited Pfam registry. Four policy alternatives retained. Ordered annotation signatures include all Pfam types and repeat occurrences. Unannotated is not absent. Single-copy Domain interval candidates require unique retained accession occurrence in both models and conservative candidate flags; no residue confidence, PAE, domain alignment, ancestral gain/loss/rearrangement or biological duplication inference.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
