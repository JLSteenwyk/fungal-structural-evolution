#!/usr/bin/env python3
"""Join all duplication/reference models and pairs to policy-preserving domain annotations."""
import argparse,csv,hashlib,json,sqlite3
from collections import Counter,defaultdict
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha
from validate_duplication_reference_coordinates import load_additional_models,model_map

POLICIES=('alignment_evalue','alignment_bitscore','envelope_evalue','envelope_bitscore')


def architecture_class(left,right):
    if not left and not right:return 'neither_annotated'
    if not left or not right:return 'one_unannotated'
    if left==right:return 'same_ordered_annotations'
    if Counter(left)==Counter(right):return 'same_content_different_order'
    return 'different_annotation_content'


def single_copy_matches(left,right):
    # Restrict to conservative eligible Domain hits, but count all retained
    # Domain occurrences: an ineligible second copy must not make a repeat unique.
    a=defaultdict(list);b=defaultdict(list)
    for hits,target in [(left,a),(right,b)]:
        for h in hits:
            if h['pfam_type']=='Domain':target[h['pfam_accession']].append(h)
    return [(accession,a[accession][0],b[accession][0]) for accession in sorted(a.keys()&b.keys())
            if len(a[accession])==len(b[accession])==1
            and a[accession][0]['domain_interval_candidate'] and b[accession][0]['domain_interval_candidate']]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        if sha(a.plan)!=ph:raise ValueError('Changed plan')
        for path,h in plan['pins'].items():
            if sha(path)!=h:raise ValueError('Changed source: '+path)
    verify();registry=Path(plan['registry']);rr=json.loads((registry/'receipt.json').read_text());audit=json.loads(Path(plan['registry_readback']).read_text());dbpath=registry/'structure_domains.sqlite'
    if audit['status']!='passed_full_structure_domain_registry_readback' or audit['producer_receipt_sha256']!=sha(registry/'receipt.json') or audit['database_sha256']!=rr['database_sha256'] or sha(dbpath)!=rr['database_sha256']:raise ValueError('Unbound domain registry')
    base=Path(plan['base_queue']);inventory=Path(plan['inventory']);models=model_map(base/'models.jsonl')
    extras=load_additional_models(plan)
    for m in extras:
        key=m['model_id'],m['version']
        if key in models:raise ValueError('Additional model already primary')
        models[key]=m
    pairs={};scope=defaultdict(set)
    for label,folder in [('primary',base),('reference',inventory)]:
        receipt=json.loads((folder/'receipt.json').read_text());path=folder/'model_pairs.tsv'
        if sha(path)!=receipt['artifacts']['model_pairs.tsv']:raise ValueError('Changed pair artifact')
        with path.open() as handle:
            for r in csv.DictReader(handle,delimiter='\t'):
                ends=tuple(sorted([(r['model_a'],int(r['version_a'])),(r['model_b'],int(r['version_b']))]))
                key=hashlib.sha256(json.dumps(ends,separators=(',',':')).encode()).hexdigest()
                if key!=r['pair_key'] or ends[0]==ends[1] or any(e not in models for e in ends):raise ValueError('Invalid pair')
                if label in scope[key] or (key in pairs and pairs[key]!=ends):raise ValueError('Repeated pair')
                pairs[key]=ends;scope[key].add(label)
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
    result=dict(status='complete_duplication_domain_control_inventory_pending_readback',plan_sha256=ph,models=dict(model_counts),unique_model_pairs=len(pairs),pair_policy_rows=len(pairs)*len(POLICIES),pair_scope_counts=dict(Counter('+'.join(sorted(s)) for s in scope.values())),architecture_counts=dict(counts),single_copy_domain_matches_by_policy=dict(matched),artifacts={p.name:sha(p) for p in out.iterdir()},scope='All primary/reference model pairs joined by exact model/version/sequence to the audited Pfam registry. Four policy alternatives retained. Ordered annotation signatures include all Pfam types and repeat occurrences. Unannotated is not absent. Single-copy Domain interval candidates require unique retained accession occurrence in both models and conservative candidate flags; no residue confidence, PAE, domain alignment, ancestral gain/loss/rearrangement or biological duplication inference.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
