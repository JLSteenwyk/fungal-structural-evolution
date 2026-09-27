#!/usr/bin/env python3
"""Compare exact duplicate-pair eligibility and reference choices between guides."""
import argparse,csv,json
from collections import Counter
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha


def reference_relation(a,b):
    if a is None:return 'mafft_only_candidate'
    if b is None:return 'profile_only_candidate'
    if a['status']!='provisional_reference_available' or b['status']!='provisional_reference_available':return 'not_provisional_in_both'
    x=set(json.loads(a['nearest_reference_genes']));y=set(json.loads(b['nearest_reference_genes']))
    if not x or not y:raise ValueError('Provisional reference without candidate')
    if x==y:return 'same_nearest_reference_set'
    if x&y:return 'overlapping_nearest_reference_sets'
    return 'disjoint_nearest_reference_sets'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['inventory','plan','review','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    rp=a.inventory/'receipt.json';r=json.loads(rp.read_text());plan=json.loads(a.plan.read_text());rr=json.loads((a.review/'receipt.json').read_text())
    if r['status']!='complete_duplication_sister_reference_inventory' or r['plan_sha256']!=sha(a.plan):raise ValueError('Incomplete/mismatched inventory')
    if Path(plan['review']).resolve()!=a.review.resolve() or plan['pins'][str(a.review/'receipt.json')]!=sha(a.review/'receipt.json'):raise ValueError('Candidate review differs')
    source_pins={str(rp):sha(rp),str(a.plan):sha(a.plan),str(a.review/'receipt.json'):sha(a.review/'receipt.json')};tables={}
    for g in ['profile','mafft']:
        path=a.inventory/(g+'_sister_references.tsv');cp=a.review/(g+'_candidate_tree_checks.tsv')
        if sha(path)!=r['artifacts'][path.name] or sha(cp)!=rr['artifacts'][cp.name]:raise ValueError('Changed input table')
        source_pins[str(path)]=sha(path);source_pins[str(cp)]=sha(cp)
        with cp.open() as f:expected={tuple(sorted([x['gene_a'],x['gene_b']])):x for x in csv.DictReader(f,delimiter='\t')}
        records={}
        with path.open() as f:
            for x in csv.DictReader(f,delimiter='\t'):
                key=tuple(sorted([x['gene_a'],x['gene_b']]))
                if key in records or key not in expected or any(x[k]!=v for k,v in expected[key].items()):raise ValueError('Original candidate fields differ')
                ties=json.loads(x['nearest_reference_genes'])
                if ties!=sorted(set(ties)) or (ties and x['chosen_reference_gene']!=ties[0]) or (not ties and x['chosen_reference_gene']):raise ValueError('Inconsistent tie representation')
                records[key]=x
        if set(records)!=set(expected):raise ValueError('Incomplete candidate universe')
        tables[g]=records
    a.output.mkdir(parents=True);summary=Counter();matrix=Counter();shared=0;both=0;same_choice=0;same_model=0
    with (a.output/'pair_reference_comparison.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['gene_a','gene_b','profile_status','mafft_status','relation','profile_reference','mafft_reference','same_chosen_reference','same_chosen_model'])
        for key in sorted(set(tables['profile'])|set(tables['mafft'])):
            x=tables['profile'].get(key);y=tables['mafft'].get(key);relation=reference_relation(x,y);summary[relation]+=1
            sx=x['status'] if x else 'not_candidate';sy=y['status'] if y else 'not_candidate';matrix[sx,sy]+=1
            eligible=sx==sy=='provisional_reference_available';choice=model=''
            if x is not None and y is not None:shared+=1
            if eligible:
                both+=1;choice=int(x['chosen_reference_gene']==y['chosen_reference_gene']);model=int((x['reference_model'],x['reference_version'])==(y['reference_model'],y['reference_version']));same_choice+=choice;same_model+=model
                if choice and not model:raise ValueError('Same reference gene has different frozen model')
            w.writerow([*key,sx,sy,relation,x['chosen_reference_gene'] if x else '',y['chosen_reference_gene'] if y else '',choice,model])
    if shared!=rr['overlap']['shared_pairs']:raise ValueError('Shared pair count differs from tree review')
    with (a.output/'eligibility_matrix.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['profile_status','mafft_status','pairs']);w.writerows((*key,n) for key,n in sorted(matrix.items()))
    for path,h in source_pins.items():
        if sha(path)!=h:raise ValueError('Source changed during comparison')
    result=dict(status='complete_duplication_reference_guide_comparison',source_pins=source_pins,shared_duplicate_pairs=shared,provisional_in_both=both,same_chosen_reference_gene=same_choice,same_chosen_reference_model=same_model,relations=dict(summary),script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in a.output.iterdir()},scope='Every original candidate field and full pair universe rechecked; eligibility and nearest-reference sets compared by protein identities, preserving ties. Does not independently reconstruct sister-clade selection or path distances. Agreement is a guide sensitivity result, not validation of reference orthology, ancestral state, structure quality or an asymmetry effect.')
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
