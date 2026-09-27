#!/usr/bin/env python3
"""Inventory all oriented whole-protein duplicate/reference triads without selecting references."""
import csv,json,hashlib
from pathlib import Path
from collections import Counter


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rows(p):
    with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))


def digest(value):return hashlib.sha256(json.dumps(value,separators=(',',':')).encode()).hexdigest()


def main():
    inventory=Path('results/structural_comparisons/duplication-reference-comparison-inventory-20260926-v1')
    base=Path('results/structural_comparisons/duplication-model-pair-queue-20260926-v1')
    ir=json.loads((inventory/'receipt.json').read_text());br=json.loads((base/'receipt.json').read_text())
    proof=Path('results/structural_comparisons/duplication-reference-comparison-readback-20260926-v1.json')
    assert json.loads(proof.read_text())['producer_receipt_sha256']==sha(inventory/'receipt.json')
    for folder,receipt,names in [(inventory,ir,['event_reference_comparisons.tsv','model_pairs.tsv','models.jsonl']),(base,br,['model_pairs.tsv'])]:
        for name in names:assert sha(folder/name)==receipt['artifacts'][name]
    primary={x['pair_key'] for x in rows(base/'model_pairs.tsv')}
    reference={x['pair_key']:x['work_disposition'] for x in rows(inventory/'model_pairs.tsv')}
    models={}
    with (inventory/'models.jsonl').open() as f:
        for line in f:
            m=json.loads(line);key=(m['model_id'],int(m['version']));assert key not in models;models[key]=m
    groups={};identity=['guide','family','gene_node','gene_a','gene_b','reference_gene']
    source=rows(inventory/'event_reference_comparisons.tsv')
    for r in source:
        key=tuple(r[x] for x in identity);sides=groups.setdefault(key,{})
        assert r['focal_side'] not in sides;sides[r['focal_side']]=r
    triads={};links=[]
    for key,sides in sorted(groups.items()):
        assert set(sides)=={'a','b'};a,b=sides['a'],sides['b']
        assert (a['reference_model'],a['reference_version'],a['lexical_representative'])==(b['reference_model'],b['reference_version'],b['lexical_representative'])
        ids=[(a['focal_model'],int(a['focal_version'])),(b['focal_model'],int(b['focal_version'])),(a['reference_model'],int(a['reference_version']))]
        assert all(x in models for x in ids)
        tid=digest(ids);record=dict(triad_id=tid)
        for role,model in zip(['a','b','reference'],ids):record.update({role+'_model':model[0],role+'_version':model[1],role+'_sequence_sha256':models[model]['sequence_sha256']})
        for label,(i,j) in {'ab':(0,1),'a_reference':(0,2),'b_reference':(1,2)}.items():
            pair=digest(sorted([ids[i],ids[j]]))
            if ids[i]==ids[j]:status='identical_model'
            elif pair in primary:status='primary_pair'
            else:assert pair in reference and reference[pair]=='additional_pair';status='additional_reference_pair'
            if label=='ab':assert status in ['primary_pair','identical_model']
            if label!='ab':
                source_row=a if i==0 else b
                assert pair==source_row['pair_key']
                assert {'identical_model':'identical_model','primary_pair':'existing_duplicate_pair','additional_reference_pair':'additional_pair'}[status]==source_row['work_disposition']
            record[label+'_pair_key']=pair;record[label+'_source']=status
        record['distinct_models']=len(set(ids))
        if tid in triads:assert triads[tid]==record
        else:triads[tid]=record
        links.append(dict(zip(identity,key),lexical_representative=a['lexical_representative'],triad_id=tid))
    assert len(links)*2==len(source)==73888
    out=Path('results/structural_comparisons/whole-protein-reference-triads-20260927-v1');out.mkdir(parents=True,exist_ok=False)
    for name,data in [('model_triads.tsv',[triads[k] for k in sorted(triads)]),('event_reference_triads.tsv',links)]:
        with (out/name).open('w') as f:
            w=csv.DictWriter(f,list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    r=dict(status='complete_full_whole_protein_reference_triad_inventory_pending_readback',event_reference_links=len(links),distinct_oriented_model_triads=len(triads),distinct_model_count_distribution=dict(Counter(x['distinct_models'] for x in triads.values())),pair_source_counts=dict(Counter(label+':'+x[label+'_source'] for x in triads.values() for label in ['ab','a_reference','b_reference'])),future_two_mask_eight_order_dispositions=16*len(triads),inventory_receipt_sha256=sha(inventory/'receipt.json'),primary_queue_receipt_sha256=sha(base/'receipt.json'),source_readback_sha256=sha(proof),script_sha256=sha(Path(__file__)),artifacts={p.name:sha(p) for p in out.iterdir()},scope='All eligible event and tied-reference links mapped to exact oriented A/B/reference model/version triples and pair-work sources. Identity models explicit without invented alignments. No residue intersection, fit, ancestral change or biological asymmetry inferred. Full original no-reference dispositions remain upstream. Future all-order mapping requires primary/reference full audits and explicit identity/mask handling.')
    (out/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))


if __name__=='__main__':main()
