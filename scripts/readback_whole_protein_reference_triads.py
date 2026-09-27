#!/usr/bin/env python3
"""Verify the full oriented protein-triad inventory using independent joins."""
import json,hashlib
from pathlib import Path
from collections import Counter
import pandas as pd


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read(p):return pd.read_csv(p,sep='\t',dtype=str,keep_default_na=False)


def main():
    root=Path('results/structural_comparisons/whole-protein-reference-triads-20260927-v1');r=json.loads((root/'receipt.json').read_text())
    inventory=Path('results/structural_comparisons/duplication-reference-comparison-inventory-20260926-v1');base=Path('results/structural_comparisons/duplication-model-pair-queue-20260926-v1')
    for folder,key in [(inventory,'inventory_receipt_sha256'),(base,'primary_queue_receipt_sha256')]:assert sha(folder/'receipt.json')==r[key]
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    ir=json.loads((inventory/'receipt.json').read_text());br=json.loads((base/'receipt.json').read_text())
    for folder,receipt,names in [(inventory,ir,['event_reference_comparisons.tsv','model_pairs.tsv','models.jsonl']),(base,br,['model_pairs.tsv'])]:
        for name in names:assert sha(folder/name)==receipt['artifacts'][name]
    source=read(inventory/'event_reference_comparisons.tsv');links=read(root/'event_reference_triads.tsv');triads=read(root/'model_triads.tsv')
    keys=['guide','family','gene_node','gene_a','gene_b','reference_gene']
    assert not links.duplicated(keys).any() and not triads.triad_id.duplicated().any()
    expanded=links.merge(triads,on='triad_id',validate='many_to_one')
    for side in ['a','b']:
        original=source[source.focal_side==side];joined=expanded.merge(original,on=keys,validate='one_to_one',suffixes=('','_source'))
        assert len(joined)==len(original)==len(links)
        for field,target in [('focal_model',side+'_model'),('focal_version',side+'_version'),('reference_model','reference_model'),('reference_version','reference_version'),('lexical_representative','lexical_representative')]:
            sourcefield=field+'_source' if field in expanded.columns else field
            assert (joined[sourcefield]==joined[target]).all()
        assert (joined[side+'_reference_pair_key']==joined.pair_key).all()
    assert set(links.triad_id)==set(triads.triad_id)
    with (inventory/'models.jsonl').open() as f:models={(m['model_id'],int(m['version'])):m for m in map(json.loads,f)}
    primary=set(read(base/'model_pairs.tsv').pair_key);refs=read(inventory/'model_pairs.tsv').set_index('pair_key').work_disposition.to_dict()
    count=Counter();source_counts=Counter()
    for row in triads.to_dict('records'):
        ordered=[[row[role+'_model'],int(row[role+'_version'])] for role in ['a','b','reference']]
        assert hashlib.sha256(json.dumps(ordered,separators=(',',':')).encode()).hexdigest()==row['triad_id']
        distinct=len({tuple(x) for x in ordered});assert int(row['distinct_models'])==distinct;count[str(distinct)]+=1
        for role,model in zip(['a','b','reference'],ordered):assert row[role+'_sequence_sha256']==models[tuple(model)]['sequence_sha256']
        for label,i,j in [('ab',0,1),('a_reference',0,2),('b_reference',1,2)]:
            endpoints=sorted([ordered[i],ordered[j]]);digest=hashlib.sha256(json.dumps(endpoints,separators=(',',':')).encode()).hexdigest()
            assert row[label+'_pair_key']==digest
            expected='identical_model' if endpoints[0]==endpoints[1] else 'primary_pair' if digest in primary else 'additional_reference_pair'
            if expected=='additional_reference_pair':assert refs[digest]=='additional_pair'
            assert row[label+'_source']==expected;source_counts[label+':'+expected]+=1
    assert dict(count)==r['distinct_model_count_distribution'] and dict(source_counts)==r['pair_source_counts']
    assert len(triads)==r['distinct_oriented_model_triads'] and len(links)==r['event_reference_links'] and len(source)==2*len(links)
    assert r['future_two_mask_eight_order_dispositions']==16*len(triads)
    result=dict(status='passed_full_whole_protein_reference_triad_independent_join',source_side_links=len(source),event_reference_links=len(links),distinct_oriented_model_triads=len(triads),future_mask_order_dispositions=16*len(triads),source_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(Path(__file__)),scope='Full independent source-side joins, event coverage, exact model/version/sequence hashes, triad/pair identities and work partitions checked. No residue mapping, fitting or biological acceptance.')
    with Path('metadata/whole_protein_reference_triads_readback_20260927.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
