#!/usr/bin/env python3
"""Preserve both audited alignment orders and expose numerical order sensitivity."""
import argparse,csv,json,math,re
from collections import Counter
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha

METRICS=['aligned_length','rmsd_recomputed','sequence_identity_exact','tm_a','tm_b','coverage_a','coverage_b','joint_plddt70_fraction']


def summarize_pair(pair,mask,statuses,numeric):
    if set(statuses)!={0,1}:raise ValueError('Missing input order')
    result=dict(pair_key=pair,mask=mask,order0_status=statuses[0],order1_status=statuses[1])
    success=[order for order in [0,1] if statuses[order]=='aligned']
    if set(numeric)!=set(success):raise ValueError('Numeric rows differ from successful dispositions')
    result['order_summary_status']='both_orders_aligned' if len(success)==2 else 'one_order_aligned' if success else 'neither_order_aligned'
    values={}
    for order in [0,1]:
        row=numeric.get(order)
        if row:
            left,right=('left','right') if order==0 else ('right','left')
            values[order]={k:float(row[k]) for k in ['aligned_length','rmsd_recomputed','sequence_identity_exact','joint_plddt70_fraction']}
            values[order].update(tm_a=float(row['tm_'+left+'_native']),tm_b=float(row['tm_'+right+'_native']),coverage_a=float(row['coverage_'+left]),coverage_b=float(row['coverage_'+right]))
            if not all(math.isfinite(v) and v>=0 for v in values[order].values()):raise ValueError('Invalid metric')
            if any(values[order][k]>1 for k in ['sequence_identity_exact','joint_plddt70_fraction','tm_a','tm_b','coverage_a','coverage_b']):raise ValueError('Invalid fraction/score')
            if values[order]['aligned_length']<1 or not values[order]['aligned_length'].is_integer():raise ValueError('Invalid aligned length')
        for metric in METRICS:result[f'order{order}_{metric}']=values[order][metric] if row else ''
    for metric in METRICS:
        result[metric+'_order_absolute_difference']=abs(values[0][metric]-values[1][metric]) if len(success)==2 else ''
    return result


def summarize(producer,readback,output):
    producer,readback,output=map(Path,(producer,readback,output))
    rp=producer/'receipt.json';ap=readback/'receipt.json'
    r=json.loads(rp.read_text());a=json.loads(ap.read_text())
    if r['status'] not in ['complete_duplication_alignment_dispositions_pending_readback','complete_reference_alignment_dispositions_pending_readback']:raise ValueError('Incomplete producer')
    if a['status']!='passed_full_duplication_alignment_mapping_rmsd_identity_readback' or a['producer_receipt_sha256']!=sha(rp):raise ValueError('Unbound numeric readback')
    cm=producer/'checkpoint_manifest.tsv';nt=readback/'numeric_readback.tsv'
    if sha(cm)!=r['artifacts']['checkpoint_manifest.tsv'] or sha(nt)!=a['artifacts']['numeric_readback.tsv']:raise ValueError('Changed table')
    pins={str(p):sha(p) for p in [rp,ap,cm,nt]}
    statuses={};counts=Counter();pairs=set()
    with cm.open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            match=re.fullmatch(r'pairs/([0-9a-f]{2})/([0-9a-f]{64})-(full|plddt70)-([01])\.json',row['path'])
            if not match:raise ValueError('Malformed checkpoint path')
            prefix,pair,mask,order=match.groups();order=int(order)
            if prefix!=pair[:2] or row['status'] not in ['aligned','input_unavailable','native_error','parse_error','timeout']:raise ValueError('Invalid disposition')
            key=pair,mask,order
            if key in statuses:raise ValueError('Duplicate disposition')
            statuses[key]=row['status'];pairs.add(pair);counts[mask+':'+row['status']]+=1
    if len(pairs)!=r['distinct_model_pairs'] or len(statuses)!=r['directed_dispositions'] or dict(counts)!=r['counts'] or dict(counts)!=a['counts'] or len(statuses)!=a['directed_dispositions']:raise ValueError('Scope/count mismatch')
    expected={(pair,mask,order) for pair in pairs for mask in ['full','plddt70'] for order in [0,1]}
    if set(statuses)!=expected:raise ValueError('Incomplete pair/mask/order grid')
    numeric={}
    with nt.open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            key=row['pair_key'],row['mask'],int(row['order'])
            if key in numeric or statuses.get(key)!='aligned':raise ValueError('Unexpected/repeated numeric row')
            numeric[key]=row
    if set(numeric)!={k for k,v in statuses.items() if v=='aligned'} or len(numeric)!=a['numerically_checked_alignments']:raise ValueError('Incomplete numeric rows')
    output.mkdir(parents=True,exist_ok=False);table=output/'pair_mask_order_summary.tsv';totals=Counter();maxima={k:0. for k in METRICS}
    with table.open('w') as f:
        writer=None
        for pair in sorted(pairs):
            for mask in ['full','plddt70']:
                row=summarize_pair(pair,mask,{o:statuses[pair,mask,o] for o in [0,1]},
                    {o:numeric[pair,mask,o] for o in [0,1] if (pair,mask,o) in numeric})
                if writer is None:writer=csv.DictWriter(f,fieldnames=list(row),delimiter='\t');writer.writeheader()
                writer.writerow(row);totals[mask+':'+row['order_summary_status']]+=1
                if row['order_summary_status']=='both_orders_aligned':
                    for k in METRICS:maxima[k]=max(maxima[k],row[k+'_order_absolute_difference'])
    for path,h in pins.items():
        if sha(path)!=h:raise ValueError('Source changed while summarizing')
    result=dict(status='complete_audited_duplication_alignment_order_summary',source_sha256=pins,script_sha256=sha(__file__),mode=a['mode'],model_pairs=len(pairs),pair_mask_rows=2*len(pairs),counts=dict(totals),maximum_order_differences=maxima,artifacts={table.name:sha(table)},scope='Both native alignment orders retained; a/b normalization refers to order-0 input endpoints. No favorable-order selection, averaging, metric additivity, pair-specific significance, structural asymmetry or biological inference. Missing/failed orders keep blank metrics. Maximum differences initialized to zero; consult successful-order counts before interpreting them.')
    (output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['producer','readback','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();print(json.dumps(summarize(a.producer,a.readback,a.output),indent=2))

if __name__=='__main__':main()
