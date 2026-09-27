#!/usr/bin/env python3
"""Project audited reference numerical availability onto every event/reference link."""
import csv,json,hashlib
from collections import Counter
from pathlib import Path


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rows(p):
    with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    inventory=Path('results/structural_comparisons/duplication-reference-comparison-inventory-20260926-v1')
    summary=Path('results/structural_comparisons/duplication-reference-usable-orders-20260927-v1')
    ir=json.loads((inventory/'receipt.json').read_text());sr=json.loads((summary/'receipt.json').read_text())
    ip=Path('results/structural_comparisons/duplication-reference-comparison-readback-20260926-v1.json');sp=Path('metadata/duplication_reference_usable_orders_readback_20260927.json')
    assert json.loads(ip.read_text())['producer_receipt_sha256']==sha(inventory/'receipt.json')
    assert json.loads(sp.read_text())['source_receipt_sha256']==sha(summary/'receipt.json')
    for folder,receipt in [(inventory,ir),(summary,sr)]:
        for name,h in receipt['artifacts'].items():assert sha(folder/name)==h
    lookup={(x['pair_key'],x['mask']):x for x in rows(summary/'pair_mask_order_summary.tsv')}
    source=rows(inventory/'event_reference_comparisons.tsv');assert len(source)==73888
    links=[];groups={};counts=Counter()
    identity=['guide','family','gene_node','gene_a','gene_b','reference_gene']
    for row in source:
        assert row['focal_side'] in ['a','b'] and row['focal_gene']==row['gene_'+row['focal_side']]
        for mask in ['full','plddt70']:
            work=row['work_disposition'];number=''
            if work=='additional_pair':
                record=lookup[row['pair_key'],mask];number=sum(record[f'order{o}_status']=='aligned' for o in [0,1]);status=f'{number}_directions_numerically_usable'
            elif work=='existing_duplicate_pair':status='awaiting_primary_pair_audit'
            else:
                assert work=='identical_model' and (row['focal_model'],row['focal_version'])==(row['reference_model'],row['reference_version']);status='identical_model_not_independent_comparison'
            link=dict(row,mask=mask,numerical_availability=status,usable_directions=number);links.append(link)
            key=tuple(row[x] for x in identity)+(mask,);group=groups.setdefault(key,{})
            assert row['focal_side'] not in group;group[row['focal_side']]=link
            counts[row['guide']+':'+mask+':'+status]+=1
    events=[];combinations=Counter()
    for key,sides in sorted(groups.items()):
        assert set(sides)=={'a','b'}
        a,b=sides['a'],sides['b'];assert a['lexical_representative']==b['lexical_representative']
        event=dict(zip(identity+['mask'],key),lexical_representative=a['lexical_representative'],side_a_status=a['numerical_availability'],side_b_status=b['numerical_availability'],both_sides_all_orders_numerically_usable=a['usable_directions']==b['usable_directions']==2)
        events.append(event);combinations[(event['guide'],event['mask'],event['side_a_status'],event['side_b_status'])]+=1
    out=Path('results/structural_comparisons/duplication-reference-event-coverage-20260927-v1');out.mkdir(parents=True,exist_ok=False)
    combined=[dict(guide=k[0],mask=k[1],side_a_status=k[2],side_b_status=k[3],event_reference_rows=v) for k,v in sorted(combinations.items())]
    for name,data in [('comparison_link_coverage.tsv',links),('event_reference_coverage.tsv',events),('coverage_summary.tsv',combined)]:
        with (out/name).open('w') as f:
            w=csv.DictWriter(f,list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    r=dict(status='complete_full_reference_event_coverage_pending_independent_readback',source_links=len(source),masked_links=len(links),masked_event_references=len(events),link_status_counts=dict(counts),both_sides_all_orders=sum(x['both_sides_all_orders_numerically_usable'] for x in events),inventory_receipt_sha256=sha(inventory/'receipt.json'),summary_receipt_sha256=sha(summary/'receipt.json'),inventory_readback_sha256=sha(ip),summary_readback_sha256=sha(sp),script_sha256=sha(Path(__file__)),artifacts={p.name:sha(p) for p in out.iterdir()},scope='All provisionally eligible event/tied-reference/duplicate-side links retained for both masks. Completed additional comparisons projected; reused primary pairs remain pending audit and identical-model links remain explicit without invented metrics. Numerical availability is not orthology, coverage/confidence suitability or structural asymmetry. This inventory excludes original events lacking a provisional reference, whose dispositions remain in the upstream full event ledger.')
    (out/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))


if __name__=='__main__':main()
