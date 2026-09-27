#!/usr/bin/env python3
"""Independently verify all reference-event coverage links with dataframe joins."""
import json,hashlib
from pathlib import Path
import pandas as pd


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read(p):return pd.read_csv(p,sep='\t',dtype=str,keep_default_na=False)


def main():
    root=Path('results/structural_comparisons/duplication-reference-event-coverage-20260927-v1');r=json.loads((root/'receipt.json').read_text())
    inventory=Path('results/structural_comparisons/duplication-reference-comparison-inventory-20260926-v1');summary=Path('results/structural_comparisons/duplication-reference-usable-orders-20260927-v1')
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    for folder,key in [(inventory,'inventory_receipt_sha256'),(summary,'summary_receipt_sha256')]:
        assert sha(folder/'receipt.json')==r[key]
        for name,h in json.loads((folder/'receipt.json').read_text())['artifacts'].items():assert sha(folder/name)==h
    source=read(inventory/'event_reference_comparisons.tsv');expected=source.merge(pd.DataFrame({'mask':['full','plddt70']}),how='cross')
    metrics=read(summary/'pair_mask_order_summary.tsv')
    lookup=metrics[['pair_key','mask','order0_status','order1_status']].copy()
    lookup['expected_directions']=((lookup.order0_status=='aligned').astype(int)+(lookup.order1_status=='aligned').astype(int)).astype(str)
    expected=expected.merge(lookup[['pair_key','mask','expected_directions']],on=['pair_key','mask'],how='left',validate='many_to_one')
    additional=expected.work_disposition=='additional_pair';assert expected.loc[additional,'expected_directions'].notna().all()
    expected['usable_directions']=expected.expected_directions.where(additional,'')
    expected['numerical_availability']=''
    expected.loc[additional,'numerical_availability']=expected.loc[additional,'usable_directions']+'_directions_numerically_usable'
    expected.loc[expected.work_disposition=='existing_duplicate_pair','numerical_availability']='awaiting_primary_pair_audit'
    expected.loc[expected.work_disposition=='identical_model','numerical_availability']='identical_model_not_independent_comparison'
    assert (expected.numerical_availability!='').all()
    links=read(root/'comparison_link_coverage.tsv');keys=['guide','family','gene_node','gene_a','gene_b','reference_gene','focal_side','mask']
    assert not links.duplicated(keys).any() and not expected.duplicated(keys).any()
    left=links.sort_values(keys).reset_index(drop=True);right=expected.sort_values(keys).reset_index(drop=True)
    pd.testing.assert_frame_equal(left,right[left.columns],check_dtype=False)
    eventkeys=[x for x in keys if x!='focal_side'];events=read(root/'event_reference_coverage.tsv')
    assert len(expected.groupby(eventkeys))==len(events) and not events.duplicated(eventkeys).any()
    a=expected[expected.focal_side=='a'];b=expected[expected.focal_side=='b']
    joined=a.merge(b,on=eventkeys,suffixes=('_a','_b'),validate='one_to_one')
    candidate=joined[eventkeys].copy();candidate['lexical_representative']=joined.lexical_representative_a
    assert (joined.lexical_representative_a==joined.lexical_representative_b).all()
    candidate['side_a_status']=joined.numerical_availability_a;candidate['side_b_status']=joined.numerical_availability_b
    candidate['both_sides_all_orders_numerically_usable']=((joined.usable_directions_a=='2')&(joined.usable_directions_b=='2')).map(str)
    pd.testing.assert_frame_equal(events.sort_values(eventkeys).reset_index(drop=True),candidate[events.columns].sort_values(eventkeys).reset_index(drop=True),check_dtype=False)
    summarykeys=['guide','mask','side_a_status','side_b_status']
    grouped=candidate.groupby(summarykeys).size().rename('event_reference_rows').reset_index();grouped.event_reference_rows=grouped.event_reference_rows.astype(str)
    observed=read(root/'coverage_summary.tsv')
    pd.testing.assert_frame_equal(observed.sort_values(summarykeys).reset_index(drop=True),grouped.sort_values(summarykeys).reset_index(drop=True),check_dtype=False)
    assert len(source)==r['source_links'] and len(links)==r['masked_links'] and len(events)==r['masked_event_references']
    count=int((candidate.both_sides_all_orders_numerically_usable=='True').sum());assert count==r['both_sides_all_orders']
    result=dict(status='passed_full_reference_event_coverage_independent_join',source_links=len(source),masked_links=len(links),masked_event_references=len(events),summary_rows=len(observed),both_sides_all_orders=count,source_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(Path(__file__)),scope='All original link fields, two-mask expansion, many-to-one numerical-status joins, both duplicate-side pivots and all coverage summaries independently reconstructed. Identical-model and unfinished reused comparisons explicit. No structural asymmetry or biological suitability claim.')
    with Path('metadata/duplication_reference_event_coverage_readback_20260927.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
