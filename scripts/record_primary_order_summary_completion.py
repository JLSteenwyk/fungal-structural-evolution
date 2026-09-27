#!/usr/bin/env python3
"""Record audited primary pair eligibility without hiding missing/one-order cases."""
import csv,json,subprocess
from collections import Counter
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha


def main():
    pp=Path('metadata/primary_usable_orders_plan_20260927.json');plan=json.loads(pp.read_text())
    paths=[str(pp),'metadata/primary_usable_orders_launch_20260927.json','metadata/primary_usable_orders_readback_launch_20260927.json','metadata/primary_usable_orders_readback_20260927.json']
    states={}
    for path in paths[1:3]:
        launch=json.loads(Path(path).read_text());assert launch['plan_sha256']==sha(pp)
        unit=launch['unit'];state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',unit,'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
        assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0'),state
        states[unit]=state
    for path,h in plan['pins'].items():assert sha(path)==h,path
    root=Path(plan['output']);rp=root/'receipt.json';receipt=json.loads(rp.read_text());audit=json.loads(Path(paths[-1]).read_text())
    assert receipt['status']=='complete_primary_numerically_usable_order_summary_pending_independent_readback'
    assert audit['status']=='passed_full_primary_usable_order_summary_readback'
    assert receipt['plan_sha256']==audit['plan_sha256']==sha(pp)
    assert audit['source_receipt_sha256']==sha(rp)
    assert audit['geometry_audit_sha256']==receipt['geometry_audit_sha256']==sha(plan['geometry_audit'])
    assert receipt['counts']==audit['counts'] and receipt['exclusion_reason_combinations']==audit['exclusions']
    for name,h in receipt['artifacts'].items():assert sha(root/name)==h
    counts=Counter();dispositions=Counter();summary=Counter();seen=set();numerical_exclusions=0;available=0;missing=0
    path=root/'pair_mask_order_summary.tsv'
    with path.open() as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            key=row['pair_key'],row['mask'];assert key not in seen;seen.add(key)
            statuses=[row[f'order{i}_status'] for i in [0,1]]
            counts[row['mask']+':'+row['order_summary_status']]+=1
            dispositions[(row['mask'],*statuses)]+=1
            available+=statuses.count('aligned');missing+=statuses.count('input_unavailable');numerical_exclusions+=statuses.count('excluded_numerically')
            if all(v=='aligned' for v in statuses):category='both_orders_numerically_usable'
            elif all(v=='input_unavailable' for v in statuses):category='both_orders_input_unavailable'
            elif all(v=='excluded_numerically' for v in statuses):category='both_orders_excluded_numerically'
            elif statuses.count('aligned')==1:category='one_order_numerically_usable'
            else:category='other_retained_disposition'
            summary[row['mask']+':'+category]+=1
    assert len(seen)==receipt['pair_mask_rows']==audit['pair_mask_rows']==206400
    assert dict(counts)==receipt['counts']
    assert available==audit['numerically_usable_directions']==387319
    assert numerical_exclusions==audit['excluded_directions']==327
    assert available+missing+numerical_exclusions==412800
    assert missing==25154
    paths += [str(rp),str(path),plan['geometry_audit']]
    result=dict(status='complete_verified_primary_order_summary',pair_mask_rows=len(seen),distinct_pairs=103200,numerically_usable_directions=available,missing_input_directions=missing,numerically_excluded_directions=numerical_exclusions,pair_categories=dict(summary),ordered_dispositions=[dict(mask=k[0],order0_status=k[1],order1_status=k[2],pairs=v) for k,v in sorted(dispositions.items())],numeric_values_independently_checked=audit['numeric_values_checked'],terminal_states=states,source_hashes={p:sha(p) for p in paths},script_sha256=sha(__file__),scope='Complete audited eligibility ledger, including input-unavailable and asymmetric order qualification. Pair counts are numerical disposition counts, not independent events, confidence/coverage-qualified comparisons or biological effects. Full-mask and confidence-mask results are dependent.')
    with Path('metadata/primary_order_summary_completed_20260927.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
