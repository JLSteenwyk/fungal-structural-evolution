#!/usr/bin/env python3
"""Report complete order/weighting sensitivity after the full record-summary audit."""
import json,csv,hashlib,time,subprocess
from pathlib import Path
from collections import Counter
import psutil,pandas as pd,numpy as np

sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    launch=Path('metadata/matched_domain_record_summary_readback_launch_20260927.json');dep=json.loads(launch.read_text());lh=sha(launch)
    while True:
        try:
            p=psutil.Process(dep['pid'])
            if p.create_time()!=dep['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==dep['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    assert sha(launch)==lh
    state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',dep['unit'],'-p','ActiveState','-p','ExecMainStatus'],text=True).splitlines());assert state=={'ActiveState':'inactive','ExecMainStatus':'0'}
    root=Path('results/structural_comparisons/matched-domain-record-summaries-20260927-v1');proof=Path('metadata/matched_domain_record_summary_readback_20260927.json');a=json.loads(proof.read_text());r=json.loads((root/'receipt.json').read_text())
    assert a['status']=='passed_full_matched_domain_record_summary_readback' and a['source_receipt_sha256']==sha(root/'receipt.json')
    source=root/'record_summary.tsv';assert sha(source)==r['artifacts']['record_summary.tsv'];data=pd.read_csv(source,sep='\t')
    base=['guide','policy','scenario_id','boundary','mask','cohort','screen'];orders=['target_order','background_order'];weights=['record','family_equal','taxon_equal'];metrics=['rmsd_difference','identity_difference','original_coverage_difference','log_aligned_length_ratio','confidence_fraction_difference'];tol=1e-10
    assert len(data)==82944 and not data.duplicated(base+orders).any();out=Path('results/structural_comparisons/matched-record-sensitivity-20260927-v1');out.mkdir(parents=True,exist_ok=False)
    rows=[];counts=Counter();weighting_flips=0;checked_values=0
    for key,group in data.groupby(base,sort=True):
        assert len(group)==4 and set(zip(group.target_order,group.background_order))=={(0,0),(0,1),(1,0),(1,1)}
        assert group.matched_records.nunique()==1;count=int(group.matched_records.iloc[0]);allvalues=group[[f'rmsd_difference_{w}_mean' for w in weights]].to_numpy()
        flip=bool(np.any((allvalues.min(axis=1)<-tol)&(allvalues.max(axis=1)>tol))) if count else False;weighting_flips+=flip
        for weight in weights:
            row=dict(zip(base,key));row.update(weighting=weight,matched_records=count,weighting_sign_flip_within_any_order=int(flip))
            vals=group[f'rmsd_difference_{weight}_mean'].to_numpy()
            status='no_matched_records' if not count else 'positive_all_orders' if np.all(vals>tol) else 'negative_all_orders' if np.all(vals < -tol) else 'within_numeric_zero_all_orders' if np.all(abs(vals)<=tol) else 'mixed_sign_or_numeric_zero'
            row['rmsd_order_sign_status']=status;counts[weight+':'+status]+=1
            for metric in metrics:
                values=group[f'{metric}_{weight}_mean'];lo=values.min();hi=values.max()
                if count:
                    assert np.isfinite(values).all();ordered=sorted(float(v) for v in values);assert lo==ordered[0] and hi==ordered[-1]
                    row[metric+'_order_min']=lo;row[metric+'_order_max']=hi;row[metric+'_order_range']=hi-lo;checked_values+=8
                else:
                    assert values.isna().all();row[metric+'_order_min']=row[metric+'_order_max']=row[metric+'_order_range']=''
            rows.append(row)
    assert len(rows)==62208
    with (out/'order_weighting_sensitivity.tsv').open('w') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    assert sha(source)==r['artifacts']['record_summary.tsv'] and a['source_receipt_sha256']==sha(root/'receipt.json')
    result=dict(status='complete_full_matched_record_sensitivity',source_summary_rows=len(data),order_groups=len(data)//4,weighting_rows=len(rows),sign_numeric_tolerance=tol,rmsd_sign_counts=dict(counts),groups_with_weighting_sign_flip_within_order=weighting_flips,range_endpoint_source_values_checked=checked_values,source_receipt_sha256=sha(root/'receipt.json'),source_readback_sha256=sha(proof),script_sha256=sha(__file__),artifacts={'order_weighting_sensitivity.tsv':sha(out/'order_weighting_sensitivity.tsv')},scope='Every completed record-summary cell included. All four order identities and fixed cohort counts checked; ranges independently cross-checked by sorted scalar endpoints. Sign categories use explicit 1e-10 numerical tolerance. Order ranges are sensitivity ranges, not uncertainty intervals; positive contrast is not a significant or causal duplication effect. No preferred setting chosen.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
