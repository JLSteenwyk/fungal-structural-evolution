#!/usr/bin/env python3
"""Verify the complete background_domain order summary against audited source tables."""
import csv,json,math,hashlib
from pathlib import Path
from collections import Counter


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rows(p):
    with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    pp=Path('metadata/background_domain_usable_orders_plan_20260927.json');plan=json.loads(pp.read_text());out=Path(plan['output'])
    r=json.loads((out/'receipt.json').read_text())
    assert r['status']=='complete_background_domain_numerically_usable_order_summary_pending_independent_readback' and r['plan_sha256']==sha(pp)
    for p,h in plan['pins'].items():assert sha(p)==h
    for name,h in r['artifacts'].items():assert sha(out/name)==h
    for field in ['native','diagnostic','geometry']:
        assert sha(Path(plan[field])/'receipt.json')==r[field+'_receipt_sha256']
    assert sha(plan['geometry_audit'])==r['geometry_audit_sha256']
    audit=json.loads(Path(plan['geometry_audit']).read_text())
    assert audit['status']=='passed_full_background_domain_geometry_readback' and audit['alignments_checked']==267246
    sources={}
    for field in ['native','diagnostic','geometry']:
        folder=Path(plan[field]);receipt=json.loads((folder/'receipt.json').read_text())
        for name,h in receipt['artifacts'].items():assert sha(folder/name)==h
        sources[field]=receipt
    native={}
    for x in rows(Path(plan['native'])/'checkpoint_manifest.tsv'):
        pair,mask,order=Path(x['path']).stem.split('-');key=pair,mask,int(order)
        assert key not in native;native[key]=x['status']
    tables=[]
    for field,name in [('diagnostic','numeric_readback.tsv'),('geometry','alignment_geometry.tsv')]:
        raw=rows(Path(plan[field])/name);index={(x['pair_key'],x['mask'],int(x['order'])):x for x in raw}
        assert len(index)==len(raw)==267246 and set(index)=={k for k,v in native.items() if v=='aligned'}
        tables.append(index)
    numeric,geometry=tables
    data=rows(out/'pair_mask_order_summary.tsv');seen=set();counts=Counter();exclusions=Counter();directions=0;metrics_checked=0
    metric_names=['aligned_length','rmsd_recomputed','sequence_identity_exact','tm_a','tm_b','coverage_a','coverage_b','joint_plddt70_fraction']
    for row in data:
        pair,mask=row['pair_key'],row['mask'];assert (pair,mask) not in seen;seen.add((pair,mask));values={}
        for order in [0,1]:
            key=pair,mask,order;raw=native[key];assert row[f'order{order}_native_status']==raw
            excluded=[];expected={}
            if raw=='aligned':
                n=numeric[key];g=geometry[key]
                assert row[f'order{order}_rmsd_status']==n['rmsd_status'] and row[f'order{order}_geometry_status']==g['geometry_status']
                if float(n['rmsd_rounding_error'])>.00501:excluded.append('rmsd_discrepancy')
                if int(n['aligned_length'])<3:excluded.append('fewer_than_three_pairs')
                if g['geometry_status']!='unique_at_numeric_tolerance':excluded.append('nonunique_rotation')
                if not excluded:
                    directions+=1
                    expected={x:float(n[x]) for x in ['aligned_length','rmsd_recomputed','sequence_identity_exact','joint_plddt70_fraction']}
                    for endpoint in ['a','b']:
                        source=('left' if endpoint=='a' else 'right') if order==0 else ('right' if endpoint=='a' else 'left')
                        expected['tm_'+endpoint]=float(n['tm_'+source+'_native']);expected['coverage_'+endpoint]=float(n['coverage_'+source])
                    values[order]=expected
            else:
                assert row[f'order{order}_rmsd_status']==row[f'order{order}_geometry_status']==''
            why=';'.join(excluded);assert row[f'order{order}_numerical_exclusion_reasons']==why
            assert row[f'order{order}_status']==('excluded_numerically' if excluded else raw)
            if excluded:exclusions[why]+=1
            for metric in metric_names:
                actual=row[f'order{order}_{metric}']
                if expected:assert math.isfinite(float(actual)) and float(actual)==expected[metric];metrics_checked+=1
                else:assert actual==''
        status={0:'neither_order_numerically_usable',1:'one_order_numerically_usable',2:'both_orders_numerically_usable'}[len(values)]
        assert row['order_summary_status']==status
        counts[mask+':'+status]+=1
        for metric in metric_names:
            actual=row[metric+'_order_absolute_difference']
            if len(values)==2:assert float(actual)==abs(values[0][metric]-values[1][metric]);metrics_checked+=1
            else:assert actual==''
    pairs={k[0] for k in native}
    assert len(native)==267716 and len(pairs)==66929
    assert seen=={(p,m) for p in pairs for m in ['full','plddt70']}
    assert len(data)==r['pair_mask_rows']==133858 and r['pairs']==len(pairs)
    assert dict(counts)==r['counts'] and dict(exclusions)==r['exclusion_reason_combinations']
    result=dict(status='passed_full_background_domain_usable_order_summary_readback',pair_mask_rows=len(data),numerically_usable_directions=directions,excluded_directions=sum(exclusions.values()),numeric_values_checked=metrics_checked,counts=dict(counts),exclusions=dict(exclusions),source_receipt_sha256=sha(out/'receipt.json'),geometry_audit_sha256=sha(plan['geometry_audit']),script_sha256=sha(Path(__file__)),scope='All pair/mask/order identities, native statuses, numerical exclusions, blank excluded metrics, endpoint-normalized values and order differences independently reconstructed from audited source tables. No new coordinate fits, coverage qualification or scientific acceptance.')
    with Path('metadata/background_domain_usable_orders_readback_20260927.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
