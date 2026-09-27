#!/usr/bin/env python3
"""Join audited target/background coverage and geometry without dropping pairs."""
import argparse,csv,json,time
from pathlib import Path
from collections import Counter
import psutil
from screen_duplication_domain_alignment_coverage import sha


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args()
    plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify();dep=plan['producer']
    while True:
        try:
            p=psutil.Process(dep['pid'])
            if p.create_time()!=dep['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==dep['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    verify();out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);summaries={};bindings={}
    for label,spec in plan['cohorts'].items():
        coverage=Path(spec['coverage']);geometry=Path(spec['geometry'])
        cr=json.loads((coverage/'receipt.json').read_text());gr=json.loads((geometry/'receipt.json').read_text())
        ca=json.loads(Path(spec['coverage_audit']).read_text());ga=json.loads(Path(spec['geometry_audit']).read_text())
        assert ca['status']==spec['coverage_audit_status'] and ga['status']==spec['geometry_audit_status']
        assert ca['producer_receipt_sha256']==sha(coverage/'receipt.json') and ga['producer_receipt_sha256']==sha(geometry/'receipt.json')
        assert gr['diagnostic_receipt_sha256']==cr['diagnostic_receipt_sha256']
        for root,r in [(coverage,cr),(geometry,gr)]:
            for name,h in r['artifacts'].items():assert sha(root/name)==h
        geos={}
        for row in csv.DictReader((geometry/'alignment_geometry.tsv').open(),delimiter='\t'):
            k=row['pair_key'],row['mask'],row['order'];assert k not in geos;geos[k]=row
        seen=set();geo_seen=set();counts=Counter();pairs={};screenids=plan['screens']
        fields=['pair_key','mask','interval_a','interval_b','order0_geometry_status','order1_geometry_status','both_rotations_unique']
        for sid in screenids:fields += [sid+'_pass',sid+'_exclusions']
        with (out/(label+'.tsv')).open('w') as f:
            w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader()
            for c in csv.DictReader((coverage/'pair_mask_coverage.tsv').open(),delimiter='\t'):
                key=c['pair_key'],c['mask'];assert key not in seen;seen.add(key)
                row={k:c[k] for k in ['pair_key','mask','interval_a','interval_b']};unique=True;degenerate=False
                for order in ['0','1']:
                    k=(*key,order);g=geos.get(k);raw=c['order'+order+'_status']
                    assert (g is None)==(raw=='input_unavailable')
                    if g:
                        geo_seen.add(k);assert g['aligned_length']==c['order'+order+'_aligned_length'] and g['rmsd_status']==raw
                    status=g['geometry_status'] if g else 'input_unavailable'
                    row['order'+order+'_geometry_status']=status
                    unique &= status=='unique_at_numeric_tolerance';degenerate |= status=='degenerate_at_numeric_tolerance'
                row['both_rotations_unique']=int(unique)
                flags={}
                for sid in screenids:
                    why=c[sid+'_exclusions'].split(';') if c[sid+'_exclusions'] else []
                    assert (c[sid+'_pass']=='1')==(not why)
                    if degenerate:why.append('nonunique_rotation')
                    passed=(c[sid+'_pass']=='1') and unique
                    assert passed==(not why)
                    row[sid+'_pass']=int(passed);row[sid+'_exclusions']=';'.join(why);flags[sid]=passed
                    counts[key[1]+':'+sid]+=passed
                pairs.setdefault(key[0],{})[key[1]]=flags;w.writerow(row)
        assert geo_seen==set(geos) and len(geos)==ga['alignments_checked']==gr['alignments']
        assert len(seen)==ca['rows_checked']==cr['pair_mask_rows'] and len(pairs)==cr['distinct_pairs']
        assert all(set(v)=={'full','plddt70'} for v in pairs.values())
        summaries[label]=dict(pair_mask_rows=len(seen),geometry_rows=len(geos),pass_counts=dict(counts),both_masks_pass_counts={sid:sum(v['full'][sid] and v['plddt70'][sid] for v in pairs.values()) for sid in screenids})
        bindings[label]={str(p):sha(p) for p in [coverage/'receipt.json',geometry/'receipt.json',Path(spec['coverage_audit']),Path(spec['geometry_audit'])]}
    verify()
    for b in bindings.values():
        for p,h in b.items():assert sha(p)==h
    r=dict(status='complete_domain_coverage_geometry_join_pending_readback',plan_sha256=ph,cohorts=summaries,source_bindings=bindings,artifacts={p.name:sha(p) for p in out.glob('*.tsv')},scope='Identical coverage thresholds plus numerical rotation uniqueness for target and background domains. All pairs and exclusions retained. Not matched-event eligibility, prediction error calibration, or duplication-effect inference.')
    (out/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2),flush=True)

if __name__=='__main__':main()
