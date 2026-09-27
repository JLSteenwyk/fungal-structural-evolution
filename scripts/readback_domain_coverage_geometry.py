#!/usr/bin/env python3
"""Independently rebuild coverage/geometry joins using dataframe merges."""
import argparse,json,hashlib
from pathlib import Path
import pandas as pd


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def table(p):return pd.read_csv(p,sep='\t',dtype=str,keep_default_na=False)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    plan=json.loads(a.plan.read_text());root=Path(plan['output']);r=json.loads((root/'receipt.json').read_text())
    assert r['status']=='complete_domain_coverage_geometry_join_pending_readback' and r['plan_sha256']==sha(a.plan)
    for p,h in plan['pins'].items():assert sha(p)==h
    for b in r['source_bindings'].values():
        for p,h in b.items():assert sha(p)==h
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    checked={}
    for label,spec in plan['cohorts'].items():
        croot=Path(spec['coverage']);groot=Path(spec['geometry'])
        for folder in [croot,groot]:
            receipt=json.loads((folder/'receipt.json').read_text())
            for name,h in receipt['artifacts'].items():assert sha(folder/name)==h
        c=table(croot/'pair_mask_coverage.tsv');g=table(groot/'alignment_geometry.tsv');actual=table(root/(label+'.tsv'))
        keys=['pair_key','mask'];assert not c.duplicated(keys).any() and not g.duplicated(keys+['order']).any()
        expected=c[keys+['interval_a','interval_b']].copy();unique=pd.Series(True,index=c.index);degenerate=pd.Series(False,index=c.index)
        for order in ['0','1']:
            subset=g[g.order==order][keys+['geometry_status','aligned_length','rmsd_status']]
            joined=c.merge(subset,on=keys,how='left',validate='one_to_one',indicator=True)
            assert len(joined)==len(c) and len(subset)==sum(joined['_merge']=='both')
            missing=joined['_merge']=='left_only';assert (missing==(c['order'+order+'_status']=='input_unavailable')).all()
            for col,source in [('aligned_length','aligned_length'),('rmsd_status','status')]:
                assert (joined.loc[~missing,col].values==c.loc[~missing,'order'+order+'_'+source].values).all()
            status=joined.geometry_status.fillna('input_unavailable');expected['order'+order+'_geometry_status']=status
            unique &= status.eq('unique_at_numeric_tolerance');degenerate |= status.eq('degenerate_at_numeric_tolerance')
        expected['both_rotations_unique']=unique.astype(int).astype(str);counts={};joint={}
        for sid in plan['screens']:
            passed=c[sid+'_pass'].eq('1') & unique
            expected[sid+'_pass']=passed.astype(int).astype(str)
            why=c[sid+'_exclusions'].copy()
            expected[sid+'_exclusions']=[(';'.join(filter(None,[x,'nonunique_rotation']))) if d else x for x,d in zip(why,degenerate)]
            for mask in ['full','plddt70']:counts[mask+':'+sid]=int((passed & c['mask'].eq(mask)).sum())
            joint[sid]=int(pd.DataFrame({'pair':c.pair_key,'pass':passed}).groupby('pair')['pass'].all().sum())
        assert set(actual.columns)==set(expected.columns)
        pd.testing.assert_frame_equal(actual[expected.columns].sort_values(keys).reset_index(drop=True),expected.sort_values(keys).reset_index(drop=True))
        summary=r['cohorts'][label];assert summary['pair_mask_rows']==len(c) and summary['geometry_rows']==len(g)
        assert summary['pass_counts']==counts and summary['both_masks_pass_counts']==joint
        checked[label]=dict(rows=len(c),geometry_rows=len(g),pass_counts=counts,both_masks_pass_counts=joint)
    result=dict(status='passed_full_domain_coverage_geometry_readback',source_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(__file__),cohorts=checked,scope='Every joined field, missing disposition, numerical-status binding, six screen flags, exclusion list and aggregate independently rebuilt by dataframe joins. No new coordinate fits or biological inference.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
