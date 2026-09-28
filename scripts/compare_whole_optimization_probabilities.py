#!/usr/bin/env python3
"""Measure probability changes across optimization stages, retaining every site."""
import csv,gzip,json,subprocess,time
from pathlib import Path
from collections import defaultdict
import numpy as np
import psutil
from prepare_case_ancestral_neighborhoods import sha

AA='ARNDCQEGHILKMFPSTWYV'

def main():
    pp=Path('metadata/whole_optimization_probability_comparison_plan_20260927.json');plan=json.loads(pp.read_text())
    def verify():
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify();collections={};cache={}
    mappings={}
    for label,config in plan['collections'].items():
        identity=config.get('auditor')
        if identity:
            while True:
                s=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',identity['unit'],'--property=ActiveState,MainPID,Result,ExecMainStatus'],text=True).splitlines())
                if s['ActiveState']=='inactive':assert s['Result']=='success' and s['ExecMainStatus']=='0';break
                assert s['ActiveState'] in ['active','activating','deactivating'];p=psutil.Process(identity['pid']);assert p.create_time()==identity['created'] and p.cmdline()==identity['cmdline'] and int(s['MainPID'])==p.pid
                print('waiting_for_posterior_audit',label,flush=True);time.sleep(30)
        root=Path(config['root']);ar=Path(config['audit_receipt']);audit=json.loads(ar.read_text());assert audit['status']==config['expected_status'] and audit['source_receipt_sha256']==sha(root/'receipt.json')
        for f,h in audit['artifacts'].items():assert sha(ar.parent/f)==h
        fit_audit_path=Path(config['fit_audit']); fit_audit=json.loads(fit_audit_path.read_text())
        r=json.loads((root/'receipt.json').read_text())
        assert r['audit_receipt_sha256']==sha(fit_audit_path)
        assert fit_audit['source_receipt_sha256']==sha(config['fits_receipt'])
        assert fit_audit['artifacts']['fit_readback.tsv']==sha(config['fit_readback'])
        jobs={x['job']['job_id']:x['job'] for x in json.loads(Path(config['fits_receipt']).read_text())['results']}
        ll={x['job_id']:float(x['checkpoint_log_likelihood']) for x in csv.DictReader(Path(config['fit_readback']).open(),delimiter='\t')}
        assert set(jobs)==set(ll)
        mapping_path=root/'candidate_mapping.tsv';assert sha(mapping_path)==r['artifacts']['candidate_mapping.tsv']
        mappings[label]={(v['job_id'],int(v['level'])):v['partition_signature_sha256'] for v in csv.DictReader(mapping_path.open(),delimiter='\t') if v['status']=='matched_unrooted_vertex'}
        assert len(jobs)==config['expected_fits']
        for job in jobs:
            path=root/(job+'.npz');assert sha(path)==r['artifacts'][path.name]
            with np.load(path,allow_pickle=False) as z:
                assert np.array_equal(z['levels'],[0,1,2]) and ''.join(z['amino_acids'])==AA;cache[label,job]=z['posterior'].copy()
        collections[label]=(jobs,ll)
    comparisons=[]
    for name,job in collections['refined'][0].items():comparisons.append(('refined_versus_baseline','refined',name,'baseline',job['base_job_id']))
    pairs=defaultdict(dict)
    for name,job in collections['refined'][0].items():pairs[job['base_job_id']][job['alpha_min']]=name
    for base,pair in pairs.items():assert set(pair)=={.02,.005};comparisons.append(('refined_lower_versus_original_bound','refined',pair[.005],'refined',pair[.02]))
    assert len(comparisons)==234
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);summary=[];nrows=0
    fields=['comparison','base_job_id','left_collection','left_job','right_collection','right_job','level','column','left_map','right_map','left_maximum_probability','right_maximum_probability','total_variation','map_disagrees','both_at_least_090_disagree']
    with gzip.open(out/'site_comparisons.tsv.gz','wt') as h:
        writer=csv.DictWriter(h,fields,delimiter='\t',lineterminator='\n');writer.writeheader()
        for kind,lc,lj,rc,rj in comparisons:
            a,b=cache[lc,lj],cache[rc,rj];left,right=collections[lc][0][lj],collections[rc][0][rj]
            assert left['alignment_sha256']==right['alignment_sha256'] and a.shape==b.shape
            base=left.get('base_job_id',lj);assert base==right.get('base_job_id',rj)
            assert all(mappings[lc][lj,level]==mappings[rc][rj,level] for level in range(3))
            tv=abs(a-b).sum(axis=2)/2;am,bm=a.argmax(axis=2),b.argmax(axis=2);ap,bp=a.max(axis=2),b.max(axis=2);diff=am!=bm;high=diff&(ap>=.9)&(bp>=.9)
            for level in range(3):
                for col in range(a.shape[1]):
                    writer.writerow(dict(comparison=kind,base_job_id=base,left_collection=lc,left_job=lj,right_collection=rc,right_job=rj,level=level,column=col+1,left_map=AA[am[level,col]],right_map=AA[bm[level,col]],left_maximum_probability=float(ap[level,col]),right_maximum_probability=float(bp[level,col]),total_variation=float(tv[level,col]),map_disagrees=int(diff[level,col]),both_at_least_090_disagree=int(high[level,col])));nrows+=1
                summary.append(dict(comparison=kind,base_job_id=base,left_collection=lc,left_job=lj,right_collection=rc,right_job=rj,level=level,columns=a.shape[1],left_minus_right_log_likelihood=collections[lc][1][lj]-collections[rc][1][rj],map_disagreements=int(diff[level].sum()),both_at_least_090_disagreements=int(high[level].sum()),mean_total_variation=float(tv[level].mean()),maximum_total_variation=float(tv[level].max())))
    assert len(summary)==702 and nrows==653157
    with (out/'comparison_summary.tsv').open('w') as h:w=csv.DictWriter(h,list(summary[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(summary)
    agg=defaultdict(lambda:dict(sites=0,map_disagreements=0,both_at_least_090_disagreements=0,tv_sum=0.,maximum_total_variation=0.))
    with gzip.open(out/'site_comparisons.tsv.gz','rt') as h:
        for r in csv.DictReader(h,delimiter='\t'):
            v=agg[r['comparison']];v['sites']+=1;v['map_disagreements']+=int(r['left_map']!=r['right_map']);v['both_at_least_090_disagreements']+=int(r['left_map']!=r['right_map'] and min(float(r['left_maximum_probability']),float(r['right_maximum_probability']))>=.9);v['tv_sum']+=float(r['total_variation']);v['maximum_total_variation']=max(v['maximum_total_variation'],float(r['total_variation']))
    for kind,v in agg.items():
        rs=[r for r in summary if r['comparison']==kind];assert v['sites']==sum(r['columns'] for r in rs) and v['map_disagreements']==sum(r['map_disagreements'] for r in rs) and v['both_at_least_090_disagreements']==sum(r['both_at_least_090_disagreements'] for r in rs)
        assert abs(v['tv_sum']-sum(r['columns']*r['mean_total_variation'] for r in rs))<1e-7
        v['mean_total_variation']=v.pop('tv_sum')/v['sites']
    verify();r=dict(status='complete_whole_optimization_probability_comparisons',fit_pairs=234,node_comparisons=702,site_comparisons=nrows,descriptive_aggregates=dict(agg),plan_sha256=sha(pp),input_posterior_audit_hashes={k:sha(v['audit_receipt']) for k,v in plan['collections'].items()},artifacts={p.name:sha(p) for p in out.iterdir()},scope='Same-alignment coordinate comparisons across audited optimization stages. All sites from baseline/refined and bound comparisons retained; alternate-start probabilities are not included; these are dependent diagnostics, not biological replicates or ensemble weights. Historical state support, indel uncertainty and model adequacy remain separate.')
    (out/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)

if __name__=='__main__':main()
