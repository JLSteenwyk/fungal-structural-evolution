#!/usr/bin/env python3
"""Apply the declared compatible-gap surrogate to the complete candidate set."""
import csv,gzip,json
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
from compatible_gap_distribution import CompatibleGaps
from prepare_case_ancestral_neighborhoods import sha


def main():
    pp=Path('metadata/compatible_gap_application_plan_20260927.json');plan=json.loads(pp.read_text())
    for p,h in plan['pins'].items():assert sha(p)==h,p
    root=Path(plan['posterior_output']);pr=json.loads((root/'receipt.json').read_text());ar=json.loads(Path(plan['posterior_audit']).read_text());assert ar['source_receipt_sha256']==sha(root/'receipt.json')
    assert ar['nonempty_models']==306 and ar['status']=='all_selected_fit_gap_probabilities_independently_verified'
    coords=defaultdict(list)
    with Path(plan['coordinates']).open() as h:
        for r in csv.DictReader(h,delimiter='\t'):coords[r['input_id']+'-'+r['terminal_policy']].append(r)
    mapping=defaultdict(list)
    with Path(plan['mapping']).open() as h:
        for r in csv.DictReader(h,delimiter='\t'):
            if r['guide']=='profile' and r['status']=='matched_unrooted_vertex':mapping[r['job_id']].append(r)
    jobs=json.loads(Path(plan['fit_plan']).read_text())['jobs'];out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    rows=[];character_count=0
    with gzip.open(out/'marginal_changes.tsv.gz','wt') as h:
        w=csv.writer(h,delimiter='\t',lineterminator='\n');w.writerow(['job_id','level','character','raw_probability','numerically_bounded_probability','conditional_probability','conditional_map_includes_gap'])
        for job in jobs:
            jid=job['job_id'];rp=root/jid/'receipt.json';assert sha(rp)==pr['job_receipts'][jid+'/receipt.json'];r=json.loads(rp.read_text())
            for name,digest in r['artifacts'].items():assert sha(root/jid/name)==digest
            encoding=jid.removesuffix('-all_taxa').removesuffix('-observed_mask');cs=sorted(coords[encoding],key=lambda x:int(x['character']));intervals=[(int(c['start_column']),int(c['end_column'])) for c in cs]
            assert len(intervals)==job['character_count'];targets=sorted(mapping[jid],key=lambda r:int(r['level']));assert len(targets)==3
            data=np.load(root/jid/'posterior.npz',allow_pickle=False) if intervals else None
            for target in targets:
                row=dict(job_id=jid,level=int(target['level']),characters=len(intervals),status='empty_no_gap_model',log_retained_mass=None,roundoff_bounds_adjusted=0,maximum_marginal_change=None,mean_marginal_change=None,map_state_switches=None,selected_exact_runs=None)
                if data is not None:
                    names=list(data['node_names']);raw=data['probability_gap'][names.index(target['posterior_node'])]
                    assert raw.min()>=-1e-12 and raw.max()<=1+1e-12
                    # Repair only already-audited floating-point excursions outside [0,1].
                    p=np.clip(raw,0,1);row['roundoff_bounds_adjusted']=int((p!=raw).sum())
                    model=CompatibleGaps(intervals,p)
                    if model.feasible:
                        q=model.marginals();selected=model.configuration();chosen=sorted([intervals[i] for i in selected]);assert all(a[1]+1<b[0] for a,b in zip(chosen,chosen[1:]))
                        assert np.all(q>=-1e-12) and np.all(q<=1+1e-12)
                        # Independent sweep checks every exclusion set, not just MAP validity.
                        for start in sorted({a for a,b in intervals}):
                            idx=[i for i,(a,b) in enumerate(intervals) if a<=start<=b+1]
                            assert q[idx].sum()<=1+1e-10
                        row.update(status='compatible_surrogate_computed',log_retained_mass=model.log_compatibility_probability,maximum_marginal_change=float(abs(q-raw).max()),mean_marginal_change=float(abs(q-raw).mean()),map_state_switches=int(((raw>=.5)!=(q>=.5)).sum()),selected_exact_runs=len(selected))
                        for i in range(len(p)):w.writerow([jid,target['level'],i+1,float(raw[i]),float(p[i]),float(q[i]),int(i in selected)])
                        character_count+=len(p)
                    else:row['status']='zero_compatible_mass_explicitly_unresolved'
                rows.append(row)
            if data is not None:data.close()
    assert len(rows)==936
    with (out/'node_summary.tsv').open('w') as h:w=csv.DictWriter(h,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    feasible=[r for r in rows if r['status']=='compatible_surrogate_computed']
    result=dict(status='complete_936_candidate_dispositions_under_compatible_surrogate',candidate_nodes=len(rows),status_counts=dict(Counter(r['status'] for r in rows)),character_marginals=character_count,roundoff_bounds_adjusted=sum(r['roundoff_bounds_adjusted'] for r in rows),minimum_log_retained_mass=min(r['log_retained_mass'] for r in feasible),maximum_marginal_change=max(r['maximum_marginal_change'] for r in feasible),map_state_switches=sum(r['map_state_switches'] for r in feasible),plan_sha256=sha(pp),artifacts={p.name:sha(p) for p in out.iterdir()},scope=plan['scope'])
    (out/'receipt.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print({k:v for k,v in result.items() if k not in ['artifacts','scope']})

if __name__=='__main__':main()
