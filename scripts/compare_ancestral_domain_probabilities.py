#!/usr/bin/env python3
"""Compare all fitted ancestral contexts only at exact full-clade residue maps."""
import csv,gzip,json,itertools
from collections import defaultdict
from pathlib import Path
import numpy as np
from Bio import SeqIO
from prepare_case_ancestral_neighborhoods import sha,read

AA='ARNDCQEGHILKMFPSTWYV'

def main():
    pp=Path('metadata/ancestral_domain_probability_comparison_plan_20260927.json');plan=json.loads(pp.read_text())
    def verify():
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify();root=Path(plan['probabilities']);receipt=json.loads((root/'receipt.json').read_text())
    for p,h in receipt['artifacts'].items():assert sha(root/p)==h
    jobs=[r['job'] for r in json.loads(Path(plan['fits_receipt']).read_text())['results']];assert len(jobs)==156
    positions=defaultdict(lambda:defaultdict(list));count=0
    with gzip.open(plan['positions'],'rt') as h:
        for r in csv.DictReader(h,delimiter='\t'):
            key=(r['family'],r['boundary'],r['method']);positions[key][int(r['column'])].append((r['gene'],int(r['protein_position'])));count+=1
    assert count==637340
    maps={};genes={};columns=[]
    for key,cols in positions.items():
        maps[key]={tuple(sorted(values)):col for col,values in cols.items()};assert len(maps[key])==len(cols)
        genes[key]={g for values in cols.values() for g,p in values}
        for sig,col in maps[key].items():columns.append(dict(family=key[0],boundary=key[1],method=key[2],column=col,signature=json.dumps(sig,separators=(',',':'))))
    assert len(maps)==52 and len(columns)==9349
    arrays={};byfamily=defaultdict(list)
    for job in jobs:
        key=(job['family'],job['boundary'],job['method']);assert len(maps[key])==job['columns']
        assert sha(job['alignment'])==job['alignment_sha256'];assert {r.id for r in SeqIO.parse(job['alignment'],'fasta')}==genes[key]
        with np.load(root/(job['job_id']+'.npz'),allow_pickle=False) as z:
            assert np.array_equal(z['levels'],[0,1,2]) and ''.join(z['amino_acids'])==AA
            arrays[job['job_id']]=z['posterior'].copy()
        byfamily[job['family']].append(job)
    nodes=read(plan['nodes']);partitions={(r['job_id'],int(r['level'])):r['partition_signature_sha256'] for r in nodes if r['guide']=='profile' and r['status']=='matched_unrooted_vertex'}
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);summaries=[];nrows=0;nshared=0
    fields=['family','job_a','job_b','level','column_a','column_b','status','map_a','map_b','map_disagrees','total_variation','maximum_probability_a','maximum_probability_b']
    with gzip.open(out/'site_comparisons.tsv.gz','wt') as h:
        w=csv.DictWriter(h,fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for family,fjobs in sorted(byfamily.items()):
            assert len(fjobs)==12
            for a,b in itertools.combinations(sorted(fjobs,key=lambda x:x['job_id']),2):
                ka=(family,a['boundary'],a['method']);kb=(family,b['boundary'],b['method']);assert genes[ka]==genes[kb]
                ma,mb=maps[ka],maps[kb];union=sorted(set(ma)|set(mb));common=set(ma)&set(mb)
                for level in range(3):
                    assert partitions[a['job_id'],level]==partitions[b['job_id'],level]
                    tvs=[];disagreements=0;both_high=0
                    for sig in union:
                        ca,cb=ma.get(sig),mb.get(sig);status='matched_full_clade_column' if ca and cb else ('only_a' if ca else 'only_b')
                        row=dict(family=family,job_a=a['job_id'],job_b=b['job_id'],level=level,column_a=ca or '',column_b=cb or '',status=status)
                        if ca and cb:
                            pa,pb=arrays[a['job_id']][level,ca-1],arrays[b['job_id']][level,cb-1]
                            tv=float(abs(pa-pb).sum()/2);assert 0<=tv<=1+1e-12
                            da=int(pa.argmax()!=pb.argmax());disagreements+=da;tvs.append(tv);both_high+=int(da and pa.max()>=.9 and pb.max()>=.9)
                            row.update(map_a=AA[pa.argmax()],map_b=AA[pb.argmax()],map_disagrees=da,total_variation=tv,maximum_probability_a=float(pa.max()),maximum_probability_b=float(pb.max()));nshared+=1
                        w.writerow(row);nrows+=1
                    summaries.append(dict(family=family,job_a=a['job_id'],job_b=b['job_id'],level=level,boundary_changed=int(a['boundary']!=b['boundary']),aligner_changed=int(a['method']!=b['method']),model_changed=int(a['model']!=b['model']),columns_a=len(ma),columns_b=len(mb),matched_columns=len(common),unmatched_a=len(ma)-len(common),unmatched_b=len(mb)-len(common),map_disagreements=disagreements,both_at_least_090_disagreements=both_high,mean_total_variation=float(np.mean(tvs)) if tvs else '',maximum_total_variation=max(tvs) if tvs else ''))
            print(family,'all_66_context_pairs_complete',flush=True)
    assert len(summaries)==2574
    for name,rows in [('comparison_summary.tsv',summaries),('column_signatures.tsv',columns)]:
        with (out/name).open('w') as h:
            w=csv.DictWriter(h,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    # Independently aggregate the saved site table, including unmatched records.
    agg=defaultdict(lambda:dict(shared=0,only_a=0,only_b=0,disagree=0,high=0,tv_sum=0.,tv_max=0.));seen=set()
    with gzip.open(out/'site_comparisons.tsv.gz','rt') as h:
        for r in csv.DictReader(h,delimiter='\t'):
            key=(r['job_a'],r['job_b'],int(r['level']));identity=key+(r['column_a'],r['column_b']);assert identity not in seen;seen.add(identity);v=agg[key]
            if r['status']=='matched_full_clade_column':
                v['shared']+=1;v['disagree']+=int(r['map_a']!=r['map_b']);v['high']+=int(r['map_a']!=r['map_b'] and min(float(r['maximum_probability_a']),float(r['maximum_probability_b']))>=.9);v['tv_sum']+=float(r['total_variation']);v['tv_max']=max(v['tv_max'],float(r['total_variation']))
            else:v[r['status']]+=1;assert r['total_variation']=='' and r['map_disagrees']==''
    assert len(seen)==nrows and len(agg)==2574
    for r in summaries:
        v=agg[r['job_a'],r['job_b'],r['level']]
        assert [v[k] for k in ['shared','only_a','only_b','disagree','high']]==[r[k] for k in ['matched_columns','unmatched_a','unmatched_b','map_disagreements','both_at_least_090_disagreements']]
        if v['shared']:assert abs(v['tv_sum']/v['shared']-r['mean_total_variation'])<1e-12 and v['tv_max']==r['maximum_total_variation']
    verify();r=dict(status='complete_full_context_ancestral_probability_comparison',families=13,context_pairs=858,node_comparisons=len(summaries),site_union_rows=nrows,matched_site_comparisons=nshared,plan_sha256=sha(pp),artifacts={p.name:sha(p) for p in out.iterdir()},scope='Descriptive sensitivity conditional on fitted models. Exact full-clade original-protein residue signatures; unmatched columns explicitly retained. MAP ties use fixed AA order, so near ties can switch. Context pairs and sites are dependent, not replicates. No model averaging, convergence, indel inference or validated historical sequence claim.')
    (out/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)

if __name__=='__main__':main()
