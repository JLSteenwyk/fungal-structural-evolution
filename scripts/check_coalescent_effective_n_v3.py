#!/usr/bin/env python3
"""Check native effective-N behavior around 0.001 and preserve strict tolerances."""
import argparse
import json
from pathlib import Path
import subprocess

import numpy as np

from coalescent_quartet_audit_v2 import annotated_branches,colored_quartets,packed,read_tree
from coalescent_quartet_numeric_v3 import numerical_check,effective_n
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def balanced(labels):
    if len(labels)==1:return labels[0]
    mid=len(labels)//2
    return '('+balanced(labels[:mid])+','+balanced(labels[mid:])+')'


def check_branch(branch,genes,tips):
    counts=np.zeros(3);available_genes=0
    for gene in genes:
        c,available=colored_quartets(*packed(gene,tips),branch['colors'])
        if available:counts+=c/available;available_genes+=1
    summary=numerical_check(branch['native'],branch['length'],counts,available_genes)
    return counts,available_genes,summary


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=False)
    plan_path='metadata/species_coalescent_native_plan_20261002.json';plan=json.loads(Path(plan_path).read_text())
    jar=plan['jar'];bindings={}
    for p in [str(Path(__file__)),'scripts/coalescent_quartet_numeric_v3.py','scripts/coalescent_quartet_audit_v2.py',
              plan_path,jar,plan['java'],'results/software-checks/astral-effective-n-rule-source-20261002-v1/provenance.json']:
        bind(bindings,p)
    source=Path('results/software-checks/astral-effective-n-rule-source-20261002-v1')
    provenance=json.loads((source/'provenance.json').read_text())
    assert provenance['installed_jar_sha256']==sha(jar) and provenance['installed_class_contains_0_001_double']
    for name,row in provenance['source_files'].items():bind(bindings,source/name,row['sha256'])
    existing=Path('results/software-checks/coalescent-quartet-numeric-contracts-20261002-v2/receipt.json')
    prior=json.loads(existing.read_text());bind(bindings,existing)
    assert prior['status']=='passed_exhaustive_quartet_dp_and_native_coalescent_numeric_contracts'
    for p,d in prior['source_hashes'].items():bind(bindings,p,d)
    rejected=branches=0
    for row in prior['native_results']:
        for p,d in row['artifacts'].items():bind(bindings,p,d)
        gp=Path(row['command'][row['command'].index('-i')+1]);sp=Path(row['command'][row['command'].index('-o')+1])
        genes=[read_tree(data=n) for n in gp.read_text().splitlines()];tips=sorted('ABCDEF')
        for branch in annotated_branches(read_tree(path=sp),tips):
            counts,available,summary=check_branch(branch,genes,tips);branches+=1
            for key in ['f1','f2','f3','pp1','pp2','pp3','EN']:
                altered=dict(branch['native']);altered[key]+=.0001
                try:numerical_check(altered,branch['length'],counts,available)
                except AssertionError:rejected+=1
                else:raise AssertionError('Accepted altered value '+key)
            try:numerical_check(branch['native'],branch['length']+.0001,counts,available)
            except AssertionError:rejected+=1
            else:raise AssertionError('Accepted altered length')
    boundary_rows=[]
    for denominator in [999,1000,1001]:
        folder=out/str(denominator);folder.mkdir()
        alabels=['A'+str(i).zfill(4) for i in range(denominator)]
        reference='(('+balanced(alabels)+',B),C,D);'
        gene='(('+','.join(alabels[1:]+['B'])+'),'+alabels[0]+',C,D);'
        gp,rp,sp,lp=[folder/name for name in ['genes.tree','reference.tree','species.tree','native.log']]
        gp.write_text(gene+'\n');rp.write_text(reference+'\n')
        cmd=[plan['java'],'-XX:ActiveProcessorCount=1','-Xmx1G','-jar',jar,'-i',str(gp.resolve()),
             '-q',str(rp.resolve()),'-o',str(sp.resolve()),'-t','2','-s','20261002']
        with lp.open('x') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=120)
        tips=sorted(alabels+['B','C','D'])
        branch=next(b for b in annotated_branches(read_tree(path=sp),tips) if b['canonical_side']==['C','D'])
        exact,available_quartets=colored_quartets(*packed(read_tree(data=gene),tips),branch['colors'])
        assert available_quartets==denominator and sorted(exact)==[0,0,denominator-1]
        counts,available,summary=check_branch(branch,[read_tree(data=gene)],tips)
        assert available==1
        assert summary['native_effective_n_adjustment_applied']==(abs(float(counts.sum())-1)>.001)
        if denominator==1001:assert branch['native']['EN']==1. and float(counts.sum())<1
        else:assert abs(branch['native']['EN']-(denominator-1)/denominator)<1e-12
        for p in [gp,rp,sp,lp]:bind(bindings,p)
        boundary_rows.append(dict(denominator=denominator,unresolved_quartets=1,available_genes=available,
            exact_counts=exact.tolist(),native=branch['native'],native_map_length=branch['length'],**summary))
        print('native_effective_n_boundary_passed',denominator,branch['native']['EN'],float(counts.sum()),flush=True)
    # A claimed available-gene count cannot replace independent fractional evidence arbitrarily.
    assert effective_n([.8,.1,.1],1)==1.
    verify(bindings)
    result=dict(status='passed_native_effective_n_threshold_and_strict_numeric_contracts',
        previous_native_cases=5,previous_branches=branches,boundary_native_cases=3,boundary_rows=boundary_rows,
        altered_native_values_rejected=rejected,comparison_tolerance=2e-8,native_algorithm_threshold=.001,
        source_hashes=bindings,scientific_eligibility=False,
        scope='Pinned installed ASTRAL5.7.8 jar: three fixed-reference contracts with one unresolved quartet among999/1000/1001 possibilities; below/near/above native .001 substitution boundary. Recheck all five previous native resolved/fractional/missing/star/discordant cases. All0.0001 count/EN/PP/length mutations rejected; comparison tolerance unchanged2e-8. Only effective-N algorithm semantics corrected. Full production readback and scientific qualification remain required.')
    with (out/'receipt.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if isinstance(v,(str,int,float,bool))}),flush=True)


if __name__=='__main__':main()
