#!/usr/bin/env python3
"""Reconstruct raw loadings, latent-space Grams and independent LU benchmarks."""
import argparse
from collections import Counter,defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy import sparse,linalg
from full_entity_operator_sources import load,digest,rhs,benchmark_variances,KINDS,MODES,KERNELS,SUMMARY_FIELDS
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def lu(core):
    factors,pivots=linalg.lu_factor(core,check_finite=True)
    diagonal=np.diag(factors);sign=int(np.prod(np.sign(diagonal)))*(-1)**int(np.count_nonzero(pivots!=np.arange(len(pivots))))
    assert sign==1 and np.all(diagonal!=0)
    return (factors,pivots),float(np.log(abs(diagonal)).sum())


def independent_benchmark(labels,operators,factor,values):
    levels=sorted(set(labels));parts=[np.flatnonzero(labels==v) for v in levels];lowers=[];logdet=0.;variances=benchmark_variances()
    for rows in parts:
        base=np.eye(len(rows))
        for name,z in operators.items():
            local=z[rows];local=local[:,np.unique(local.indices)].tocsc()
            for column in range(local.shape[1]):
                start,stop=local.indptr[column:column+2];indices=local.indices[start:stop];loadings=local.data[start:stop]
                base[np.ix_(indices,indices)]+=variances[name]*np.outer(loadings,loadings)
        lower,ld=lu(base);lowers.append(lower);logdet+=ld
    def base_solve(v):
        out=np.empty_like(v)
        for rows,lower in zip(parts,lowers):out[rows]=linalg.lu_solve(lower,v[rows],check_finite=False)
        return out
    bi=base_solve(values);bf=base_solve(factor)
    core=np.eye(factor.shape[1])+.25*factor.T@bf;lower,ld=lu(core);logdet+=ld
    solved=bi-.25*bf@linalg.lu_solve(lower,factor.T@bi,check_finite=False)
    return solved,logdet


def run(path,output):
    path=Path(path);plan=json.loads(path.read_text());source,bindings=load(plan,path);original=dict(bindings)
    root=Path(plan['output']);rp=root/'receipt.json';receipt=json.loads(rp.read_text())
    assert receipt['status']=='complete_full_entity_operator_bank_pending_independent_readback' and receipt['plan_sha256']==sha(path)
    assert receipt['source_hashes']==original and receipt['scientific_eligibility'] is False
    bind(bindings,rp)
    for name,d in receipt['artifacts'].items():bind(bindings,root/name,d)
    verify(bindings);assert json.loads((root/'stage_plan.json').read_text())==dict(plan_sha256=sha(path),schema='expanded-entity-operators-v1')
    cases=source['cases'];n=len(cases);index={r['case_id']:i for i,r in enumerate(cases)}
    assert json.loads((root/'case_ids.json').read_text())==[r['case_id'] for r in cases]
    labels=np.load(root/'block_labels.npy');np.testing.assert_array_equal(labels,np.asarray([r['family_component'] for r in cases],dtype='S64'))
    labelmap=json.loads((root/'entity_labels.json').read_text());assert set(labelmap)==set(KINDS)
    for k in KINDS:assert labelmap[k]==sorted(set(labelmap[k])) and len(labelmap[k])==source['completion']['unique_entities'][k]
    positions={k:{v:i for i,v in enumerate(labelmap[k])} for k in KINDS};cells=defaultdict(Counter);seen=Counter();occurrences=Counter()
    with gzip.open(source['incidence'],'rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            k=r['entity_kind'];token=json.loads(r['entity_value']);side=r['side'];assert k in KINDS and side in ['target','background']
            assert r['entity_id']==hashlib.sha256(json.dumps(['expanded-covariance-'+k+'-v1',token],separators=(',',':')).encode()).hexdigest()
            magnitude=.5 if k in ['gene','model'] else 1.;signed=magnitude if side=='target' else -magnitude
            assert float(r['signed_loading'])==signed and float(r['unsigned_loading'])==magnitude
            assert r['endpoint'] in ['a','b'] if k in ['gene','model'] else r['endpoint']==''
            if k.endswith('_node'):assert k==side+'_node'
            row=index[r['case_id']];column=positions[k][r['entity_id']]
            cells['signed',k][row,column]+=signed;cells['unsigned',k][row,column]+=magnitude
            seen[r['case_id']]+=1;occurrences[k]+=1
    assert set(seen)==set(index) and set(seen.values())=={14} and sum(occurrences.values())==source['completion']['entity_occurrences']
    entries=json.loads((root/'operator_manifest.json').read_text());expected_keys={(m,k) for m in MODES for k in KINDS}|{('contrast','family_intercept')}
    assert len(entries)==len(expected_keys)==13 and {(r['mode'],r['kind']) for r in entries}==expected_keys
    operators={}
    for r in entries:
        m,k=r['mode'],r['kind'];expected_path='operators/'+(m+'-'+k if m!='contrast' else 'family_intercept')+'.npz'
        assert r['path']==expected_path and r['sha256']==receipt['artifacts'][r['path']]
        z=sparse.load_npz(root/r['path']).tocsr();assert z.has_canonical_format and np.isfinite(z.data).all() and np.all(z.data!=0)
        assert list(z.shape)==r['shape']==[n,len(labelmap[k if k!='family_intercept' else 'family'])] and z.nnz==r['nnz']
        saved={(i,int(j)):float(v) for i in range(n) for j,v in zip(z.indices[z.indptr[i]:z.indptr[i+1]],z.data[z.indptr[i]:z.indptr[i+1]])}
        if m=='contrast':
            expected={(i,positions['family'][digest(['expanded-covariance-family-v1',[c['target_family']]])]):1. for i,c in enumerate(cases)}
        else:expected={key:v for key,v in cells[m,k].items() if v!=0}
        assert saved==expected;operators[m,k]=z
        # Block factorization must never ignore a shared entity across components.
        bycolumn={}
        for (i,j) in saved:
            if j in bycolumn:assert bycolumn[j]==labels[i]
            else:bycolumn[j]=labels[i]
    family=operators['contrast','family_intercept'];grams=json.loads((root/'kernel_gram_manifest.json').read_text());benchmarks=json.loads((root/'benchmark_manifest.json').read_text())
    assert [(r['mode'],r['tree']) for r in grams]==[(m,t) for t in plan['trees'] for m in MODES]
    assert [(r['mode'],r['tree']) for r in benchmarks]==[(m,t) for t in plan['trees'] for m in MODES]
    fixed={}
    for m in MODES:
        matrices=[operators[m,k] for k in KINDS]+[family];g=np.zeros((len(KERNELS),len(KERNELS)));g[0,0]=n
        for i,z in enumerate(matrices,1):
            g[0,i]=g[i,0]=float(z.multiply(z).sum())
            for j,w in enumerate(matrices,1):
                cross=z.T@w;g[i,j]=float(cross.multiply(cross).sum())
        fixed[m]=g
    patterns=np.asarray([int(c['species_pattern_row']) for c in cases]);maximum_gram_error=maximum_benchmark_error=0.
    for tree in plan['trees']:
        with np.load(source['root']/(tree+'.npz')) as a:original_factor=a['factor'];factor=original_factor[patterns]
        multiplicity=np.bincount(patterns,minlength=len(original_factor));core=original_factor.T@(multiplicity[:,None]*original_factor)
        for m in MODES:
            e=next(r for r in grams if r['mode']==m and r['tree']==tree);g=fixed[m].copy()
            g[-1,-1]=float(np.sum(core*core));g[0,-1]=g[-1,0]=float(np.sum(multiplicity*np.sum(original_factor**2,axis=1)))
            for k,z in enumerate([operators[m,k] for k in KINDS]+[family],1):
                projected=z.T@factor;g[k,-1]=g[-1,k]=float(np.sum(projected*projected));del projected
            assert e['path']==m+'-'+tree+'-kernel-gram.npz' and e['sha256']==receipt['artifacts'][e['path']]
            with np.load(root/e['path']) as a:assert a.files==['gram'];saved=a['gram']
            np.testing.assert_allclose(saved,g,rtol=2e-10,atol=1e-7);maximum_gram_error=max(maximum_gram_error,float(np.max(abs(saved-g))))
            diagonal=np.diag(g);active=np.flatnonzero(diagonal>0);normalized=g[np.ix_(active,active)]/np.sqrt(diagonal[active,None]*diagonal[None,active])
            singular=linalg.svd(normalized,compute_uv=False,lapack_driver='gesvd');tol=len(active)*np.finfo(float).eps*singular[0]
            assert e['kernel_names']==KERNELS and e['zero_kernels']==[KERNELS[i] for i in np.flatnonzero(diagonal==0)] and e['active_kernel_indices']==active.tolist()
            np.testing.assert_allclose(e['normalized_gram_singular_values'],singular,rtol=1e-8,atol=1e-12)
            np.testing.assert_allclose(e['rank_tolerance'],tol,rtol=1e-8,atol=0)
            assert e['normalized_kernel_gram_rank']==int((singular>tol).sum())
            b=next(r for r in benchmarks if r['mode']==m and r['tree']==tree)
            assert b['path']=='benchmarks/'+m+'-'+tree+'.npz' and b['sha256']==receipt['artifacts'][b['path']]
            assert b['entity_variances']==benchmark_variances() and b['species_variance']==.25 and b['residual_diagonal']==1.
            reference,logdet=independent_benchmark(labels,{**{k:operators[m,k] for k in KINDS},'family_intercept':family},factor,rhs(n))
            with np.load(root/b['path']) as a:assert a.files==['solution'];solution=a['solution']
            np.testing.assert_allclose(solution,reference,rtol=2e-8,atol=2e-10)
            np.testing.assert_allclose(b['logdet'],logdet,rtol=2e-9,atol=2e-7)
            assert b['maximum_solve_residual']<=1e-9 and b['elapsed_seconds']>=0
            maximum_benchmark_error=max(maximum_benchmark_error,float(np.max(abs(solution-reference))))
            print('independent_full_entity_operator_benchmark',m,tree,n,flush=True)
        del factor
    counts=Counter(labels.tolist())
    summary=dict(logical_cases=n,entity_occurrences=sum(occurrences.values()),unique_entities={k:len(v) for k,v in labelmap.items()},matrices=len(entries),
        operator_nnz={r['mode']+'|'+r['kind']:r['nnz'] for r in entries},family_components=len(counts),largest_component=max(counts.values()),sum_component_squared_sizes=sum(v*v for v in counts.values()),
        trees=plan['trees'],kernel_gram_ranks={r['mode']+'|'+r['tree']:r['normalized_kernel_gram_rank'] for r in grams},zero_kernels={r['mode']+'|'+r['tree']:r['zero_kernels'] for r in grams},numerical_benchmarks=len(benchmarks))
    assert all(receipt[k]==summary[k] for k in SUMMARY_FIELDS)
    expected_artifacts={'stage_plan.json','operator_manifest.json','entity_labels.json','case_ids.json','block_labels.npy','kernel_gram_manifest.json','benchmark_manifest.json'}|{r['path'] for r in entries+grams+benchmarks}
    assert set(receipt['artifacts'])==expected_artifacts;verify(bindings)
    result=dict(status='passed_full_entity_operator_raw_loading_and_gram_readback',plan_sha256=sha(path),producer_receipt_sha256=sha(rp),**summary,
        maximum_absolute_kernel_gram_error=maximum_gram_error,maximum_absolute_benchmark_solve_error=maximum_benchmark_error,source_hashes=bindings,scientific_eligibility=False,scope=plan['scope'])
    with Path(output).open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(summary),flush=True);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.plan,a.output)
