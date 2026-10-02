#!/usr/bin/env python3
"""Export all shared-entity loadings and audit complete fixed numerical kernels."""
import argparse
from collections import Counter,defaultdict
import csv
import fcntl
import gzip
import json
from pathlib import Path
import shutil
import time
import numpy as np
from scipy import sparse
from full_entity_operator_sources import load,digest,gram_diagnostics,rhs,benchmark_variances,KINDS,MODES,KERNELS,SUMMARY_FIELDS
from shared_entity_covariance import SharedEntityCovariance
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def run(path,stop_after_exports=False):
    path=Path(path);plan=json.loads(path.read_text());source,bindings=load(plan,path);root=Path(plan['output'])
    assert shutil.disk_usage(root.parent).free>=plan['resources']['minimum_free_disk_gib']*2**30
    root.mkdir(exist_ok=True);lock=(root/'stage.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    assert not (root/'receipt.json').exists(),'Completed operator bank cannot restart'
    marker=root/'stage_plan.json';state=dict(plan_sha256=sha(path),schema='expanded-entity-operators-v1')
    if marker.exists():assert json.loads(marker.read_text())==state
    else:marker.write_text(json.dumps(state,indent=2)+'\n')
    cases=source['cases'];index={r['case_id']:i for i,r in enumerate(cases)};n=len(cases)
    rows={k:[] for k in KINDS};labels={k:[] for k in KINDS};values={(m,k):[] for m in MODES for k in KINDS};seen=Counter();counts=Counter()
    with gzip.open(source['incidence'],'rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            k=r['entity_kind'];assert k in KINDS;rows[k].append(index[r['case_id']]);labels[k].append(r['entity_id'])
            for m in MODES:values[m,k].append(float(r[m+'_loading']))
            seen[r['case_id']]+=1;counts[k]+=1
    assert sum(counts.values())==source['completion']['entity_occurrences'] and set(seen.values())=={14} and set(seen)==set(index)
    folder=root/'operators';folder.mkdir(exist_ok=True);operators={};entries=[];labelmap={}
    for k in KINDS:
        unique,inverse=np.unique(labels[k],return_inverse=True);labelmap[k]=unique.tolist()
        assert len(unique)==source['completion']['unique_entities'][k]
        for m in MODES:
            z=sparse.coo_matrix((values[m,k],(rows[k],inverse)),shape=(n,len(unique))).tocsr();z.sum_duplicates();z.eliminate_zeros();z.sort_indices()
            assert np.isfinite(z.data).all();operators[m,k]=z;fp=folder/(m+'-'+k+'.npz');sparse.save_npz(fp,z)
            entries.append(dict(mode=m,kind=k,path=str(fp.relative_to(root)),sha256=sha(fp),shape=list(z.shape),nnz=z.nnz))
    family_positions={v:i for i,v in enumerate(labelmap['family'])}
    famids=[digest(['expanded-covariance-family-v1',[c['target_family']]]) for c in cases]
    family=sparse.csr_matrix((np.ones(n),(np.arange(n),[family_positions[v] for v in famids])),shape=(n,len(family_positions)))
    assert operators['signed','family'].nnz==0
    assert (operators['unsigned','family']-2*family).nnz==0
    fp=folder/'family_intercept.npz';sparse.save_npz(fp,family)
    entries.append(dict(mode='contrast',kind='family_intercept',path=str(fp.relative_to(root)),sha256=sha(fp),shape=list(family.shape),nnz=family.nnz))
    (root/'entity_labels.json').write_text(json.dumps(labelmap,indent=2)+'\n')
    (root/'case_ids.json').write_text(json.dumps([c['case_id'] for c in cases],indent=2)+'\n')
    blocklabels=np.asarray([c['family_component'] for c in cases],dtype='S64');np.save(root/'block_labels.npy',blocklabels)
    (root/'operator_manifest.json').write_text(json.dumps(entries,indent=2)+'\n')
    if stop_after_exports:raise InterruptedError('Software interruption contract')
    codes=np.unique(blocklabels,return_inverse=True)[1];order=np.argsort(codes,kind='stable');parts=np.split(order,np.flatnonzero(np.diff(codes[order]))+1)
    fixed={};traces={}
    for m in MODES:
        gram=np.zeros((len(KERNELS),len(KERNELS)));gram[0,0]=n
        matrices=[operators[m,k] for k in KINDS]+[family]
        for part in parts:
            kernels=np.asarray([(z[part]@z[part].T).toarray() for z in matrices])
            flat=kernels.reshape(len(matrices),-1);gram[1:-1,1:-1]+=flat@flat.T
            diagonal=np.trace(kernels,axis1=1,axis2=2);gram[0,1:-1]+=diagonal;gram[1:-1,0]+=diagonal
        fixed[m]=gram;traces[m]=matrices
    allgrams=[];benchmarks=[];folder=root/'benchmarks';folder.mkdir(exist_ok=True);v=benchmark_variances()
    patterns=np.asarray([int(c['species_pattern_row']) for c in cases]);rhS=rhs(n)
    for tree in plan['trees']:
        with np.load(source['root']/(tree+'.npz')) as a:factor=a['factor'][patterns]
        core=factor.T@factor
        for m in MODES:
            gram=fixed[m].copy();gram[-1,-1]=float(np.sum(core*core));gram[0,-1]=gram[-1,0]=float(np.sum(factor*factor))
            for k,z in enumerate(traces[m],1):
                cross=z.T@factor;gram[k,-1]=gram[-1,k]=float(np.sum(cross*cross));del cross
            diagnostics=gram_diagnostics(gram);gp=root/(m+'-'+tree+'-kernel-gram.npz');np.savez_compressed(gp,gram=gram)
            allgrams.append(dict(mode=m,tree=tree,path=gp.name,sha256=sha(gp),**diagnostics))
            started=time.perf_counter();operator=SharedEntityCovariance(blocklabels,{**{k:operators[m,k] for k in KINDS},'family_intercept':family},factor,np.ones(n),v,.25)
            solved=operator.solve(rhS);residual=float(np.max(abs(operator.apply(solved)-rhS)));elapsed=time.perf_counter()-started
            out=folder/(m+'-'+tree+'.npz');np.savez_compressed(out,solution=solved)
            benchmarks.append(dict(mode=m,tree=tree,path=str(out.relative_to(root)),sha256=sha(out),entity_variances=v,species_variance=.25,residual_diagonal=1.,
                logdet=operator.logdet,maximum_solve_residual=residual,solve_diagnostic=operator.last_solve_diagnostic,elapsed_seconds=elapsed))
            print('full_entity_operator_numeric_benchmark',m,tree,'cases',n,'seconds',elapsed,flush=True)
            del operator
        del factor
    (root/'kernel_gram_manifest.json').write_text(json.dumps(allgrams,indent=2)+'\n')
    (root/'benchmark_manifest.json').write_text(json.dumps(benchmarks,indent=2)+'\n')
    summary=dict(logical_cases=n,entity_occurrences=sum(counts.values()),unique_entities={k:len(v) for k,v in labelmap.items()},matrices=len(entries),
        operator_nnz={r['mode']+'|'+r['kind']:r['nnz'] for r in entries},family_components=len(parts),largest_component=max(len(p) for p in parts),
        sum_component_squared_sizes=sum(len(p)**2 for p in parts),trees=plan['trees'],kernel_gram_ranks={r['mode']+'|'+r['tree']:r['normalized_kernel_gram_rank'] for r in allgrams},
        zero_kernels={r['mode']+'|'+r['tree']:r['zero_kernels'] for r in allgrams},numerical_benchmarks=len(benchmarks))
    assert summary['matrices']==13 and summary['family_components']==plan['expected']['family_components']
    verify(bindings)
    names=['stage_plan.json','operator_manifest.json','entity_labels.json','case_ids.json','block_labels.npy','kernel_gram_manifest.json','benchmark_manifest.json']+[r['path'] for r in entries+allgrams+benchmarks]
    receipt=dict(status='complete_full_entity_operator_bank_pending_independent_readback',plan_sha256=sha(path),**summary,source_hashes=bindings,
        artifacts={name:sha(root/name) for name in names},scientific_eligibility=False,scope=plan['scope'])
    with (root/'receipt.json').open('x') as f:f.write(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(summary),flush=True);return receipt


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);run(p.parse_args().plan)
