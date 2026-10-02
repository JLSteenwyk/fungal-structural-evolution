"""Closed source and schema definitions for every expanded entity operator."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
from background_measurement_union_sources import closed_source
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha

KINDS=['target_node','background_node','model_pair','family','gene','model']
MODES=['signed','unsigned']
KERNELS=['residual_identity',*KINDS,'family_intercept','species']
SUMMARY_FIELDS=['logical_cases','entity_occurrences','unique_entities','matrices','operator_nnz',
    'family_components','largest_component','sum_component_squared_sizes','trees',
    'kernel_gram_ranks','zero_kernels','numerical_benchmarks']


def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def table(path):
    with gzip.open(path,'rt') as f:return list(csv.DictReader(f,delimiter='\t'))


def load(plan,path):
    bindings=dict(plan['pins']);bind(bindings,path)
    c=closed_source(plan['covariance_completion'],'complete_verified_full_expanded_covariance',
        'complete_verified_full_expanded_covariance_archive',2,bindings)
    config=json.loads(Path(plan['covariance_plan']).read_text());bind(bindings,plan['covariance_plan']);root=Path(config['output'])
    rp=root/'receipt.json';r=json.loads(rp.read_text())
    assert c['producer_receipt']==str(rp) and c['producer_receipt_sha256']==bindings[str(rp)]==sha(rp)
    assert r['status']=='complete_full_expanded_covariance_pending_independent_readback' and r['plan_sha256']==sha(plan['covariance_plan'])
    assert r['scientific_eligibility'] is c['scientific_eligibility'] is False
    for k in ['logical_cases','entity_occurrences','family_components']:assert c[k]==r[k]==plan['expected'][k]
    assert c['trees']==config['trees']==plan['trees']
    assert set(c['unique_entities'])==set(KINDS)
    for name,d in r['artifacts'].items():bind(bindings,root/name,d)
    del r
    cases=table(root/'case_covariance_index.tsv.gz')
    assert len(cases)==len({r['case_id'] for r in cases})==c['logical_cases']
    assert all(r['target_family']==r['background_family'] for r in cases)
    assert all(len(r['case_id'])==len(r['family_component'])==64 for r in cases)
    verify(bindings)
    return dict(root=root,cases=cases,completion=c,incidence=root/'entity_incidence.tsv.gz'),bindings


def gram_diagnostics(matrix):
    diagonal=np.diag(matrix);assert np.isfinite(matrix).all() and np.all(diagonal>=0)
    active=np.flatnonzero(diagonal>0);normalized=matrix[np.ix_(active,active)]/np.sqrt(np.outer(diagonal[active],diagonal[active]))
    singular=np.linalg.svd(normalized,compute_uv=False)
    tolerance=len(active)*np.finfo(float).eps*(float(singular[0]) if len(singular) else 0)
    return dict(kernel_names=KERNELS,zero_kernels=[KERNELS[i] for i in np.flatnonzero(diagonal==0)],
        active_kernel_indices=active.tolist(),normalized_gram_singular_values=singular.tolist(),rank_tolerance=tolerance,
        normalized_kernel_gram_rank=int((singular>tolerance).sum()),
        scope='Full unfiltered-case Frobenius Gram diagnostic. Gram conditioning squares covariance-basis conditioning; cohort-specific identifiability and variance fitting remain separate. No eigenvalue clipping or model qualification.')


def rhs(n):return np.column_stack([np.ones(n),np.linspace(-1,1,n),np.sin(np.arange(n,dtype=float))])


def benchmark_variances():return {k:(.5 if k=='family_intercept' else .125) for k in [*KINDS,'family_intercept']}
