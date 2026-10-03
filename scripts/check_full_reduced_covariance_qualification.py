#!/usr/bin/env python3
"""Check dense covariance equivalence, inherited envelopes and full export contracts."""
import argparse
from copy import deepcopy
from datetime import datetime, timezone
import csv
import gzip
import itertools
import json
from pathlib import Path

import numpy as np
from scipy import sparse, linalg

from ancestral_chain_attempt import sha
from covariance_basis_audit import audit_covariance_bases, _diagnostics
from covariance_exact_folds_v2 import certificate, variance_map
from full_covariance_qualification_sources import LINK_EXTRA
from full_expanded_model_design_sources import SETTING_FIELDS, digest
from full_reduced_covariance_sources import load
from reduced_covariance_basis import OLD_NAMES, reduce_record, readback_record, independent_diagnostics
from run_full_reduced_covariance_qualification import run, SUMMARY_FIELDS


def write(path, value):
    with Path(path).open('x') as f: f.write(json.dumps(value, indent=2, allow_nan=False) + '\n')


def rejected(function):
    try: function()
    except (AssertionError, ValueError, ArithmeticError, KeyError, FileExistsError, FileNotFoundError): return
    raise AssertionError('Malformed covariance evidence accepted')


def fixture(pair_exception, mode, seed):
    rng = np.random.default_rng(seed); n = 40
    labels = np.repeat(np.arange(8), 5)
    target = sparse.eye(n, format='csr')
    background = sparse.csr_matrix((np.ones(n), (np.arange(n), np.arange(n)//2)), shape=(n, n))
    # Keep every entity column within one family component.
    background = sparse.csr_matrix((np.ones(n), (np.arange(n), labels*3 + (np.arange(n)%5)//2)), shape=(n, n))
    pair_target = target.copy()
    if pair_exception:
        indices = np.arange(n); indices[1] = indices[0]
        pair_target = sparse.csr_matrix((np.ones(n), (np.arange(n), indices)), shape=(n, n))
    pair = sparse.hstack([pair_target, background], format='csr')
    gene = .5 * sparse.hstack([target, background, target, background], format='csr')
    model = .5 * sparse.hstack([pair, pair], format='csr')
    intercept = sparse.csr_matrix((np.ones(n), (np.arange(n), labels)), shape=(n, 8))
    family = sparse.csr_matrix(intercept.shape) if mode == 'signed' else 2 * intercept
    operators = dict(target_node=target, background_node=background, model_pair=pair, gene=gene, model=model, family=family)
    c = certificate(operators, intercept, mode, np.arange(n))
    cert = dict(cohort_id=digest(['synthetic', pair_exception]), certificate=c, variance_map=variance_map(c['relations'], mode),
        cohort_rows_sha256=digest(list(range(n))), ordered_case_ids_sha256=digest(['synthetic-case', n]),
        raw_reml_basis_qualification_complete=False, scientific_eligibility=False)
    assert len(cert['variance_map']['retained_names']) == (5 if pair_exception else 4)
    factor = rng.normal(size=(n, 6)); x = np.column_stack([np.ones(n), rng.normal(size=(n, 2))])
    incidence = {k: operators[k] for k in ['background_node', 'model_pair', 'gene', 'model']}
    incidence['family_intercept'] = intercept
    value = audit_covariance_bases(labels, incidence, factor, np.ones(n), x)
    old = dict(audit_id=digest(['synthetic-old', pair_exception, mode, seed]), source_contract='synthetic-closed-seven-kernel-source',
        cohort_id=cert['cohort_id'], design_id=digest(['synthetic-design', seed]), loading_mode=mode, tree='mafft_guide',
        residual_diagonal='uniform_one', records=n, source_design_disposition='full_rank_design',
        folded_terms={'target_node':'combined_with_uniform_residual','family':mode}, active_column_indices=[0,1,2],
        raw_design_sha256=digest(x.tolist()), disposition='covariance_basis_requires_review', numerical_audit=value)
    kernels = [np.eye(n)] + [(z @ z.T).toarray() for z in incidence.values()] + [factor @ factor.T]
    q = linalg.qr(x, mode='economic')[0]; h = np.eye(n) - q @ q.T
    direct_raw = [[float(np.sum(a*b)) for b in kernels] for a in kernels]
    projected = [h @ a @ h for a in kernels]
    direct_reml = [[float(np.sum(a*b)) for b in projected] for a in projected]
    np.testing.assert_allclose(value['raw_gram'], direct_raw, rtol=1e-12, atol=1e-10)
    np.testing.assert_allclose(value['projected_gram'], direct_reml, rtol=1e-12, atol=1e-10)
    original_kernels = [np.eye(n), np.eye(n), *[(operators[k]@operators[k].T).toarray()
        for k in ['background_node','model_pair','gene','model','family']], (intercept@intercept.T).toarray(), factor@factor.T]
    retained = [kernels[OLD_NAMES.index(k)] for k in cert['variance_map']['retained_names']]
    for _ in range(20):
        theta = rng.random(9); eta = np.asarray(cert['variance_map']['forward']) @ theta
        np.testing.assert_allclose(sum(t*k for t,k in zip(theta,original_kernels)), sum(t*k for t,k in zip(eta,retained)), rtol=1e-13, atol=1e-13)
    return old, cert


def closed(root, name, status, paths, summary):
    hashes = {str(p): sha(p) for p in paths}
    archive = root / (name + '-archive.json')
    write(archive, dict(status=status+'_archive', source_hashes=hashes, summary=summary,
        services=[{'synthetic_original_journal_fixture':True}]*2))
    completion = root / (name + '-completed.json')
    write(completion, dict(status=status, **summary, exact_process_journals_checked=2,
        scientific_eligibility=False, full_hash_archive=str(archive), full_hash_archive_sha256=sha(archive), bound_source_hashes=len(hashes)))
    return completion


def full_fixture(root, fixtures, trees):
    oldroot = root / 'original'; exactroot = root / 'exact'; oldroot.mkdir(); exactroot.mkdir()
    oldplan = root / 'original-plan.json'; exactplan = root / 'exact-plan.json'
    write(oldplan, dict(output=str(oldroot))); write(exactplan, dict(output=str(exactroot)))
    for folder in [oldroot, exactroot]:
        for name in ['receipt.json', 'readback.json']: write(folder/name, {'synthetic_closed_parent':True})
    write(oldroot/'stage_plan.json', {'synthetic_closed_parent':True})
    certificate_rows = [fixtures[p,m][1] for p,m in itertools.product([False,True],['signed','unsigned'])]
    with (exactroot/'cohort_certificates.jsonl').open('x') as f:
        for c in certificate_rows: f.write(json.dumps(c)+'\n')
    designs = []; audits = []
    for exception in [False,True]:
        for d in range(30):
            did = digest(['full-synthetic-design', exception, d]); designs.append((exception, did))
            for mode, tree in itertools.product(['signed','unsigned'], trees):
                old = deepcopy(fixtures[exception,mode][0]); old.update(design_id=did, tree=tree, audit_id=digest(['old-synthetic-audit',did,mode,tree]))
                audits.append(old)
    with gzip.open(oldroot/'design_covariance_audits.jsonl.gz','wt') as f:
        for a in audits: f.write(json.dumps(a)+'\n')
    index = {(a['design_id'],a['loading_mode'],a['tree']): a for a in audits}
    with gzip.open(oldroot/'setting_audit_links.tsv.gz','wt') as f:
        writer = csv.DictWriter(f,SETTING_FIELDS+LINK_EXTRA,delimiter='\t',lineterminator='\n');writer.writeheader()
        for exception,did in designs:
            for outcome in ['rmsd_delta','native_tm_dissimilarity_delta']:
                setting = {k:'synthetic_setting' for k in SETTING_FIELDS}
                setting.update(design_id=did,cohort_id=fixtures[exception,'signed'][1]['cohort_id'],
                    fit_input_id=digest([did,outcome]),outcome=outcome,records='40',nominal_tree_fits='5',
                    disposition='constant_response' if outcome=='rmsd_delta' else 'ready_for_working_covariance_fit')
                for mode,tree in itertools.product(['signed','unsigned'],trees):
                    old = index[did,mode,tree]
                    writer.writerow({**setting,'loading_mode':mode,'tree':tree,'audit_id':old['audit_id'],
                        'covariance_disposition':old['disposition'],'combined_disposition':setting['disposition']
                        if setting['disposition']!='ready_for_working_covariance_fit' else old['disposition']})
    expected=dict(logical_cases=80,model_setting_rows=120,unique_cohorts=2,unique_designs=60,audit_rows=600,setting_audit_links=1200)
    oldsummary=dict(**expected,trees=trees,loading_modes=['signed','unsigned'])
    qcompletion=closed(root,'original','complete_verified_full_uniform_covariance_qualification',
        [oldplan,*[oldroot/n for n in ['receipt.json','readback.json','stage_plan.json','design_covariance_audits.jsonl.gz','setting_audit_links.tsv.gz']]],oldsummary)
    ecompletion=closed(root,'exact','complete_verified_full_exact_uniform_covariance_folds_v2',
        [exactplan,*[exactroot/n for n in ['receipt.json','readback.json','cohort_certificates.jsonl']]],dict(logical_cases=80,cohorts=2,certificates=4))
    plan=dict(qualification_plan=str(oldplan),qualification_completion=str(qcompletion),exact_plan=str(exactplan),
        exact_completion=str(ecompletion),expected=expected,trees=trees,pins={},output=str(root/'new'),
        resources={'minimum_free_disk_gib':1},scope='Synthetic software fixture; parent journals are synthetic and no biological inference is performed.')
    path=root/'new-plan.json'; write(path,plan)
    producer=run(path); reader=run(path,reader=True)
    assert all(producer[k]==reader[k] for k in SUMMARY_FIELDS)
    assert producer['audit_rows']==600 and producer['setting_audit_links']==1200
    assert producer['retained_basis_audit_counts']=={'4':300,'5':300}
    assert producer['link_status_counts']['constant_response']==600
    assert producer['audit_status_counts']=={'qualified_exact_retained_uniform_covariance_basis':600}
    rejected(lambda:run(path)); rejected(lambda:run(path,reader=True))
    return path, qcompletion


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(exist_ok=False);assert not a.receipt.exists()
    fixtures={}; cases=[]; malformed=[]
    for exception,mode in itertools.product([False,True],['signed','unsigned']):
        old,cert=fixture(exception,mode,100+int(exception)); fixtures[exception,mode]=(old,cert)
        reduced=reduce_record(old,cert,'synthetic-contract');readback_record(reduced,old,cert,'synthetic-contract')
        assert reduced['disposition']=='qualified_exact_retained_uniform_covariance_basis'
        cases.append(dict(pair_exception=exception,mode=mode,retained=len(reduced['retained_kernel_names'])))
    # Numerically unresolved and boundary cases retain explicit review states.
    for gram,error in [(np.diag([1.,1.,0.,1.]),np.eye(4)*1e-10),
        (np.ones((4,4)),np.ones((4,4))*1e-8),
        (np.eye(4),np.eye(4)*.2)]:
        names=['a','b','c','d'];recorded=_diagnostics(gram,error,names);independent=independent_diagnostics(gram,error,names)
        for k in ['disposition','rank','rank_at_one_tenth_tolerance','rank_at_ten_times_tolerance']:
            assert recorded[k]==independent[k]
        assert recorded['disposition']!='numerically_independent_covariance_bases'
        cases.append(dict(numerical_review=recorded['disposition']))
    old,cert=fixtures[True,'unsigned'];saved=reduce_record(old,cert,'synthetic-contract')
    def altered(name,change):
        row=deepcopy(saved);change(row);rejected(lambda:readback_record(row,old,cert,'synthetic-contract'));malformed.append(name)
    altered('changed_principal_gram',lambda r:r['numerical_audit']['raw_gram'][0].__setitem__(0,-1))
    altered('narrowed_inherited_error',lambda r:r['numerical_audit']['raw_roundoff_envelope'][0].__setitem__(0,0))
    altered('lost_fifth_kernel',lambda r:r['retained_kernel_names'].remove('model_pair'))
    altered('changed_source_audit_identity',lambda r:r.__setitem__('source_audit_id','foreign'))
    altered('changed_fixed_design',lambda r:r.__setitem__('raw_design_sha256','foreign'))
    altered('promoted_scientific_eligibility',lambda r:r.__setitem__('scientific_eligibility',True))
    altered('promoted_component_attribution',lambda r:r.__setitem__('component_variance_attribution_accepted',True))
    altered('changed_original_status',lambda r:r.__setitem__('source_covariance_disposition','qualified_uniform_working_covariance_basis'))
    altered('false_numerical_rank',lambda r:r['numerical_audit']['reml_diagnostics'].__setitem__('rank',0))
    changed=deepcopy(cert);changed['variance_map']['forward'][0][0]=-1
    rejected(lambda:reduce_record(old,changed,'synthetic-contract'));malformed.append('negative_variance_map')
    changed=deepcopy(old);changed['residual_diagonal']='nonuniform'
    rejected(lambda:reduce_record(changed,cert,'synthetic-contract'));malformed.append('nonuniform_residual')
    changed=deepcopy(old);changed['cohort_id']='foreign'
    rejected(lambda:reduce_record(changed,cert,'synthetic-contract'));malformed.append('foreign_cohort')
    for status in ['empty_design','rank_deficient_design','numerical_covariance_qualification_requires_review']:
        changed=deepcopy(old);changed.update(disposition=status,numerical_audit=None)
        reduced=reduce_record(changed,cert,'synthetic-contract');readback_record(reduced,changed,cert,'synthetic-contract')
        assert reduced['disposition']==status and reduced['numerical_audit'] is None
        bad=deepcopy(reduced);bad['disposition']='qualified_exact_retained_uniform_covariance_basis'
        rejected(lambda:readback_record(bad,changed,cert,'synthetic-contract'))
        malformed.append('promoted_'+status)
    trees=['mafft_guide','pmsf_mafft_profile','pmsf_profile_mafft','pmsf_profile_profile','profile_guide']
    path,completion=full_fixture(a.output,fixtures,trees)
    plan=json.loads(path.read_text());parent=json.loads(completion.read_text())
    # A source cannot be consumed before its original independent closure.
    completion.write_text(json.dumps({**parent,'status':'pending_independent_readback'}))
    rejected(lambda:load(plan,path));completion.write_text(json.dumps(parent,indent=2)+'\n');malformed.append('unclosed_original_source')
    real=Path('metadata/full_exact_covariance_folds_completed_20261003_v2.json');summary=json.loads(real.read_text())
    assert summary['status']=='complete_verified_full_exact_uniform_covariance_folds_v2'
    assert (summary['cohorts'],summary['certificates'],summary['logical_cases'])==(4340,8680,75188)
    from full_exact_covariance_sources import closed_subset
    bindings={};cfg=Path('metadata/full_exact_covariance_folds_plan_20261003_v2.json');root=Path(json.loads(cfg.read_text())['output'])
    closed_subset(real,summary['status'],[cfg,root/'cohort_certificates.jsonl'],bindings)
    from full_covariance_qualification_sources import jsonl
    realcerts=list(jsonl(root/'cohort_certificates.jsonl'))
    from reduced_covariance_basis import validate_certificate
    counts={}
    for c in realcerts:
        n=len(validate_certificate(c));counts[str(n)]=counts.get(str(n),0)+1
    assert counts=={'4':7232,'5':1448}
    own=['check_full_reduced_covariance_qualification','reduced_covariance_basis','full_reduced_covariance_sources',
        'run_full_reduced_covariance_qualification','covariance_exact_folds_v2','full_exact_covariance_sources','covariance_basis_audit']
    for name in own:bindings['scripts/'+name+'.py']=sha('scripts/'+name+'.py')
    artifacts={str(p):sha(p) for p in a.output.rglob('*') if p.is_file()}
    result=dict(status='passed_full_exact_retained_covariance_qualification_software_contracts',checked_utc=datetime.now(timezone.utc).isoformat(),
        dense_and_review_scenarios=cases,dense_nonnegative_cone_roundtrips=80,malformed_cases_rejected=malformed,
        synthetic_pipeline_audits=600,synthetic_pipeline_setting_links=1200,all_five_trees_and_both_modes=True,
        inherited_numerical_envelopes_unchanged=True,nonready_setting_states_retained=True,completed_restart_refused=True,
        source_and_journal_fixtures_synthetic=True,full_real_certificates_checked=8680,full_real_cohorts=4340,
        full_real_logical_cases=75188,retained_real_certificate_counts=counts,source_hashes=bindings,artifacts=artifacts,
        scientific_eligibility=False,scope='Software and full closed certificate checks. Direct dense raw/REML products, covariance-cone equality, unchanged envelopes, independent gesvd, every synthetic audit/link and explicit malformed cases; no biological pilot or fit. Full production qualification remains gated on closed original full-grid arithmetic.')
    write(a.receipt,result);print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']},indent=2))


if __name__=='__main__':main()
