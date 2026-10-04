#!/usr/bin/env python3
"""Dense and original-model contracts for dependent Gaussian calibration draws."""
import argparse
from copy import deepcopy
from datetime import datetime,timezone
import itertools
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.stats import binomtest

from ancestral_chain_attempt import sha
from check_retained_shared_entity_candidates import inputs
from check_shared_entity_covariance import dense_covariance
from full_expanded_model_design_sources import digest
from full_weighted_fit_admission import FULL
from full_weighted_fit_exports import atomic
from reference_measurement_union_sources import bind,verify
from shared_entity_covariance import SharedEntityCovariance
from weighted_shared_entity_simulation import GaussianSharedEntitySimulation,coverage_counts,replay_response


def rejected(call):
    try:call()
    except (ValueError,ArithmeticError,AssertionError,TypeError,KeyError):return
    raise AssertionError('Invalid simulation contract accepted')


def reference_response(model,mean,scale,ratios,latents):
    # Deliberately explicit latent columns, independent of generator's sparse
    # matrix products. Signed/coalesced entries and shared columns stay intact.
    response=mean.copy()+np.sqrt(scale*model.diagonal)*latents['residual']
    for j,(name,z) in enumerate(model.incidence.items()):
        for k in range(z.shape[1]):response+=np.sqrt(scale*ratios[j])*z[:,k].toarray().ravel()*latents[name][k]
    for k in range(model.factor.shape[1]):response+=np.sqrt(scale*ratios[-1])*model.factor[:,k]*latents['species'][k]
    return response


def model_for(case):
    source,rows,x,y,ops,*_=inputs(case['pair_exception'],case['loading_mode'])
    if case['policy']!='uniform':ops={'target_node':sparse.eye(len(rows),format='csr'),**ops}
    model=GaussianSharedEntitySimulation(source['labels'],ops,source['factors']['mafft_guide'],case['diagonal'],x,case['candidate']['fit']['parameter_names'])
    assert ['residual',*model.parameter_names]==case['audit']['retained_kernel_names']
    return model


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True);a=p.parse_args()
    assert not a.receipt.exists();a.output.mkdir(parents=True,exist_ok=False)
    gate=Path('metadata/weighted_shared_entity_candidate_software_validation_20261004_v1.json')
    parent=json.loads(gate.read_text());verify(parent['source_hashes']);bindings=dict(parent['source_hashes']);bind(bindings,gate)
    case_path=Path('data/software_audits/weighted-shared-entity-candidates-20261004-v1/actual_weighted_candidate_cases.json')
    cases=json.loads(case_path.read_text());assert len(cases)==64;bind(bindings,case_path)
    full_plan=Path('metadata/full_weighted_shared_entity_fit_draft_plan_20261004_v1.json')
    full=json.loads(full_plan.read_text());verify(full['pins']);assert full['expected']==FULL;bind(bindings,full_plan)
    for name in ['weighted_shared_entity_simulation','check_weighted_shared_entity_simulation','matched_calibration_intervals',
        'shared_entity_covariance','check_shared_entity_covariance','check_retained_shared_entity_candidates']:bind(bindings,'scripts/'+name+'.py')
    dimensions=set();all_records=[];mean_error=covariance_error=response_error=0.;empirical=[]
    for number,case in enumerate(cases):
        model=model_for(case);dimensions.add(1+len(model.parameter_names));beta=np.array([.2,-.3,.4])
        source_id=case['candidate']['candidate_id'];records=[];responses=[]
        for scenario,scale in [('zero_components',.25),('positive_components',1.),('sparse_boundary',2.)]:
            ratios=(np.zeros(len(model.parameter_names)) if scenario=='zero_components' else
                np.full(len(model.parameter_names),.7) if scenario=='positive_components' else
                np.array([1e4]+[0.]*(len(model.parameter_names)-1)))
            mean,residual,components,description=model.scenario(beta,scale,ratios)
            expected=scale*dense_covariance(model.diagonal,model.incidence,model.factor,dict(zip(model.incidence,ratios[:-1])),ratios[-1])
            covariance=SharedEntityCovariance(np.repeat(np.arange(8),5),model.incidence,model.factor,model.diagonal,dict(zip(model.incidence,ratios[:-1])),ratios[-1])
            actual=scale*covariance.apply(np.eye(model.n));np.testing.assert_allclose(actual,expected,rtol=2e-12,atol=2e-11)
            covariance_error=max(covariance_error,float(np.max(abs(actual-expected))))
            for replicate in range(3):
                response,record,latents=model.draw(beta,scale,ratios,4002026,source_id,scenario,replicate,keep_latents=True)
                independent=reference_response(model,mean,scale,ratios,latents)
                np.testing.assert_allclose(response,independent,rtol=2e-12,atol=2e-11)
                response_error=max(response_error,float(np.max(abs(response-independent))))
                np.testing.assert_array_equal(response,replay_response(model,record))
                assert record['scientific_eligibility'] is False and record['fits_computed']==0
                records.append(record);responses.append(response)
        # Response order and worker order do not determine the random stream.
        for record,response in reversed(list(zip(records,responses))):np.testing.assert_array_equal(replay_response(model,record),response)
        assert len({r['response_sha256'] for r in records})==9
        atomic(a.output/('case-'+str(number)+'.json'),dict(original_saved_case=number,software_fixture_only=True,records=records))
        np.savez_compressed(a.output/('case-'+str(number)+'.npz'),responses=np.asarray(responses),design=model.design,
            diagonal=model.diagonal,factor=model.factor,beta=beta)
        all_records+=records
        if case['outcome']=='rmsd_delta' and case['method']=='ml':
            # Seeded Gaussian moment smoke checks supplement exact algebra;
            # this is not an empirical coverage or biological calibration run.
            ratios=np.full(len(model.parameter_names),.7);m=512
            draws=np.asarray([model.draw(beta,1.,ratios,4012026,source_id,'moment-smoke',r)[0] for r in range(m)])
            expected=dense_covariance(model.diagonal,model.incidence,model.factor,dict(zip(model.incidence,ratios[:-1])),ratios[-1])
            error_mean=np.max(abs(draws.mean(axis=0)-model.design@beta)/np.sqrt(np.diag(expected)/m))
            se=np.sqrt((expected**2+np.outer(np.diag(expected),np.diag(expected)))/(m-1))
            error_cov=np.max(abs(np.cov(draws,rowvar=False)-expected)/se)
            assert error_mean<8 and error_cov<8
            empirical.append(dict(case=number,replicates=m,maximum_mean_standard_error_units=float(error_mean),maximum_covariance_standard_error_units=float(error_cov)))
    assert dimensions=={4,5,6} and len(all_records)==576 and len(empirical)==16
    axes={(c['pair_exception'],c['loading_mode'],c['policy'],c['outcome'],c['method']) for c in cases}
    assert len(axes)==64 and len({r['response_sha256'] for r in all_records})==576
    atomic(a.output/'moment-smoke.json',dict(cases=empirical,software_only=True))
    # Signed duplicates cancel before drawing. Species rank zero and explicit
    # zero operators are retained; target I stays distinct from nonuniform D.
    labels=['A']*3+['B']*3
    z=sparse.coo_matrix(([.5,-.5,1.,-1.,1.,-1.],([0,0,1,2,3,4],[0,0,0,0,1,1])),shape=(6,2))
    ops={'target_node':sparse.eye(6,format='csr'),'signed_shared':z,'zero':sparse.csr_matrix((6,0))}
    d=np.linspace(.5,2.,6);x=np.column_stack([np.ones(6),np.arange(6)])
    special=GaussianSharedEntitySimulation(labels,ops,np.zeros((6,0)),d,x,[*ops,'species'])
    np.testing.assert_array_equal(special.incidence['signed_shared'].toarray()[0],np.zeros(2))
    ratios=np.array([.3,2.,1.,4.]);response,r,latents=special.draw([.1,.2],1.,ratios,40,'signed-software','zero-rank',0,True)
    np.testing.assert_allclose(response,reference_response(special,x@np.array([.1,.2]),1.,ratios,latents),atol=1e-12)
    cov=dense_covariance(d,special.incidence,special.factor,dict(zip(special.incidence,ratios[:-1])),ratios[-1]);assert cov[1,2]<0
    atomic(a.output/'signed-cancellation-zero-rank.json',dict(record=r,covariance=cov.tolist(),software_fixture_only=True))
    bad=[];baseline=model_for(cases[0]);names=baseline.parameter_names
    for name in ['zero_diagonal','negative_diagonal','nan_diagonal','wrong_diagonal','nan_factor','nan_design','rank_design','zero_column',
        'missing_name','wrong_name_order','duplicate_name','entity_crosses_blocks','nan_operator','negative_ratio','nan_ratio','wrong_ratios',
        'wrong_beta','nan_beta','zero_scale','negative_scale','infinite_scale','boolean_scale','underflow_variance','overflow_variance',
        'empty_candidate','empty_scenario','negative_seed','boolean_seed','negative_replicate']:
        if name in ['zero_diagonal','negative_diagonal','nan_diagonal','wrong_diagonal','nan_factor','nan_design','rank_design','zero_column',
            'missing_name','wrong_name_order','duplicate_name','entity_crosses_blocks','nan_operator']:
            d=baseline.diagonal.copy();f=baseline.factor.copy();x=baseline.design.copy();op={k:z.copy() for k,z in baseline.incidence.items()};ns=names.copy()
            if name in ['zero_diagonal','negative_diagonal','nan_diagonal']:d[0]={'zero_diagonal':0.,'negative_diagonal':-1.,'nan_diagonal':np.nan}[name]
            elif name=='wrong_diagonal':d=d[:-1]
            elif name=='nan_factor':f[0,0]=np.nan
            elif name=='nan_design':x[0,0]=np.nan
            elif name=='rank_design':x[:,1]=x[:,0]
            elif name=='zero_column':x[:,1]=0.
            elif name=='missing_name':ns.pop()
            elif name=='wrong_name_order':ns.reverse()
            elif name=='duplicate_name':ns[1]=ns[0]
            elif name=='entity_crosses_blocks':op[names[0]]=sparse.csr_matrix(np.ones((40,1)))
            else:op[names[0]]=sparse.csr_matrix(np.full((40,1),np.nan))
            rejected(lambda:GaussianSharedEntitySimulation(np.repeat(np.arange(8),5),op,f,d,x,ns))
        else:
            b=[.2,-.3,.4];ratios=np.ones(len(names));scale=1.;seed=40;candidate='software';scenario='invalid';replicate=0
            if name=='negative_ratio':ratios[0]=-1
            elif name=='nan_ratio':ratios[0]=np.nan
            elif name=='wrong_ratios':ratios=ratios[:-1]
            elif name=='wrong_beta':b=b[:-1]
            elif name=='nan_beta':b[0]=np.nan
            elif name in ['zero_scale','negative_scale','infinite_scale','boolean_scale']:scale={'zero_scale':0.,'negative_scale':-1.,'infinite_scale':np.inf,'boolean_scale':True}[name]
            elif name=='underflow_variance':scale=1e-200;ratios[0]=1e-200
            elif name=='overflow_variance':scale=1e300;ratios[0]=1e300
            elif name=='empty_candidate':candidate=''
            elif name=='empty_scenario':scenario=''
            elif name=='negative_seed':seed=-1
            elif name=='boolean_seed':seed=True
            else:replicate=-1
            rejected(lambda:baseline.draw(b,scale,ratios,seed,candidate,scenario,replicate))
        bad.append(name)
    # Never discard unresolved replicates or combine scenarios. Coverage states
    # below are synthetic accounting fixtures, with no optimizer/refit claim.
    accounting=[]
    for i in range(4):
        _,record,_=baseline.draw([.2,-.3,.4],1.,np.ones(len(names)),45,'accounting-software','fixed-four',i)
        if i<2:
            record['status']='refit_independently_checked_pending_calibration';record['nominal_interval_covers_generating_beta']=[[True,False,True],[False,True,True]][i]
        elif i==3:record['status']='refit_error_requires_review'
        accounting.append(record)
    counted=coverage_counts(accounting,3,4);assert counted['attempted']==4 and counted['unresolved']==2 and counted['covered_qualified']==[1,1,2]
    for covered,interval in zip(counted['covered_qualified'],counted['marginal_monte_carlo_intervals']):
        intervals=[binomtest(k,4).proportion_ci(method='exact') for k in range(covered,covered+3)]
        np.testing.assert_allclose([interval['lower'],interval['upper']],[min(v.low for v in intervals),max(v.high for v in intervals)],rtol=1e-10,atol=1e-12)
    account_bad=[]
    for name in ['missing_replicate','duplicate_replicate','boolean_replicate','foreign_candidate','foreign_scenario','foreign_input','foreign_truth','promoted_science','unknown_status','nonboolean_coverage','missing_coefficient']:
        records=deepcopy(accounting)
        if name=='missing_replicate':records.pop()
        elif name=='duplicate_replicate':records[3]['replicate']=0
        elif name=='boolean_replicate':records[0]['replicate']=False
        elif name in ['foreign_candidate','foreign_scenario','foreign_input','foreign_truth']:records[0][{'foreign_candidate':'candidate_id','foreign_scenario':'scenario_id','foreign_input':'input_contract','foreign_truth':'scenario_contract'}[name]]='foreign'
        elif name=='promoted_science':records[0]['scientific_eligibility']=True
        elif name=='unknown_status':records[0]['status']='unknown'
        elif name=='nonboolean_coverage':records[0]['nominal_interval_covers_generating_beta'][0]=1
        else:records[0]['nominal_interval_covers_generating_beta'].pop()
        rejected(lambda:coverage_counts(records,3,4));account_bad.append(name)
    response,saved,_=baseline.draw([.2,-.3,.4],1.,np.ones(len(names)),45,'replay-software','fixed',0)
    replay_bad=[]
    for name in ['response_hash','changed_seed','changed_truth','changed_component_stream','changed_names','promoted_science','invented_refit']:
        record=deepcopy(saved)
        if name=='response_hash':record['response_sha256']='foreign'
        elif name=='changed_seed':record['master_seed']+=1
        elif name=='changed_truth':record['scenario']['beta'][0]+=1
        elif name=='changed_component_stream':record['component_seed_entropy']['species'][0]+=1
        elif name=='changed_names':record['scenario']['parameter_names'].reverse()
        elif name=='promoted_science':record['scientific_eligibility']=True
        else:record['fits_computed']=1
        rejected(lambda:replay_response(baseline,record));replay_bad.append(name)
    atomic(a.output/'coverage-accounting.json',dict(records=accounting,summary=counted,synthetic_accounting_only=True))
    for q in a.output.rglob('*'):
        if q.is_file():bind(bindings,q)
    verify(bindings)
    result=dict(status='passed_weighted_shared_entity_gaussian_simulation_and_fixed_sample_accounting_v1',checked_utc=datetime.now(timezone.utc).isoformat(),
        original_software_model_cases=64,basis_dimensions=sorted(dimensions),exact_latent_response_replays=576,
        dense_covariance_cases=192,maximum_absolute_dense_covariance_error=covariance_error,maximum_absolute_latent_response_error=response_error,
        gaussian_moment_smoke_cases=16,gaussian_moment_smoke_replicates=8192,signed_duplicate_cancellation_and_zero_rank_passed=True,
        input_and_scenario_rejections=bad,coverage_accounting_rejections=account_bad,response_replay_rejections=replay_bad,
        scheduling_order_independent=True,original_full_expected=FULL,full_data_simulation_launched=False,
        fits_computed=0,scientific_eligibility=False,gpu=False,new_cost_usd=0,source_hashes=bindings,
        scope='Software Gaussian-response generator and fixed-sample unresolved-coverage accounting for all64original synthetic weighted models, q4/q5/q6 and all control/mode/outcome/method axes. No refit mocked as actual fit; coverage statuses are explicitly synthetic accounting fixtures. Actual original closed source/fit admission, full scenarios/refits/calibration, global dependence/multiple-testing and adequacy remain required; full biological scope unchanged, all eight aims incomplete.')
    atomic(a.receipt,result)
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
