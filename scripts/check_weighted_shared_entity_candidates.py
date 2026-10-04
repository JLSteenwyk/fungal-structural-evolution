#!/usr/bin/env python3
"""Qualify actual-D fits and complete closed-grid source/candidate identities."""
import argparse
from collections import Counter
from copy import deepcopy
from datetime import datetime, timezone
import itertools
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from ancestral_chain_attempt import sha
from check_full_reduced_covariance_qualification import closed, rejected
from check_full_shared_entity_fits_v2 import settings
from check_retained_shared_entity_candidates import inputs
from check_shared_entity_likelihood import dense
from full_expanded_model_design_sources import array_digest, digest
from full_weighted_covariance_qualification import SUMMARY, disposition
from full_weighted_covariance_sources_v2 import POLICIES
from full_weighted_shared_entity_fit_sources import load, cohorts, cases
from independent_positive_diagonal_basis_context import IndependentPositiveDiagonalBasisContext
from nonuniform_covariance_cone import record as cone_record
from positive_diagonal_basis_context import PositiveDiagonalBasisContext
from readback_full_covariance_qualification import numeric as audit_qualification
from readback_weighted_shared_entity_candidate import numeric
from reference_measurement_union_sources import verify
from shared_entity_likelihood import SharedEntityLikelihood
from weighted_shared_entity_candidate import READY, candidate, backend_guard


def write(path, value):
    with path.open('x') as f: json.dump(value, f, indent=2, allow_nan=False); f.write('\n')


def source_grids(root):
    # Reuse completed, frozen full synthetic grids without editing any parent.
    gate_path = Path('metadata/full_weighted_covariance_qualification_software_validation_20261004_v2.json')
    gate = json.loads(gate_path.read_text()); verify(gate['source_hashes'])
    plans = [Path('data/software_audits/full-weighted-covariance-qualification-20261004-v2/fixture/numerical.plan.json')]
    plans += [Path(g['source_plan']) for g in gate['additional_complete_qualified_grids']]
    exports = []; bindings = {str(gate_path):sha(gate_path)}; corruptions = []
    for number, qp in enumerate(plans):
        qplan = json.loads(qp.read_text()); qr = Path(qplan['output'])
        producer = json.loads((qr / 'receipt.json').read_text()); reader = json.loads((qr / 'readback.json').read_text())
        paths = [qp, qr / 'receipt.json', qr / 'readback.json', *[qr / n for n in producer['artifacts']]]
        summary = {k:reader[k] for k in SUMMARY}
        completion = closed(root, 'numerical-' + str(number),
            'complete_verified_full_four_control_covariance_numerical_qualification_v1', paths, summary)
        value = json.loads(completion.read_text()); value.update(producer_receipt=str(qr / 'receipt.json'),
            producer_receipt_sha256=sha(qr / 'receipt.json'), independent_readback=str(qr / 'readback.json'),
            independent_readback_sha256=sha(qr / 'readback.json'))
        completion.write_text(json.dumps(value, indent=2) + '\n')
        plan = dict(qualification_plan=str(qp), qualification_completion=str(completion),
            methods=['ml','reml'], loading_modes=['signed','unsigned'], policies=POLICIES, trees=qplan['trees'],
            pins={}, expected=dict(candidate_rows=24000, setting_fit_links=48000), **settings())
        pp = root / ('fit-plan-' + str(number) + '.json'); write(pp, plan)
        source, hashes = load(plan, pp); statuses = Counter(); identities = set(); count = 0
        policies = Counter(); dims = set(); outcomes = set()
        for cohort, rows, entries in cohorts(source, plan):
            for identity, x, y, audit, route, diagonal in cases(source, plan, cohort, entries):
                from weighted_shared_entity_candidate import validate_source
                validate_source(identity, audit, route, diagonal)
                assert identity['candidate_id'] not in identities; identities.add(identity['candidate_id'])
                assert identity['response_sha256'] == array_digest(y, '<f8')
                statuses[identity['source_combined_disposition']] += 1
                policies[identity['control_policy']] += 1; dims.add(len(route['names'])); outcomes.add(identity['outcome']); count += 1
        assert count == 24000 and set(policies.values()) == {6000} and len(outcomes) == 2
        if number > 0: assert statuses[READY] > 0 and dims == ({4,5} if number == 1 else {5,6})
        # Reject false closure receipts and missing consumed artifacts, even
        # with a newly hashed compact completion file. Parent bytes stay intact.
        original = completion.read_bytes()
        for name in ['producer_hash', 'reader_hash', 'wrong_status', 'changed_count']:
            altered = json.loads(original)
            if name == 'producer_hash': altered['producer_receipt_sha256'] = 'foreign'
            elif name == 'reader_hash': altered['independent_readback_sha256'] = 'foreign'
            elif name == 'wrong_status': altered['status'] = 'foreign'
            else: altered['numerical_audit_rows'] += 1
            completion.write_text(json.dumps(altered) + '\n')
            rejected(lambda: load(plan, pp)); corruptions.append([number, name])
            completion.write_bytes(original)
        fresh, fresh_hashes = load(plan, pp); assert fresh['fit_contract'] == source['fit_contract']
        assert hashes == fresh_hashes; verify(fresh_hashes); bindings.update(fresh_hashes)
        exports.append(dict(grid=number, source_plan=str(qp), candidate_rows=count,
            candidate_status_counts=dict(statuses), policy_counts=dict(policies), basis_dimensions=sorted(dims)))
    return exports, bindings, corruptions


def numerical(root):
    plan = settings(); records = []; primary_counts = Counter(); reader_counts = Counter(); rejects = []
    maximum_objective = maximum_gradient = 0.; dense_cases = exclusions = exhausted = 0
    for exception, mode, policy in itertools.product([False,True], ['signed','unsigned'], POLICIES):
        source, rows, x, y, ops, old_retained, old, cert = inputs(exception, mode); n = len(rows)
        diagonal = (np.ones(n) if policy == 'uniform' else
            np.where(np.arange(n) % 2, .5, 2.) if policy != 'family_component' else np.linspace(.5,2.,n))
        uniform = policy == 'uniform'
        proof = cert if uniform else cone_record(cert, 'synthetic-weighted-fit-cone')
        names = proof['variance_map']['retained_names']
        route = dict(names=names, certificate_sha256=digest(proof), exact_uniform_one=uniform,
            route='closed_exact_uniform_named_fold' if uniform else 'closed_positive_diagonal_cone')
        operators = {} if uniform else {'target_node':sparse.eye(n,format='csr')}
        operators.update(ops)
        qualification = PositiveDiagonalBasisContext(source['labels'], operators, source['factors'][old['tree']]).design(x).audit(diagonal)
        audit = dict(audit_id=digest(['synthetic-weighted', exception, mode, policy]),
            loading_mode=mode, tree=old['tree'], control_policy=policy, records=n,
            diagonal_sha256=array_digest(diagonal,'<f8'), control_record_sha256=digest(['synthetic-control',policy]),
            exact_certificate_sha256=route['certificate_sha256'], retained_kernel_names=names,
            basis_route=route['route'], residual_diagonal_is_exact_uniform_one=uniform,
            numerical_audit=qualification, disposition=disposition(qualification), scientific_eligibility=False,
            nonuniform_weighting_accepted=False, component_variance_attribution_accepted=False, covariance_model_selected=False)
        assert audit['disposition'] == READY
        snapshot = deepcopy(audit)
        guard = backend_guard(source,rows,x,operators,audit,diagonal)
        latent = IndependentPositiveDiagonalBasisContext(source['labels'],operators,source['factors'][old['tree']]).design(x).audit(diagonal)
        review,_ = audit_qualification(guard,latent,names,n); assert review is False
        for outcome, method in itertools.product(['rmsd_delta','native_tm_dissimilarity_delta'],['ml','reml']):
            response = y if outcome == 'rmsd_delta' else np.tanh(y) + .15*x[:,1]
            identity = dict(candidate_id=digest(['synthetic-fit',audit['audit_id'],outcome,method]),
                covariance_audit_id=audit['audit_id'], covariance_audit_sha256=digest(audit),
                source_fit_disposition='ready_for_working_covariance_fit', source_covariance_disposition=READY,
                source_combined_disposition=READY, method=method, outcome=outcome, response_sha256=array_digest(response,'<f8'),
                scientific_eligibility=False, nonuniform_weighting_accepted=False, component_variance_attribution_accepted=False)
            identity.update({k:audit[k] for k in ['loading_mode','tree','control_policy','records','diagonal_sha256',
                'control_record_sha256','exact_certificate_sha256','retained_kernel_names']})
            likelihood = SharedEntityLikelihood(source['labels'],operators,source['factors'][old['tree']],diagonal,x,response)
            for scale in [0., .7]:
                ratios=np.full(len(operators)+1,scale); value=likelihood.evaluate(ratios,method)
                reference=dense(source['labels'],operators,source['factors'][old['tree']],diagonal,x,response,ratios,method)
                for k in reference: np.testing.assert_allclose(value[k],reference[k],rtol=3e-9,atol=3e-10)
                maximum_objective=max(maximum_objective,abs(value['negative_profiled_likelihood']-reference['negative_profiled_likelihood']))
                maximum_gradient=max(maximum_gradient,float(np.max(abs(value['gradient']-reference['gradient']))));dense_cases+=1
            fitted=candidate(source,plan,rows,identity,x,response,operators,audit,route,diagonal)
            replay=numeric(source,plan,rows,identity,x,response,fitted,operators,audit,route,diagonal)
            assert fitted['fit'] is not None and 'variance_ratios' in fitted['fit'],fitted
            assert fitted['fit']['parameter_names']==names[1:]
            assert fitted['fit']['source_covariance_qualification']==guard
            primary_counts[fitted['disposition']]+=1;reader_counts[replay['disposition']]+=1
            records.append(dict(pair_exception=exception, policy=policy, loading_mode=mode, outcome=outcome,
                method=method, diagonal=diagonal.tolist(), audit=audit, candidate=fitted, independent=replay))
            low=deepcopy(plan);low['optimizer'].update(max_iterations=1,max_evaluations=1)
            failed=candidate(source,low,rows,identity,x,response,operators,audit,route,diagonal)
            failed_replay=numeric(source,low,rows,identity,x,response,failed,operators,audit,route,diagonal)
            assert failed_replay['disposition']=='original_numerical_failure_reproduced_requires_review';exhausted+=1
            constant=deepcopy(identity);constant.update(source_fit_disposition='constant_response_requires_review',
                source_combined_disposition='constant_response_requires_review',response_sha256=array_digest(np.ones(n),'<f8'))
            exported=candidate(source,plan,rows,constant,x,np.ones(n),operators,audit,route,diagonal)
            assert exported['numerical_attempted'] is False
            assert numeric(source,plan,rows,constant,x,np.ones(n),exported,operators,audit,route,diagonal)['disposition']=='source_review_or_exclusion_retained';exclusions+=1
            for name in ['objective','variance','start','backend_gram','invented_failure','false_weighting_acceptance']:
                bad=deepcopy(fitted)
                if name=='objective':bad['fit']['negative_profiled_likelihood']+=1.
                elif name=='variance':bad['fit']['variance_components'][0]+=1.
                elif name=='start':bad['fit']['starts'][1]['initial_coordinates'][0]+=.1
                elif name=='backend_gram':bad['fit']['source_covariance_qualification']['raw_gram'][0][0]+=1.
                elif name=='false_weighting_acceptance':bad['nonuniform_weighting_accepted']=True
                else:bad.update(fit=None,disposition='shared_entity_fit_error_requires_review',error_type='ArithmeticError',error_message='invented')
                rejected(lambda:numeric(source,plan,rows,identity,x,response,bad,operators,audit,route,diagonal))
                rejects.append([exception,mode,policy,outcome,method,name])
            print('actual_weighted_fit_case',len(records),'/64',flush=True)
        for name in ['changed_diagonal','uniform_substitution','lost_kernel','narrowed_envelope','false_science']:
            bad=deepcopy(audit);changed=diagonal.copy()
            if name=='changed_diagonal':changed[0]*=1.01
            elif name=='uniform_substitution':
                if uniform:changed[0]=2.
                else:changed[:]=1.
            elif name=='lost_kernel':bad['retained_kernel_names'].pop(1)
            elif name=='narrowed_envelope':bad['numerical_audit']['raw_roundoff_envelope'][0][0]=0.
            else:bad['scientific_eligibility']=True
            bad_identity=deepcopy(identity);bad_identity['covariance_audit_sha256']=digest(bad)
            rejected(lambda:candidate(source,plan,rows,bad_identity,x,response,operators,bad,route,changed))
            rejects.append([exception,mode,policy,name])
        # Dependent near-uniform and constant D remain actual source reviews.
        if not uniform:
            for d in [np.full(n,2.),1.+np.arange(n)*2.**-40]:
                bad=deepcopy(audit);bad['diagonal_sha256']=array_digest(d,'<f8')
                bad['numerical_audit']=PositiveDiagonalBasisContext(source['labels'],operators,source['factors'][old['tree']]).design(x).audit(d)
                bad['disposition']=disposition(bad['numerical_audit']);assert bad['disposition']!=READY
                bad_identity=deepcopy(identity);bad_identity.update(covariance_audit_sha256=digest(bad),diagonal_sha256=bad['diagonal_sha256'],
                    source_covariance_disposition=bad['disposition'],source_combined_disposition=bad['disposition'])
                excluded=candidate(source,plan,rows,bad_identity,x,response,operators,bad,route,d)
                assert excluded['numerical_attempted'] is False
                assert numeric(source,plan,rows,bad_identity,x,response,excluded,operators,bad,route,d)['disposition']=='source_review_or_exclusion_retained';exclusions+=1
        assert audit==snapshot
    assert len(records)==64 and dense_cases==128 and exhausted==64 and exclusions==88 and len(rejects)==464
    artifact=root/'actual_weighted_candidate_cases.json';write(artifact,records)
    return dict(actual_numerical_candidates=64,dense_ml_reml_cases=128,exhausted_budget_replays=64,
        retained_source_reviews_and_constants=exclusions,altered_numeric_cases_rejected=len(rejects),
        maximum_dense_objective_error=maximum_objective,maximum_dense_gradient_error=maximum_gradient,
        producer_status_counts=dict(primary_counts),independent_status_counts=dict(reader_counts)),artifact


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(exist_ok=False);assert not a.receipt.exists()
    grids,bindings,corruptions=source_grids(a.output)
    result,artifact=numerical(a.output)
    modules=['check_weighted_shared_entity_candidates','weighted_shared_entity_candidate','readback_weighted_shared_entity_candidate',
        'full_weighted_shared_entity_fit_sources','fit_shared_entity_likelihood','shared_entity_likelihood',
        'independent_shared_entity_likelihood','independent_shared_entity_likelihood_fast','independent_shared_entity_optimizer',
        'independent_positive_diagonal_basis_context','positive_diagonal_basis_context','full_weighted_covariance_qualification',
        'check_retained_shared_entity_candidates','check_shared_entity_likelihood']
    bindings.update({'scripts/'+m+'.py':sha('scripts/'+m+'.py') for m in modules});bindings[str(artifact)]=sha(artifact)
    verify(bindings)
    result.update(status='passed_complete_weighted_fit_source_and_actual_diagonal_candidate_contracts_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(),complete_source_grids=grids,
        total_candidate_identities=72000,synthetic_closure_alterations_rejected=corruptions,
        source_and_journal_fixtures_synthetic=True,source_hashes=bindings,scientific_eligibility=False,
        production_fitting_launched=False,full_fitting_export_workflow_implemented=False,
        scope='Three complete declared five-cohort synthetic source grids/72000candidate identities;64actual q4/q5/q6 uniform/nonuniform signed/unsigned ML/REML both-response fits with unchanged frozen optimizer and independent spectral/start/curvature/search checks. Numerical fixtures use arbitrary positive D; full source grids use closed reciprocal-reuse D. All actual failures/reviews retained. No biological pilot, full real source qualification, timing or scientific acceptance.')
    write(a.receipt,result);print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','scope']},indent=2),flush=True)


if __name__=='__main__':main()
