#!/usr/bin/env python3
"""Complete retained source/grid/checkpoint/link and positive numeric contracts."""
import argparse
from collections import Counter
from copy import deepcopy
from datetime import datetime, timezone
import csv
import gzip
import itertools
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from ancestral_chain_attempt import sha
from check_full_shared_entity_fits_v2 import mutate as original_mutate, settings
from check_retained_shared_entity_candidates import inputs
from full_expanded_model_design_sources import digest, array_digest, fit_id
from full_retained_shared_entity_fit_fixtures import setup
from full_retained_shared_entity_fit_sources import load, identity, operators_for
from prepare_full_retained_shared_entity_fits import run, candidate
from readback_full_retained_shared_entity_fits import run as readback, numeric
from reference_measurement_union_sources import verify


def rejected(action):
    try:action()
    except (AssertionError,ValueError,ArithmeticError,KeyError,FileExistsError,FileNotFoundError,StopIteration):return
    raise AssertionError('Altered retained full fitting accepted')


def mutate(root,name):
    if name in ['wrong_original_audit_sha','wrong_exact_certificate_sha','wrong_retained_basis','wrong_method']:
        output=root/'output';receipt=json.loads((output/'receipt.json').read_text())
        manifest=json.loads((output/'cohort_manifest.json').read_text());part=manifest[0];fp=output/part['path']
        with gzip.open(fp,'rt') as f:rows=[json.loads(line) for line in f]
        if name=='wrong_original_audit_sha':rows[0]['original_covariance_audit_sha256']='foreign'
        elif name=='wrong_exact_certificate_sha':rows[0]['exact_covariance_certificate_sha256']='foreign'
        elif name=='wrong_retained_basis':rows[0]['retained_kernel_names'].remove('background_node')
        else:rows[0]['method']='foreign'
        with gzip.open(fp,'wt') as f:f.write(''.join(json.dumps(r)+'\n' for r in rows))
        part['sha256']=sha(fp)
        cp=output/part['receipt_path'];saved=json.loads(cp.read_text());saved['candidate_sha256']=sha(fp)
        cp.write_text(json.dumps(saved,indent=2)+'\n');part['receipt_sha256']=sha(cp)
        (output/'cohort_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        for name in [part['path'],part['receipt_path'],'cohort_manifest.json']:receipt['artifacts'][name]=sha(output/name)
        (output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    else:original_mutate(root,name)


def positive_candidates(folder):
    folder.mkdir();plan=settings();outcomes=['rmsd_delta','native_tm_dissimilarity_delta']
    exports=[];checks=Counter()
    for exception,mode in itertools.product([False,True],['signed','unsigned']):
        source,rows,x,response,operators,audit,old,cert=inputs(exception,mode)
        source.update(fit_contract='explicit-synthetic-complete-candidate-export-contract',retained_contract=audit['retained_source_contract'])
        n=len(rows);target=sparse.eye(n,format='csr');background=operators['background_node']
        index=np.arange(n)
        if exception:index[1]=index[0]
        pair_target=sparse.csr_matrix((np.ones(n),(np.arange(n),index)),shape=(n,n))
        pair=sparse.hstack([pair_target,background],format='csr');family=operators['family_intercept']
        gene=.5*sparse.hstack([target,background,target,background],format='csr')
        model=.5*sparse.hstack([pair,pair],format='csr')
        source['operators']={('contrast','family_intercept'):family}
        for m in ['signed','unsigned']:
            for name,z in dict(target_node=target,background_node=background,model_pair=pair,gene=gene,model=model,
                               family=sparse.csr_matrix(family.shape) if m=='signed' else 2*family).items():
                source['operators'][m,name]=z
        selected=operators_for(source,rows,audit)
        assert list(selected)==list(operators) and all((selected[k]-operators[k]).nnz==0 for k in selected)
        cohort=dict(cohort_id=cert['cohort_id'],ordered_case_ids_sha256=cert['ordered_case_ids_sha256'])
        design=dict(design_id=audit['design_id'],raw_design_sha256=audit['raw_design_sha256'],
            active_column_indices=[0,1,2],exactly_zero_column_indices=[],predictor_columns=['intercept','sequence1','sequence2'])
        for outcome,method in itertools.product(outcomes,['ml','reml']):
            y=response.copy() if outcome==outcomes[0] else response+1.
            response_sha=array_digest(y,'<f8')
            fit=dict(fit_input_id=fit_id(design['design_id'],outcome,response_sha),outcome=outcome,
                disposition='ready_for_working_covariance_fit',records=n,response_sha256=response_sha)
            expected=identity(source,cohort,design,fit,audit,old,cert,mode,audit['tree'],method)
            value=candidate(source,plan,rows,expected,x,y,selected,audit,old,cert)
            # Exact JSON roundtrip is the full producer's serialization boundary.
            exported=json.loads(json.dumps(value,sort_keys=True,allow_nan=False))
            checked=numeric(source,plan,rows,expected,x,y,exported,selected,audit,old,cert)
            assert exported['numerical_attempted'] is True and exported['fit'] is not None
            assert exported['fit']['parameter_names']==audit['retained_kernel_names'][1:]
            checks[checked['disposition']]+=1
            exports.append(dict(candidate=exported,independent_audit=checked,original_audit=old,retained_audit=audit,certificate=cert))
            for field in ['original_covariance_audit_sha256','exact_covariance_certificate_sha256','retained_kernel_names']:
                changed=deepcopy(exported);changed[field]='foreign'
                rejected(lambda:numeric(source,plan,rows,expected,x,y,changed,selected,audit,old,cert))
    artifact=folder/'serialized_positive_retained_candidates.json'
    artifact.write_text(json.dumps(exports,indent=2,allow_nan=False)+'\n')
    assert len(exports)==16
    return dict(cases=16,independent_dispositions=dict(checks),altered_identity_cases=48,artifact=str(artifact),artifact_sha256=sha(artifact))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--receipt',type=Path,required=True);a=p.parse_args();root=a.output.resolve();root.mkdir(exist_ok=False)
    assert not a.receipt.exists()
    positive=positive_candidates(root/'positive')
    fixture=root/'fixture';fixture.mkdir();pp=setup(fixture)
    produced=run(pp);read=readback(pp,fixture/'readback.json')
    assert (produced['unique_cohorts'],produced['model_setting_rows'],produced['candidate_rows'],produced['setting_fit_links'])==(5,600,6000,12000)
    assert produced['candidate_status_counts']==read['candidate_status_counts']
    assert produced['setting_status_counts']==read['setting_status_counts']
    assert read['independent_candidate_status_counts']=={'source_review_or_exclusion_retained':6000}
    rejected(lambda:run(pp));rejected(lambda:readback(pp,fixture/'duplicate-readback.json'))
    output=fixture/'output';snapshot={str(p.relative_to(output)):p.read_bytes() for p in output.rglob('*') if p.is_file()}
    changes=['omit_last_tree','omit_unsigned','duplicate_candidate','wrong_fit_input','wrong_response_hash',
        'wrong_covariance_audit','wrong_active_columns','false_scientific_acceptance','promote_source_review',
        'link_omit_reml','link_wrong_candidate','link_wrong_status',
        'wrong_original_audit_sha','wrong_exact_certificate_sha','wrong_retained_basis','wrong_method']
    rejected_cases=[]
    for change in changes:
        mutate(fixture,change)
        (output/'readback_stage.json').unlink();(output/'independent_readback_completed.json').unlink()
        rejected(lambda:readback(pp,fixture/'bad-readback.json'));rejected_cases.append(change)
        for p in output.rglob('*'):
            if p.is_file() and str(p.relative_to(output)) not in snapshot:p.unlink()
        for name,raw in snapshot.items():(output/name).write_bytes(raw)
    assert all((output/name).read_bytes()==raw for name,raw in snapshot.items())
    source,bindings=load(json.loads(pp.read_text()),pp);verify(bindings)
    source_cases=[];plan=json.loads(pp.read_text());reduced=Path(plan['retained_plan']);rp=json.loads(reduced.read_text())
    completion=Path(plan['retained_completion']);compact=json.loads(completion.read_text());archive=Path(compact['full_hash_archive'])
    reader=Path(rp['output'])/'readback.json';saved={p:p.read_bytes() for p in [reader,completion,archive]}
    def altered_reader(name,change):
        value=json.loads(saved[reader]);change(value);reader.write_text(json.dumps(value,indent=2)+'\n')
        proof=json.loads(saved[archive]);proof['source_hashes'][str(reader)]=sha(reader);archive.write_text(json.dumps(proof,indent=2)+'\n')
        c=dict(compact,full_hash_archive_sha256=sha(archive),independent_readback_sha256=sha(reader));completion.write_text(json.dumps(c,indent=2)+'\n')
        rejected(lambda:load(plan,pp));source_cases.append(name)
        for p,raw in saved.items():p.write_bytes(raw)
    for name,change in [
        ('unqualified_retained_reader',lambda r:r.__setitem__('status','pending')),
        ('foreign_retained_contract',lambda r:r.__setitem__('source_contract','foreign')),
        ('foreign_retained_plan',lambda r:r.__setitem__('plan_sha256','foreign')),
        ('missing_retained_tree',lambda r:r.__setitem__('trees',r['trees'][:-1])),
        ('promoted_retained_science',lambda r:r.__setitem__('scientific_eligibility',True)),
        ('foreign_retained_producer',lambda r:r.__setitem__('producer_receipt_sha256','foreign'))]:altered_reader(name,change)
    assert all(p.read_bytes()==raw for p,raw in saved.items())
    verify(produced['source_hashes']);verify(read['source_hashes'])
    interrupted=root/'interrupted';interrupted.mkdir();pp2=setup(interrupted)
    try:run(pp2,stop_after_cohorts=2)
    except InterruptedError:pass
    else:raise AssertionError('Synthetic checkpoint interruption not exercised')
    parts=list((interrupted/'output/cohorts').glob('*.jsonl.gz'));assert len(parts)==2
    hashes={str(p):sha(p) for p in parts};run(pp2);readback(pp2,interrupted/'readback.json')
    assert all(sha(p)==h for p,h in hashes.items())
    bindings={str(p):sha(p) for p in [Path(__file__),*[Path('scripts')/n for n in
        ['full_retained_shared_entity_fit_sources.py','prepare_full_retained_shared_entity_fits.py',
         'readback_full_retained_shared_entity_fits.py','full_retained_shared_entity_fit_fixtures.py',
         'retained_shared_entity_candidate.py','readback_retained_shared_entity_candidate.py',
         'full_shared_entity_fit_sources.py','full_reduced_covariance_sources_v3.py']]]}
    for p in root.rglob('*'):
        if p.is_file():bindings[str(p)]=sha(p)
    verify(bindings)
    result=dict(status='passed_full_retained_shared_entity_source_grid_checkpoint_spectral_contracts_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(),synthetic_declared_cohorts=5,synthetic_model_settings=600,
        synthetic_candidate_rows=6000,synthetic_setting_fit_links=12000,all_declared_cohorts_settings_trees_modes_outcomes_methods_retained=True,
        positive_serialized_candidate_cases=positive['cases'],positive_independent_dispositions=positive['independent_dispositions'],
        positive_rehashed_identity_cases_rejected=positive['altered_identity_cases'],
        rehashed_grid_cases_rejected=rejected_cases,rehashed_source_cases_rejected=source_cases,
        original_retained_audit_and_exact_certificate_bound=True,inherited_numerical_envelopes_never_replaced=True,
        full_q4_q5_numerical_parameter_order_preserved=True,completed_producer_and_reader_restarts_refused=True,
        synthetic_interrupted_closed_chunks_reused_without_rewrite=True,complete_positive_source_graph_rehashed=True,
        source_and_journal_fixtures_synthetic=True,full_production_source_grid_loaded=False,production_fitting_launched=False,
        full_scope_timing_completed=False,scientific_eligibility=False,source_hashes=bindings,
        scope='Complete declared five-positive-membership-cohort source fixture:600settings6000candidates12000links, '
              'all original non-ready/constant/rank/review states and five trees/both modes/outcomes/MLREML; sixteen '
              'actual numerical q4/q5 serialized candidates separately exercise both outcomes/modes/methods. '
              'Parent/source journals synthetic,48altered numerical identities16altered output grids6source closures '
              'rejected; exact restoration and positive hashes replayed. Production still requires all4340cohorts '
              'and5208000candidates12441600links, real closed sources and full-scope timing. No biological pilot or acceptance.')
    a.receipt.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
