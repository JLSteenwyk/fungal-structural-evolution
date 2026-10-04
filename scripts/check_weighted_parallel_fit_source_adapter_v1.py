#!/usr/bin/env python3
"""Qualify full synthetic candidate inputs from closed parallel numerical artifacts."""
import argparse
from copy import deepcopy
from datetime import datetime,timezone
import json
from pathlib import Path

import numpy as np

from ancestral_chain_attempt import sha
from check_full_inverse_reuse_weights import closure
from full_weighted_covariance_qualification import SUMMARY
from full_weighted_shared_entity_fit_sources import load as serial_load,cohorts as serial_cohorts,cases as serial_cases
from full_weighted_shared_entity_fit_sources_parallel_v1 import load,cohorts,cases
from prepare_baliphy_scalar_v6_preflight_v1 import project_sources
from reference_measurement_union_sources import bind,verify
from weighted_parallel_checkpoint_contracts_v1 import checkpoint_pair


def write(path,value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2);f.write('\n')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--receipt',type=Path,required=True);a=p.parse_args();assert not a.receipt.exists();a.output.mkdir(exist_ok=False)
    pp=Path('metadata/weighted_parallel_fixtures_software_preparation_20261004_v4.json')
    prepared=json.loads(pp.read_text());verify(prepared['source_hashes'])
    pins={};project_sources(pins,[Path(__file__)]);bind(pins,pp)
    originalfit=Path('metadata/full_weighted_shared_entity_fit_draft_plan_20261004_v1.json')
    config=json.loads(originalfit.read_text());verify(config['pins']);bind(pins,originalfit)
    compared=0;grids=[];first=None
    for grid in prepared['grids']:
        plans={}
        for version in ['serial','parallel']:
            qp=Path(grid[version]);q=json.loads(qp.read_text());root=Path(q['output'])
            rp=root/'receipt.json';rb=root/'readback.json';producer=json.loads(rp.read_text());reader=json.loads(rb.read_text())
            paths=set(map(Path,reader['source_hashes']))|{qp,rp,rb}
            paths.update(root/name for name in producer['artifacts'])
            paths.update(root/name for name in reader.get('reader_checkpoint_artifacts',{}))
            cp=a.output/(grid['name']+'-'+version+'.closed.json')
            closure(cp,'complete_verified_full_four_control_covariance_numerical_qualification_v1',sorted(paths),
                {**{k:reader[k] for k in SUMMARY},'producer_receipt':str(rp),'producer_receipt_sha256':sha(rp),
                 'independent_readback':str(rb),'independent_readback_sha256':sha(rb)})
            fit=dict(config,qualification_plan=str(qp),qualification_completion=str(cp),trees=q['trees'],
                output=str(a.output/(grid['name']+'-'+version+'-unlaunched-fits')),
                expected=dict(logical_cases=24,cohorts=5,designs=150,fit_inputs=300,settings=600,
                    candidate_rows=24000,setting_fit_links=48000),pins={**config['pins'],**pins})
            fp=a.output/(grid['name']+'-'+version+'.fit.plan.json');write(fp,fit);plans[version]=(fp,fit)
        ps,pfit=plans['parallel'];ss,sfit=plans['serial']
        source,bindings=load(pfit,ps);reference,rbindings=serial_load(sfit,ss)
        for mapping in [bindings,rbindings]:
            for path,digest in mapping.items():bind(pins,path,digest)
        count=0
        for (pc,pr,pe),(sc,sr,se) in zip(cohorts(source,pfit),serial_cohorts(reference,sfit)):
            assert pc==sc and np.array_equal(pr,sr)
            for left,right in zip(cases(source,pfit,pc,pe),serial_cases(reference,sfit,sc,se)):
                pi,px,py,pa,pb,pd=left;si,sx,sy,sa,sb,sd=right
                excluded={'candidate_id','source_contract','covariance_audit_id','covariance_audit_sha256'}
                assert {k:v for k,v in pi.items() if k not in excluded}=={k:v for k,v in si.items() if k not in excluded}
                assert np.array_equal(px,sx) and np.array_equal(py,sy) and np.array_equal(pd,sd)
                assert pb==sb
                assert {k:v for k,v in pa.items() if k not in ['audit_id','source_contract']}=={k:v for k,v in sa.items() if k not in ['audit_id','source_contract']}
                count+=1
        assert count==24000;compared+=count;grids.append(dict(name=grid['name'],candidate_inputs=count))
        if first is None:
            root=source['numerical_root'];entry=source['numerical_manifest'][0]
            first=(json.loads((root/'checkpoints/00000.producer.json').read_text()),json.loads((root/'checkpoints/00000.reader.json').read_text()),entry)
    assert compared==72000
    original,independent,entry=first
    checkpoint_pair(original,independent,entry,0)
    corruptions=['cohort_index','cohort_id','entry','audits','links','designs','counts','link_counts',
        'policy_counts','cache_changed','cache_unchecked','worker_as','scientific_eligibility','worker_pid','worker_created','worker_command']
    rejected=[]
    for case in corruptions:
        value=deepcopy(independent)
        if case in ['cohort_index','audits','links','designs']:value[case]+=1
        elif case=='cohort_id':value[case]='foreign'
        elif case=='entry':value[case]['audit_sha256']='foreign'
        elif case in ['counts','link_counts','policy_counts']:value[case]['foreign']=1
        elif case=='cache_changed':value['cached_numeric_inputs_preserved']=False
        elif case=='cache_unchecked':value['cached_input_array_bindings_checked']=0
        elif case=='worker_as':value['worker_address_space_limit_bytes']+=1
        elif case=='scientific_eligibility':value[case]=True
        elif case=='worker_pid':value['worker']['pid']=0
        elif case=='worker_created':value['worker']['created']=0
        elif case=='worker_command':value['worker']['cmdline']=[]
        path=a.output/(case+'.private-checkpoint.json');write(path,value);loaded=json.loads(path.read_text())
        try:checkpoint_pair(original,loaded,entry,0)
        except AssertionError:rejected.append(case)
        else:raise AssertionError('Private malformed checkpoint accepted: '+case)
        bind(pins,path)
    assert len(rejected)==16;verify(pins)
    result=dict(status='passed_complete_parallel_numerical_to_fit_source_adapter_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(),complete_grids=grids,
        serial_parallel_candidate_inputs_compared=compared,private_serialized_checkpoint_corruptions_rejected=rejected,
        matrices_responses_diagonals_and_audits_exactly_compared=True,
        complete_checkpoint_artifact_sets_checked=True,source_hashes=pins,
        synthetic_numerical_closure_fixtures=True,working_model_fits_computed=0,scientific_eligibility=False,
        scope='All candidate inputs across three complete synthetic grids compared to the '
              'frozen serial source adapter. Closed parallel producer and independent reader '
              'checkpoint sets are additionally verified; sixteen private serialized checkpoint '
              'alterations rejected. Synthetic closure containers are explicit. No full actual '
              'numerical closure, timing run, fit or biological acceptance.')
    with a.receipt.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
