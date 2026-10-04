#!/usr/bin/env python3
"""Qualify complete parallel timing against serial arithmetic and private corruptions."""
import argparse
from collections import Counter
from copy import deepcopy
from datetime import datetime,timezone
import gzip
import json
from pathlib import Path
import shutil

from ancestral_chain_attempt import sha
from full_weighted_timing_contracts_v2 import semantic,same_numeric
from full_weighted_timing_parallel_v1 import run as parallel,SUMMARY
from prepare_weighted_timing_parallel_source_reference_v1 import run as serial,sources,SCHEMA
from readback_weighted_timing_parallel_source_reference_v1 import run as serial_reader
from full_weighted_shared_entity_fit_sources_parallel_v1 import cohorts,cases
from full_weighted_shared_entity_timing import groups
from prepare_baliphy_scalar_v6_preflight_v1 import project_sources
from reference_measurement_union_sources import bind,verify
from weighted_shared_entity_candidate import READY


def write(path,value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2);f.write('\n')


def rejected(action):
    try:action()
    except (AssertionError,ValueError,ArithmeticError,KeyError,FileExistsError,StopIteration):return
    raise AssertionError('Private corruption or original completed-stage restart accepted')


def private_clone(planpath,destination,root):
    """Rebind only private serialized fixture headers; never change qualified roots."""
    plan=json.loads(planpath.read_text());original=Path(plan['output'])
    shutil.copytree(original,root)
    (root/'readback.json').unlink();(root/'reader_completed.json').unlink()
    for path in (root/'checkpoints').glob('*.reader.json'):path.unlink()
    plan.update(output=str(root),scope='Private copied serialized negative timing fixture; producer not executed for this rebound plan. Actual native outputs originate from the separately qualified parent.')
    write(destination,plan);fit,source,bindings=sources(plan,destination)
    stage=dict(schema=SCHEMA,plan_sha256=sha(destination),fit_contract=source['fit_contract'])
    (root/'stage_plan.json').write_text(json.dumps(stage)+'\n')
    manifest=json.loads((root/'cohort_manifest.json').read_text())
    for index,entry in enumerate(manifest):
        cp=root/entry['receipt_path'];checkpoint=json.loads(cp.read_text());checkpoint['stage']=stage
        cp.write_text(json.dumps(checkpoint)+'\n');entry['receipt_sha256']=sha(cp)
        wp=root/'checkpoints'/(str(index).zfill(5)+'.producer.json')
        worker=json.loads(wp.read_text());worker['manifest']=entry;wp.write_text(json.dumps(worker)+'\n')
    (root/'cohort_manifest.json').write_text(json.dumps(manifest)+'\n')
    rp=root/'receipt.json';receipt=json.loads(rp.read_text())
    receipt.update(plan_sha256=sha(destination),source_hashes=bindings,scope=plan['scope'])
    receipt['artifacts']={name:sha(root/name) for name in receipt['artifacts']}
    rp.write_text(json.dumps(receipt)+'\n')
    return manifest


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--receipt',type=Path,required=True);a=p.parse_args();a.output.mkdir(exist_ok=False);assert not a.receipt.exists()
    gate=Path('metadata/weighted_parallel_fit_source_adapter_software_validation_20261004_v1.json')
    proven=json.loads(gate.read_text());verify(proven['source_hashes'])
    assert proven['status']=='passed_complete_parallel_numerical_to_fit_source_adapter_v1'
    pins={};project_sources(pins,[Path(__file__)]);bind(pins,gate)
    fixture=Path('data/software_audits/full-weighted-parallel-fit-source-adapter-20261004-v1')
    grids=[];total=group_count=0;positive_plan=None
    for name in ['unqualified','qualified-common-pair','qualified-pair-exception']:
        old=fixture/(name+'-parallel.fit.plan.json');fit=json.loads(old.read_text())
        fit.update(output=str(a.output/(name+'-unlaunched-fits')),pins={**fit['pins'],**pins})
        fp=a.output/(name+'.fit.plan.json');write(fp,fit)
        specs={}
        for version in ['serial','parallel']:
            spec=dict(fit_plan=str(fp),fit_plan_sha256=sha(fp),output=str(a.output/(name+'-'+version)),
                scaled_variance_points=[0.,1.],pins=dict(pins),resources=dict(cpus=2,memory_gib=32,
                    workers=2,worker_address_space_gib=12,reservation_capacity_gib=24,
                    parent_headroom_gib=8,minimum_free_disk_gib=128),
                scope='Complete synthetic five-cohort/24000candidate timing grid with real native numerics; source/journal closures synthetic, no biological pilot or fit.')
            tp=a.output/(name+'-'+version+'.timing.plan.json');write(tp,spec);specs[version]=tp
        sp=serial(specs['serial']);sr=serial_reader(specs['serial'],Path(json.loads(specs['serial'].read_text())['output'])/'readback.json')
        qp=parallel(specs['parallel']);qr=parallel(specs['parallel'],True)
        for field in SUMMARY:
            if field!='conditional_budget_weighted_seconds':assert sp[field]==sr[field]==qp[field]==qr[field],field
        assert qp['candidate_rows']==24000 and qp['unique_cohorts']==5
        roots={v:Path(json.loads(tp.read_text())['output']) for v,tp in specs.items()}
        sm=json.loads((roots['serial']/'cohort_manifest.json').read_text());qm=json.loads((roots['parallel']/'cohort_manifest.json').read_text())
        for se,qe in zip(sm,qm):
            assert se['cohort_id']==qe['cohort_id']
            with gzip.open(roots['serial']/se['census_path'],'rt') as sf,gzip.open(roots['parallel']/qe['census_path'],'rt') as qf:
                assert list(sf)==list(qf)
            left=json.loads((roots['serial']/se['probes_path']).read_text());right=json.loads((roots['parallel']/qe['probes_path']).read_text())
            same_numeric(semantic(left),semantic(right))
        # Independently derive the representative tuple from every candidate.
        spec=json.loads(specs['parallel'].read_text());fit,source,bindings=sources(spec,specs['parallel'])
        for cohort,rows,entries in cohorts(source,fit):
            census,selected,counts=groups(source,fit,cohort,entries);reference={};counter=Counter()
            for identity,x,y,audit,route,diagonal in cases(source,fit,cohort,entries):
                if identity['source_combined_disposition']!=READY:continue
                key=tuple(identity[k] for k in ['control_policy','loading_mode','tree','method','outcome'])
                rank=x.shape[1],float(audit['numerical_audit']['normalized_design_condition_number']),identity['candidate_id']
                reference[key]=max(rank,reference.get(key,rank));counter[key]+=1
            assert {r['key']:r['rank'] for r in selected.values()}==reference
            assert {r['key']:counts[gid] for gid,r in selected.items()}==counter
            assert len(census)==4800
        workers={}
        for cp in (roots['parallel']/'checkpoints').glob('*.producer.json'):
            record=json.loads(cp.read_text())
            assert record['cached_numeric_inputs_preserved'] and record['submitted_numeric_inputs_preserved']
            assert record['submitted_task_array_bindings_checked']>0 and record['worker_address_space_limit_bytes']==12*2**30
            workers[record['worker']['pid']]=record['worker']['created']
        assert len(workers)==2
        for r in [sr,qr]:
            verify(r['source_hashes'])
            for path,digest in r['source_hashes'].items():bind(pins,path,digest)
        for root in roots.values():
            for path in root.rglob('*'):
                if path.is_file():bind(pins,path)
        rejected(lambda:parallel(specs['parallel']));rejected(lambda:parallel(specs['parallel'],True,a.output/(name+'-restart-reader.json')))
        total+=qp['candidate_rows'];group_count+=qp['timing_groups']
        grids.append(dict(name=name,candidates=qp['candidate_rows'],timing_groups=qp['timing_groups'],native_producer_workers=len(workers)))
        if name=='qualified-common-pair':positive_plan=specs['parallel']
    assert total==72000 and group_count==320
    # A positive rebound copy must pass before negative copies are meaningful.
    control_plan=a.output/'private-control.plan.json';control_root=a.output/'private-control'
    private_clone(positive_plan,control_plan,control_root);parallel(control_plan,True)
    rejected_cases=[]
    for case in ['omit_candidate','duplicate_candidate','control_identity','selected_candidate','eligible_count','missing_probe',
        'false_success','invented_review','changed_diagonal','changed_planning','missing_artifact','extra_artifact',
        'source_hash','stage_contract','cache_guard_false','task_guard_false','worker_limit']:
        tp=a.output/('private-'+case+'.plan.json');root=a.output/('private-'+case)
        manifest=private_clone(positive_plan,tp,root)
        part=next(entry for entry in manifest if json.loads((root/entry['probes_path']).read_text()))
        fp=root/part['census_path'];pp=root/part['probes_path'];cp=root/part['receipt_path']
        wp=next(path for path in (root/'checkpoints').glob('*.producer.json') if json.loads(path.read_text())['cohort_id']==part['cohort_id'])
        rp=root/'receipt.json';receipt=json.loads(rp.read_text())
        census=[json.loads(line) for line in gzip.decompress(fp.read_bytes()).decode().splitlines()]
        probes=json.loads(pp.read_text());checkpoint=json.loads(cp.read_text())
        if case=='omit_candidate':census.pop()
        elif case=='duplicate_candidate':census[-1]=deepcopy(census[0])
        elif case=='control_identity':census[0]['identity_sha256']='foreign'
        elif case=='selected_candidate':probes[0]['representative']['candidate_id']='foreign'
        elif case=='eligible_count':probes[0]['eligible_candidates']+=1
        elif case=='missing_probe':probes.pop()
        elif case=='false_success':probes[0]['points'][0]['maximum_coordinate_gradient_error']=1.
        elif case=='invented_review':
            probes[0]['status']='timing_group_requires_review'
            probes[0]['points'][0]=dict(scaled_variance=0.,status='timing_probe_precision_requires_review',error_type='ArithmeticError',error_message='Invented private review',elapsed_seconds=.1)
        elif case=='changed_diagonal':probes[0]['representative']['diagonal_sha256']='foreign'
        elif case=='changed_planning':
            path=root/'conditional_planning.json';v=json.loads(path.read_text());v['conditional_budget_weighted_seconds']+=1;path.write_text(json.dumps(v)+'\n')
        elif case=='missing_artifact':receipt['artifacts'].pop(part['probes_path'])
        elif case=='extra_artifact':
            path=root/'extra.json';path.write_text('{}\n');receipt['artifacts'][path.name]=sha(path)
        elif case=='source_hash':receipt['source_hashes'][str(tp)]='0'*64
        elif case=='stage_contract':checkpoint['stage']['fit_contract']='foreign'
        else:
            value=json.loads(wp.read_text())
            if case=='cache_guard_false':value['cached_numeric_inputs_preserved']=False
            elif case=='task_guard_false':value['submitted_numeric_inputs_preserved']=False
            elif case=='worker_limit':value['worker_address_space_limit_bytes']+=1
            wp.write_text(json.dumps(value)+'\n')
        fp.write_bytes(gzip.compress(''.join(json.dumps(row)+'\n' for row in census).encode(),mtime=0));pp.write_text(json.dumps(probes)+'\n')
        checkpoint.update(census_sha256=sha(fp),probes_sha256=sha(pp));cp.write_text(json.dumps(checkpoint)+'\n')
        part.update(census_sha256=sha(fp),probes_sha256=sha(pp),receipt_sha256=sha(cp))
        (root/'cohort_manifest.json').write_text(json.dumps(manifest)+'\n')
        value=json.loads(wp.read_text());value['manifest']=part;wp.write_text(json.dumps(value)+'\n')
        for path in [fp,pp,cp,wp,root/'cohort_manifest.json',root/'conditional_planning.json']:
            name=str(path.relative_to(root))
            if name in receipt['artifacts']:receipt['artifacts'][name]=sha(path)
        rp.write_text(json.dumps(receipt)+'\n')
        rejected(lambda:parallel(tp,True));assert not (root/'readback.json').exists()
        for path in root.rglob('*'):
            if path.is_file():bind(pins,path)
        bind(pins,tp);rejected_cases.append(case)
        print('parallel_timing_private_corruption_rejected',case,flush=True)
    verify(pins)
    result=dict(status='passed_complete_parallel_four_control_timing_and_numeric_readback_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(),complete_synthetic_grids=grids,
        total_candidate_rows=total,full_grid_selected_groups=group_count,
        serial_parallel_numeric_probe_groups_compared=group_count,independent_numeric_groups_replayed=group_count,
        actual_variance_points_per_group=2,private_rehashed_corruptions_rejected=rejected_cases,
        full_grid_selections_independently_recomputed=True,positive_private_rebound_control_passed=True,
        completed_producer_reader_restarts_refused=6,cached_and_submitted_inputs_preserved_all_cohorts=True,
        source_and_journal_fixtures_synthetic=True,full_grid_native_probes_mocked=False,
        source_hashes=pins,fits_computed=0,scientific_eligibility=False,
        scope='Three complete five-cohort/72000candidate grids with real serial/parallel native '
              'timing and independent numeric replay; every selected group retained. Hardware '
              'durations differ and are not numerical equality claims. Positive private copied '
              'serialized fixture passes before independent private corruptions are tested; '
              'those rebound fixtures do not claim new producer executions. Source/journal '
              'containers remain synthetic. No production timing, fit or biological acceptance.')
    with a.receipt.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','scope']},indent=2))


if __name__=='__main__':main()
