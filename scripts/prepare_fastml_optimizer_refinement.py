#!/usr/bin/env python3
"""Prepare all original indel inputs for stricter, fixed-tree multistart fitting."""
import argparse
import copy
import json
from pathlib import Path
from ancestral_chain_attempt import sha,write_json


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source-plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    source=json.loads(args.source_plan.read_text())
    for path,h in source['pins'].items():assert sha(path)==h
    groups={}
    for item in source['jobs']:groups.setdefault(item['job']['job_id'],{})[item['variant']]=item
    assert len(groups)==156
    args.output.mkdir(parents=True,exist_ok=False);jobs=[];bindings={};trees=0
    for job_id,pair in sorted(groups.items()):
        assert set(pair)=={'precision-only','precision-cache-refresh'}
        item=pair['precision-cache-refresh'];job=item['job'];assert job==pair['precision-only']['job']
        if not job['character_count']:
            for label in ['precision_parameters','refreshed_parameters','low_ratio','balanced','high_ratio']:
                jobs.append(dict(id=job_id+'__'+label,start=label,job=job,config=None))
            continue
        prior={};tree_paths=[]
        for variant,old in pair.items():
            rp=Path(source['output'])/old['id']/'readback.json';readback=json.loads(rp.read_text());bindings[str(rp)]=sha(rp)
            attempt=Path(readback['attempt_receipt']);assert sha(attempt)==readback['attempt_receipt_sha256']
            ar=json.loads(attempt.read_text());assert ar['exit_code']==0
            tp=attempt.parent/'RESULTS/TheTree.INodes.ph';assert sha(tp)==ar['artifacts']['RESULTS/TheTree.INodes.ph'];tree_paths.append(tp)
            rr=Path('results/ancestral/full-fastml-roundoff-readback-20260928-v1')/(old['id']+'.json')
            prior[variant]=json.loads(rr.read_text())['fitted_parameters'];bindings[str(rr)]=sha(rr)
        assert sha(tree_paths[0])==sha(tree_paths[1]),'Paired fitted tree changed'
        tree=args.output/(job_id+'.nwk');tree.write_bytes(tree_paths[0].read_bytes());trees+=1
        base={}
        parameter_file=Path(item['config']['command'][-1]);bindings[str(parameter_file)]=sha(parameter_file)
        for line in parameter_file.read_text().splitlines():
            if not line.strip() or line.startswith('#'):continue
            key,value=line.split(None,1);assert key not in base;base[key]=value
        base.update(_optimizationLevel='mid',_maxNumOfIterations='100',_maxNumOfIterationsModel='100',
                    _epsilonOptimizationModel='0.000001',_epsilonOptimizationIterationCycle='0.000001',
                    _correctOptimizationEpsilon='0',_performOptimizationsManyStarts='0',
                    _isInitGainLossByEmpiricalFreq='0',_isMultipleAllBranchesByFactorAtStart='0',_treeFile=str(tree.resolve()))
        starts={'precision_parameters':prior['precision-only'],'refreshed_parameters':prior['precision-cache-refresh']}
        for label,alpha,ratio in [('low_ratio',.2,.1),('balanced',1.,1.),('high_ratio',20.,10.)]:
            starts[label]=dict(alpha=alpha,gain=(1+ratio)/2,loss=(1+ratio)/(2*ratio))
        for label,parameters in starts.items():
            config=copy.deepcopy(item['config']);identifier=job_id+'__'+label
            options={**base,'_userAlphaRate':repr(parameters['alpha']),'_userGain':repr(parameters['gain']),'_userLoss':repr(parameters['loss'])}
            param=args.output/(identifier+'.txt');param.write_text(''.join(k+' '+v+'\n' for k,v in options.items()))
            old_param=config['command'][-1];config['command'][-1]=str(param.resolve());config['pins'].pop(old_param)
            config['pins'].update({str(param.resolve()):sha(param),str(tree.resolve()):sha(tree)})
            config['model_input_identity']=identifier;config['timeout_seconds']=43200
            jobs.append(dict(id=identifier,start=label,job=job,config=config,initial_parameters=parameters,
                expected_effective_options={k:options[k] for k in ['_optimizationLevel','_maxNumOfIterations','_maxNumOfIterationsModel','_epsilonOptimizationModel','_epsilonOptimizationIterationCycle','_correctOptimizationEpsilon','_isInitGainLossByEmpiricalFreq','_isMultipleAllBranchesByFactorAtStart']}))
    assert len(jobs)==780 and sum(j['config'] is not None for j in jobs)==765 and trees==153
    write_json(args.output/'jobs.json',jobs)
    result=dict(status='full_multistart_refinement_inputs_prepared_not_launched',inputs=156,starts_per_input=5,nonempty_fits=765,empty_dispositions=15,
        identical_source_tree_pairs=trees,pins={str(args.source_plan):sha(args.source_plan),str(Path(__file__)):sha(__file__),**bindings},
        resources=dict(proposed_workers=8,cpus=8,memory_gib=48,swap_gib=0,per_fit_address_space_gib=4,storage_allowance_gib=128,planning_hours=[12,336],per_fit_timeout_hours=12,gpu=False,paid_cost=0,basis='Original153 fits used approximately5.2 CPU-hours at one model iteration. Five starts and up to100 outer/model iterations can substantially increase work; broad planning range is not a reliable ETA or convergence bound.'),
        artifacts={p.name:sha(p) for p in args.output.iterdir()},
        scope='All156 input designs preserved, including empty inputs. Uses identical serialized fitted trees with initial MP rescaling disabled, normalized stationary rate, two prior parameter starts plus three fixed starts. mid prevents low-level overrides; future runner must verify effective options from native output before accepting fits. No refined fits yet run; numerical and model-adequacy checks remain.')
    write_json(args.output/'receipt.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['pins','artifacts']},indent=2))


if __name__=='__main__':main()
