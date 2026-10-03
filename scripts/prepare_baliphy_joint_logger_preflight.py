#!/usr/bin/env python3
"""Prepare every future joint-logger startup role; no posterior sampling."""
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil

import psutil

from ancestral_chain_attempt import sha, write_json
from baliphy_joint_node_logger_v3 import restore
from run_baliphy_reference_preflight import verify


def build(roles, original_jobs, binary, prlimit, api_paths):
    assert len(roles) == len(original_jobs) == 1620
    originals = {j['chain']['chain_id']: j for j in original_jobs}
    assert len(originals) == 1620
    assert {r['source_chain_id'] for r in roles} == set(originals)
    assert len({r['chain']['chain_id'] for r in roles}) == 1620
    seeds = {r['chain']['seed'] for r in roles}
    forbidden = {j['chain']['seed'] for j in original_jobs} | {j['source_seed'] for j in original_jobs} | {20267001,20267002,20267003}
    assert len(seeds) == 1620 and seeds.isdisjoint(forbidden)
    assert all(type(seed) is int and 1 <= seed <= 2**31-1 for seed in seeds)
    jobs = []; groups = defaultdict(list)
    for role in roles:
        future = role['chain']; original = originals[role['source_chain_id']]['chain']
        assert role['native_execution_launched'] is role['scientific_eligibility'] is role['posterior_qualified'] is False
        assert role['current_short_sampler_seed'] == original['seed']
        shared = ['effective_input_group','prior_label','chain','original_configuration_ids','family',
                  'proteins','alignment','alignment_sha256','tree','tree_sha256']
        assert all(future[k] == original[k] for k in shared)
        for c in [original,future]:
            assert sha(c['program']) == c['program_sha256']
        assert restore(Path(future['program']).read_text()) == Path(original['program']).read_text()
        for name in ['alignment','tree']: assert sha(future[name]) == future[name+'_sha256']
        chain = dict(future, seed=original['seed'], program=original['program'],program_sha256=original['program_sha256'])
        paths = [binary,prlimit,Path(future['program']),Path(future['alignment']),Path(future['tree']),*api_paths]
        config = dict(command=[str(prlimit),'--as='+str(12*2**30),'--cpu=600','--fsize='+str(256*2**20),
            '--',str(binary),'--seed',str(future['seed']),'run',future['program'],'--test','--log-format','json'],
            timeout_seconds=900,pins={str(p.resolve()):sha(p) for p in paths})
        jobs.append(dict(chain=chain,fresh_seed=future['seed'],source_chain_id=role['source_chain_id'],
                         program=future['program'],config=config))
        groups[future['effective_input_group']+'-'+future['prior_label']].append(future)
    assert len(groups)==405
    assert all(len(cs)==4 and {c['chain'] for c in cs}=={1,2,3,4} for cs in groups.values())
    assert Counter(j['chain']['prior_label'] for j in jobs)=={'broad':540,'centered':540,'package':540}
    assert len({j['chain']['effective_input_group'] for j in jobs})==135
    assert len({a for j in jobs for a in j['chain']['original_configuration_ids']})==324
    assert all(j['config']['command'][-3:]==['--test','--log-format','json'] and
               '--iterations' not in j['config']['command'] for j in jobs)
    return jobs


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--validation',type=Path,required=True)
    p.add_argument('--inputs',type=Path,required=True);p.add_argument('--plan',type=Path,required=True)
    a=p.parse_args();gate=json.loads(a.validation.read_text())
    assert gate['status']=='passed_full_joint_logger_preflight_software_contracts'
    assert gate['full_roles']==1620 and gate['full_quartets']==405 and gate['native_startups']==3
    assert gate['full_producer_reader_serialization_checked'] and gate['artificial_failures_retained']==2
    verify({'pins':gate['source_hashes']})
    future_path=Path('metadata/baliphy_joint_node_logger_future_models_20261003_v4.json')
    future=json.loads(future_path.read_text());verify(future)
    native_path=Path(future['source_plan']);native=json.loads(native_path.read_text());verify(native)
    roles=json.loads(Path(future['future_roles']).read_text());original_jobs=json.loads(Path(native['jobs']).read_text())
    binary=Path(original_jobs[0]['config']['command'][5]);prlimit=Path('/usr/bin/prlimit')
    api=Path('data/software_audits/baliphy-4.3-20260927/install/bali-phy-4.3/lib/bali-phy/haskell')
    apis=[api/name for name in ['Bio/Alignment.hs','Bio/Alphabet.hs','Graph.hs','Tree.hs','Data/Text.hs',
        'Data/OldList.hs','Probability/Logger.hs','Data/JSON/Encoding.hs','Data/JSON/Types/ToJSON.hs',
        'Probability/Distribution/PhyloAlignment.hs','Probability/Random.hs','MCMC.hs']]
    root=a.inputs.resolve();root.mkdir(exist_ok=False)
    jobs=build(roles,original_jobs,binary,prlimit,apis);job_path=root/'jobs.json';write_json(job_path,jobs)
    # Full closed previous startup census supplies a planning comparison only.
    closed_path=Path('metadata/baliphy_reference_startup_footer_completed_20261003.json')
    closed=json.loads(closed_path.read_text());assert closed['validated_startups']==1620
    archive=Path(closed['full_hash_archive']);assert sha(archive)==closed['full_hash_archive_sha256']
    proof=json.loads(archive.read_text());history=Path(closed['producer_receipt']).parent/'dispositions.json'
    assert sha(history)==proof['source_hashes'][str(history)]
    rows=json.loads(history.read_text());assert len(rows)==1620
    seconds=sum(r['elapsed_seconds'] for r in rows)
    own=['prepare_baliphy_joint_logger_preflight','run_baliphy_joint_logger_preflight',
         'readback_baliphy_joint_logger_preflight','check_baliphy_joint_logger_preflight',
         'launch_baliphy_joint_logger_preflight','record_baliphy_joint_logger_preflight_checkpoint',
         'baliphy_joint_node_logger_v3','baliphy_full_node_logger','baliphy_reference_startup_readback_v2',
         'baliphy_reference_initialization','ancestral_chain_attempt']
    paths=[Path('scripts',name+'.py') for name in own]+[a.validation,future_path,native_path,job_path,
        Path(future['future_roles']),closed_path,archive,history,binary,prlimit,*apis]
    pins=dict(future['pins']);pins.update({str(path):sha(path) for path in paths});verify({'pins':pins})
    scope=('All1620futurejointloggerseedroles,405quartets,135inputs,324aliases. Native--test only; '
        'source inverse, reference alignment/homology, fixed-tip/free-ancestor representation, degree-aware prior '
        'and runtime rooted tree checked. No MCMC, saved joint logger frames, memory repair, mixing, independent '
        'likelihood/root/model/predictor acceptance, GPU or charges. Every native failure/invalid role retained; '
        'full readback/source/artifact/two-journal closure required. Existing sources/attempts untouched.')
    resources=dict(checked_utc=datetime.now(timezone.utc).isoformat(),workers=2,cpus=2,memory_gib=32,swap_gib=0,
        blas_threads=1,per_worker_address_space_gib=12,per_worker_cpu_seconds_cap=600,
        per_worker_wall_seconds_cap=900,per_file_limit_mib=256,output_allowance_gib=16,minimum_free_disk_gib=64,
        historical_full_startup_worker_seconds=seconds,sensitivity_worker_seconds=[seconds/2,seconds*4],
        sum_wall_timeout_worker_hours=1620*900/3600,runtime_uncalibrated=True,finish_eta=None,
        available_memory_gib=psutil.virtual_memory().available/2**30,
        available_disk_gib=shutil.disk_usage(root).free/2**30,native_startup_only=True,posterior_sampling=False,
        gpu=False,new_cost_usd=0,caveat='Different sources/seeds and current shared host load; historical elapsed startup work is not an ETA or measured memory bound. Two12GiB AS caps plus8GiB group headroom; no swap. Output allowance is planning; native per-file limits and free-disk guard are distinct.')
    plan=dict(jobs=str(job_path),output='results/ancestral/full-baliphy-joint-logger-preflight-20261003-v1',
        completion='metadata/baliphy_joint_logger_preflight_completed_20261003.json',
        completion_plan='metadata/baliphy_joint_logger_preflight_completion_plan_20261003.json',
        launch_inventory='metadata/baliphy_joint_logger_preflight_launches_20261003.json',
        dependencies=[],resources=resources,pins=pins,scope=scope)
    with a.plan.open('x') as h:h.write(json.dumps(plan,indent=2)+'\n')
    print(json.dumps(dict(plan=str(a.plan),**resources)))


if __name__=='__main__':main()
