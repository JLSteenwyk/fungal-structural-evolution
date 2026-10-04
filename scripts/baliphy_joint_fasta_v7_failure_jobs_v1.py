"""Exact source-preserving construction of all24 original failed-role comparisons."""
from collections import Counter,defaultdict
import copy
import hashlib
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from baliphy_joint_fasta_lines_v7 import transform,reverse
from reference_measurement_union_sources import bind,verify

SUFFIX='-joint-fasta-v7-failure-comparison'


def failed_sources():
    proof_path=Path('metadata/current_analysis_error_recheck_20261004_user_1224.json')
    plan_path=Path('metadata/baliphy_scalar_v6_sampler_execution_plan_20261004_v1.json')
    proof=json.loads(proof_path.read_text());plan=json.loads(plan_path.read_text())
    jobs_path=Path(plan['jobs']);jobs=json.loads(jobs_path.read_text());assert len(jobs)==1620
    by_id={j['chain']['chain_id']:j for j in jobs};assert len(by_id)==1620
    pins={};originals=[]
    for name,h in proof['source_hashes'].items():
        if '/full-baliphy-scalar-v6-short-sampler-' not in name or not name.endswith('/receipt.json'):
            continue
        bind(pins,name,h);native=json.loads(Path(name).read_text());assert native['exit_code']==-11
        cid=Path(name).parent.parent.name;original=by_id[cid]
        encoded=json.dumps(original['config'],sort_keys=True,separators=(',',':'),allow_nan=False)
        assert hashlib.sha256(encoded.encode()).hexdigest()==native['configuration_sha256']
        for artifact,digest in native['artifacts'].items():bind(pins,Path(name).parent/artifact,digest)
        for source,digest in original['config']['pins'].items():bind(pins,source,digest)
        originals.append(original)
    for p in [proof_path,plan_path,jobs_path,Path(plan['mapping'])]:bind(pins,p)
    assert len(originals)==24 and len({j['chain']['chain_id'] for j in originals})==24
    assert {j['chain']['family'] for j in originals}=={'OG0000972'}
    assert {j['chain']['proteins'] for j in originals}=={622}
    assert Counter(j['chain']['prior_label'] for j in originals)=={'broad':8,'centered':8,'package':8}
    groups=defaultdict(list)
    for j in originals:groups[(j['chain']['effective_input_group'],j['chain']['prior_label'])].append(j['chain']['chain'])
    assert len(groups)==6 and all(sorted(v)==[1,2,3,4] for v in groups.values())
    assert len({j['chain']['effective_input_group'] for j in originals})==2
    assert all(j['memory_reservation_bytes']==48*2**30 for j in originals)
    verify(pins)
    return sorted(originals,key=lambda j:j['chain']['chain_id']),pins,plan['mapping']


def candidate_job(original,program):
    program=Path(program).resolve();old=original['chain']
    assert sha(old['program'])==old['program_sha256']
    assert reverse(program.read_text())==Path(old['program']).read_text()
    job=copy.deepcopy(original)
    job['chain'].update(chain_id=old['chain_id']+SUFFIX,program=str(program),program_sha256=sha(program))
    command=job['config']['command']
    assert command[command.index('run')+1]==old['program']
    assert command[command.index('--seed')+1]==str(old['seed'])
    assert command[command.index('--iterations')+1]=='20'
    assert '--as='+str(48*2**30) in command
    assert '--fsize='+str(2*2**30) in command
    assert not any(v.startswith('--stack=') for v in command)
    command[command.index('run')+1]=str(program)
    command.insert(command.index('--'),'--stack='+str(8*2**20)+':unlimited')
    assert job['config']['pins'].pop(old['program'])==old['program_sha256']
    job['config']['pins'][str(program)]=sha(program)
    job['source_v6_chain_id']=old['chain_id'];job['source_v6_seed']=old['seed']
    return job


def build_jobs(originals,root):
    root=Path(root).resolve();root.mkdir(exist_ok=False)
    jobs=[];models={}
    for original in originals:
        old=original['chain'];key=old['effective_input_group']+'-'+old['prior_label']
        program=root/(key+'.hs')
        text=transform(Path(old['program']).read_text())
        if key not in models:program.write_text(text);models[key]=program
        else:assert program.read_text()==text
        jobs.append(candidate_job(original,program))
    assert len(models)==6
    return jobs


def validate_jobs(jobs,originals):
    assert len(jobs)==len(originals)==24
    by_id={j['chain']['chain_id']:j for j in originals}
    assert len(by_id)==24
    assert {j['source_v6_chain_id'] for j in jobs}==set(by_id)
    assert len({j['chain']['chain_id'] for j in jobs})==24
    assert len({j['chain']['seed'] for j in jobs})==24
    for job in jobs:
        original=by_id[job['source_v6_chain_id']]
        expected=candidate_job(original,job['chain']['program'])
        assert job==expected,'Altered original scope/configuration or candidate model'
        verify(job['config']['pins'])
    assert len({j['chain']['program'] for j in jobs})==6
    return dict(original_failed_roles=24,effective_inputs=2,prior_quartets=6,
                prior_role_counts={'broad':8,'centered':8,'package':8},chains_per_quartet=4,
                original_seeds_preserved=24,original_models_reversibly_preserved=6,
                native_as_bytes=48*2**30,stack_soft_bytes=8*2**20,stack_hard_unlimited=True,
                original_cpu_wall_file_caps_preserved=True,iterations=20,
                no_original_failure_restarted=True)
