#!/usr/bin/env python3
"""Validate controlled native stack reproduction and the initial joint frame."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re

from ancestral_chain_attempt import sha
from independent_joint_ancestral_frames import decode
from independent_native_ancestral_alignment import fasta_records
from independent_native_ancestral_topology import match_trees
from independent_short_sampler_outputs_v2 import mapping_rows, strict_json
from reference_measurement_union_sources import verify


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    assert not a.output.exists()
    plan_path=Path('metadata/baliphy_joint_sampler_qualification_v3_plan_20261003.json')
    plan=json.loads(plan_path.read_text());verify(plan['pins'])
    jobs=json.loads(Path(plan['jobs']).read_text())
    cid='acfc837699b7ad7b47bda478356d778fa35864ba0a2af4bfe55096cf9e8b65a5-broad-joint-cjson-v5-chain1'
    job=next(j for j in jobs if j['chain']['chain_id']==cid);chain=job['chain'];bindings=dict(job['config']['pins'])
    rp0=Path(plan['output'])/'chains'/(cid+'.json');original=json.loads(rp0.read_text())
    assert original['exit_code']==-11 and original['saved_alignments']==0
    bindings.update({str(p):sha(p) for p in [plan_path,rp0,Path(plan['jobs']),Path(__file__)]})
    configs={};logs={};receipts={}
    for v in [1,2,3]:
        root=Path(f'data/software_audits/baliphy-native-segfault-diagnostic-20261003-v{v}')
        configs[v]=json.loads((root/'configuration.json').read_text())
        verify(configs[v]['pins']);bindings.update(configs[v]['pins'])
        folder=root/'attempt/attempt-0001';rp=folder/'receipt.json';receipts[v]=json.loads(rp.read_text())
        assert receipts[v]['status']!='timeout'
        logs[v]=(folder/'stdout.log').read_text()
        for p in [root/'configuration.json',rp,Path(f'metadata/baliphy_native_segfault_diagnostic_resources_20261003_v{v}.json')]:
            bindings[str(p)]=sha(p)
        for name,h in receipts[v]['artifacts'].items():
            fp=folder/name;assert sha(fp)==h;bindings[str(fp)]=h
    assert 'failed to set the FSIZE resource limit: Operation not permitted' in Path(
        'data/software_audits/baliphy-native-segfault-diagnostic-20261003-v1/attempt/attempt-0001/stderr.log').read_text()
    assert receipts[1]['exit_code']==1
    # A debugger can exit zero while its inferior crashed. Native outcome is
    # established from explicit debugger events and complete original outputs.
    assert 'Program received signal SIGSEGV, Segmentation fault.' in logs[2]
    assert 'reg_heap::incremental_evaluate1_changeable_' in logs[2]
    sp=int(re.search(r'^rsp\s+(0x[0-9a-f]+)',logs[2],re.M)[1],16)
    m=re.search(r'^\s*(0x[0-9a-f]+)\s+(0x[0-9a-f]+)\s+(0x[0-9a-f]+).*\[stack\]$',logs[2],re.M)
    low,high,size=(int(m[i],16) for i in [1,2,3])
    assert high-low==size==8*2**20 and 0<low-sp<4096
    before=configs[2]['command'];after=configs[3]['command']
    assert '--args' in before and '--args' in after
    native_before=before[before.index('--args')+1:]
    native_after=after[after.index('--args')+1:]
    assert native_before==[p for p in native_after if p!='--stack=67108864']
    assert native_before[native_before.index('--iterations')+1]=='0'
    assert native_before[native_before.index('--seed')+1]==str(chain['seed'])
    assert native_before[native_before.index('run')+1]==chain['program']
    assert 'SIGSEGV' not in logs[3] and re.search(r'\[Inferior \d+ \(process \d+\) exited normally\]',logs[3])
    assert receipts[3]['exit_code']==0
    folder=Path('data/software_audits/baliphy-native-segfault-diagnostic-20261003-v3/attempt/attempt-0001/independent-chain-1')
    frames=[strict_json(line) for line in (folder/'C1.P1.site-property-samples.jsonl').read_text().splitlines()]
    assert len(frames)==1 and frames[0]['iter']==0
    matched=match_trees(Path(chain['tree']).read_text(),(folder/'runtime-tree.nwk').read_text())
    observed={label:seq.replace('-','') for label,seq in fasta_records(Path(chain['alignment']).read_text().splitlines()).items()}
    assert len(observed)==chain['proteins']==622 and set(observed)==set(matched['tips'])
    bits={tip:1<<i for i,tip in enumerate(matched['tips'])};candidates={}
    for row in mapping_rows(plan['mapping'],chain):
        tips=strict_json(row['retained_set_json']);mask=sum(bits[tip] for tip in tips)
        assert mask==matched['source_labels'][row['source_node']]['mask']
        candidates[row['source_node']]=matched['runtime_index'][mask]['label']
    assert len(candidates)==len(set(candidates.values()))==4
    summary,arrays=decode(frames[0],0,observed,matched['runtime_labels'],candidates)
    assert summary['scientific_eligibility'] is summary['posterior_qualified'] is False
    bindings[str(Path(plan['mapping']))]=sha(plan['mapping'])
    failed_path=Path('metadata/baliphy_joint_sampler_failure_audit_20261003_v1.json')
    failed=json.loads(failed_path.read_text());verify(failed['source_hashes']);bindings[str(failed_path)]=sha(failed_path)
    assert failed['failed_exit_code_counts']=={'-11':12}
    assert failed['failed_input_group_counts']=={chain['effective_input_group']:12}
    assert failed['cgroup_initial_memory_events']['oom_kill']==failed['cgroup_final_memory_events']['oom_kill']==0
    assert failed['cgroup_initial_memory_events']['oom']==failed['cgroup_final_memory_events']['oom']==0
    for row in failed['failed_role_rows']:
        resource=row['resource_observations'];assert resource['live_observations']>0
        assert resource['maximum_reported_memory_bytes']['VmPeak']<48*2**30
    verify(bindings)
    result=dict(status='validated_controlled_native_stack_exhaustion_and_initial_joint_frame',
        checked_utc=datetime.now(timezone.utc).isoformat(),source_role=cid,representative_default_stack_signal='SIGSEGV',
        default_mapped_stack_bytes=size,stack_pointer_below_mapped_stack_bytes=low-sp,
        corrected_native_stack_bytes=64*2**20,corrected_native_initial_state_exit='normal',
        corrected_initial_joint_frames=1,initial_frame_summary=summary,
        projection_array_keys=sorted(arrays),original_failed_roles=12,original_failed_input_groups=1,
        all_failed_originals_sampled_below_address_cap=True,no_cgroup_oom_event_in_observed_prefix=True,
        unchanged_model_prior_alignment_tree_seed_verified=True,posterior_qualified=False,
        full_twenty_iteration_correction_qualified=False,existing_jobs_restarted=False,source_hashes=bindings,
        scope='One representative sourceinput/seed:8MiB stack fails at evaluator with SP144bytes below mapped stack;64MiB stack exits normally at iterations0. Model/priors/input/tree/nativeAS48GiB/file2GiB/CPU5674sec unchanged. Initial622tip joint frame independently decoded with strict rates,states,coordinates and clade/candidate checks. Twelve original failures and all their artifacts retained; sampled telemetry below AS cap and no prefix OOM event do not guarantee final memory or identify every failure causally. Diagnostic samples are excluded from posterior ensembles. Complete20iteration follow-up/full cohort closure remains required.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','initial_frame_summary','projection_array_keys']},indent=2))
    print('initial_frame',json.dumps(summary))


if __name__=='__main__':main()
