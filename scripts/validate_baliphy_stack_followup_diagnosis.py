#!/usr/bin/env python3
"""Validate both representative stack interventions and all original failures."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import subprocess

from ancestral_chain_attempt import sha
from independent_joint_ancestral_frames import decode
from independent_native_ancestral_alignment import fasta_records
from independent_native_ancestral_topology import match_trees
from independent_short_sampler_outputs_v2 import mapping_rows, strict_json
from reference_measurement_union_sources import verify


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); assert not a.output.exists()
    pp = Path('metadata/baliphy_joint_sampler_qualification_v3_plan_20261003.json')
    plan = json.loads(pp.read_text()); verify(plan['pins'])
    jobs = json.loads(Path(plan['jobs']).read_text()); by = {j['chain']['chain_id']:j for j in jobs}
    root = Path(plan['output']); rows = json.loads((root/'dispositions.json').read_text())
    assert len(rows) == len({r['chain_id'] for r in rows}) == len(by) == 1620
    failed = [r for r in rows if r['exit_code'] != 0]
    assert len(failed) == 24 and Counter(r['family'] for r in failed) == {'OG0000972':24}
    inputs = Counter(r['effective_input_group'] for r in failed)
    assert len(inputs) == 2 and set(inputs.values()) == {12}
    bindings = {str(p):sha(p) for p in [pp, Path(plan['jobs']), root/'dispositions.json', Path(__file__)]}
    for row in failed:
        assert row['exit_code'] == -11 and row['saved_alignments'] == 0 and row['joint_frames'] == []
        assert row['bad_alloc'] is False and row['allocation_warning_lines'] == 0
        job = by[row['chain_id']]; assert job['chain']['proteins'] == 622
        rp = Path(row['native_receipt']); r = json.loads(rp.read_text())
        assert sha(rp) == row['native_receipt_sha256']; bindings[str(rp)] = sha(rp)
        for name, h in r['artifacts'].items():
            fp = rp.parent/name; assert sha(fp) == h; bindings[str(fp)] = h
            if name.endswith(('C1.log', 'C1.P1.fastas', 'C1.P1.site-property-samples.jsonl')):
                assert fp.stat().st_size == 0
    first_path = Path('metadata/baliphy_native_stack_diagnostic_validation_20261003_v1.json')
    first = json.loads(first_path.read_text()); verify(first['source_hashes']); bindings[str(first_path)] = sha(first_path)
    assert first['status'] == 'validated_controlled_native_stack_exhaustion_and_initial_joint_frame'
    cases = [dict(input=by[first['source_role']]['chain']['effective_input_group'],
                  source_role=first['source_role'], default_stack_bytes=first['default_mapped_stack_bytes'],
                  stack_pointer_below_map_bytes=first['stack_pointer_below_mapped_stack_bytes'],
                  corrected_stack_bytes=first['corrected_native_stack_bytes'], initial_frame=first['initial_frame_summary'])]
    configs = {}; logs = {}
    for v, sid, invocation in [(4,24099,'3fd4bf04794f4afcb21bac8db61e3450'),(5,45514,'abc96c728f6b4c4798743b384e36872f')]:
        dr = Path(f'data/software_audits/baliphy-native-segfault-diagnostic-20261003-v{v}')
        cp = dr/'configuration.json'; configs[v] = json.loads(cp.read_text()); verify(configs[v]['pins'])
        rp = dr/'attempt/attempt-0001/receipt.json'; receipt = json.loads(rp.read_text())
        assert receipt['exit_code'] == 0 and receipt['status'] == 'exited_zero_pending_scientific_validation'
        logs[v] = (rp.parent/'stdout.log').read_text()
        for name,h in receipt['artifacts'].items():
            fp = rp.parent/name; assert sha(fp) == h; bindings[str(fp)] = h
        controller = json.loads((dr/'controller.json').read_text()); assert controller['invocation_id'] == invocation
        unit = f'fungal-baliphy-segfault-diagnostic-20261003-v{v}.service'
        raw = subprocess.check_output(['journalctl','--user','-u',unit,'-o','json','--no-pager'],text=True)
        jr = [json.loads(line) for line in raw.splitlines()]
        exact = [r for r in jr if r.get('_PID') == str(controller['pid']) and r.get('_CMDLINE') == ' '.join(controller['cmdline'])]
        assert exact and {r['_SYSTEMD_INVOCATION_ID'] for r in exact} == {invocation}
        resources = [r for r in jr if r.get('USER_INVOCATION_ID') == invocation and r.get('CPU_USAGE_NSEC')]
        assert resources and not any('Failed with result' in r.get('MESSAGE','') or 'Main process exited' in r.get('MESSAGE','') for r in jr if r.get('USER_INVOCATION_ID') == invocation)
        jp = dr/'original-invocation-journal.jsonl'; assert not jp.exists(); jp.write_text(raw)
        transport = dr/'transport.json'; assert not transport.exists()
        transport.write_text(json.dumps(dict(session_id=sid, actual_tool_terminal_exit_code=0, invocation_id=invocation,
            controller=controller, original_process_messages=len(exact), original_completion_resource_records=len(resources),
            native_outcome='SIGSEGV' if v == 4 else 'normal_initial_exit', debugger_exit_code=0),indent=2)+'\n')
        for fp in [cp,rp,dr/'controller.json',dr/'source_job.json',jp,transport,
                   Path(f'metadata/baliphy_native_segfault_diagnostic_resources_20261003_v{v}.json')]:
            bindings[str(fp)] = sha(fp)
        bindings.update(configs[v]['pins'])
    before = configs[4]['command']; after = configs[5]['command']
    assert before[before.index('--args')+1:] == [x for x in after[after.index('--args')+1:] if x != '--stack=67108864']
    assert 'Program received signal SIGSEGV, Segmentation fault.' in logs[4]
    assert 'reg_heap::incremental_evaluate1_changeable_' in logs[4]
    sp = int(re.search(r'^rsp\s+(0x[0-9a-f]+)',logs[4],re.M)[1],16)
    m = re.search(r'^\s*(0x[0-9a-f]+)\s+(0x[0-9a-f]+)\s+(0x[0-9a-f]+).*\[stack\]$',logs[4],re.M)
    low,high,size = (int(m[i],16) for i in [1,2,3]); assert high-low == size == 8*2**20 and 0 < low-sp < 4096
    assert 'SIGSEGV' not in logs[5] and re.search(r'\[Inferior \d+ \(process \d+\) exited normally\]',logs[5])
    source = json.loads(Path('data/software_audits/baliphy-native-segfault-diagnostic-20261003-v5/source_job.json').read_text())
    chain = source['chain']; assert source == by[chain['chain_id']]
    native = after[after.index('--args')+1:]
    assert native[native.index('--seed')+1] == str(chain['seed']) and native[native.index('run')+1] == chain['program']
    folder = Path('data/software_audits/baliphy-native-segfault-diagnostic-20261003-v5/attempt/attempt-0001/independent-chain-1')
    frames = [strict_json(line) for line in (folder/'C1.P1.site-property-samples.jsonl').read_text().splitlines()]
    assert len(frames) == 1 and frames[0]['iter'] == 0
    matched = match_trees(Path(chain['tree']).read_text(),(folder/'runtime-tree.nwk').read_text())
    observed = {label:seq.replace('-','') for label,seq in fasta_records(Path(chain['alignment']).read_text().splitlines()).items()}
    assert len(observed) == 622 and set(observed) == set(matched['tips'])
    bits = {tip:1<<i for i,tip in enumerate(matched['tips'])}; candidates = {}
    for row in mapping_rows(plan['mapping'],chain):
        mask = sum(bits[tip] for tip in strict_json(row['retained_set_json']))
        assert mask == matched['source_labels'][row['source_node']]['mask']
        candidates[row['source_node']] = matched['runtime_index'][mask]['label']
    summary,arrays = decode(frames[0],0,observed,matched['runtime_labels'],candidates)
    assert len(candidates) == 4 and summary['posterior_qualified'] is summary['scientific_eligibility'] is False
    cases.append(dict(input=chain['effective_input_group'],source_role=chain['chain_id'],default_stack_bytes=size,
                      stack_pointer_below_map_bytes=low-sp,corrected_stack_bytes=64*2**20,initial_frame=summary))
    assert {c['input'] for c in cases} == set(inputs)
    verify(bindings)
    result = dict(status='validated_two_input_stack_exhaustion_and_initial_joint_frames',checked_utc=datetime.now(timezone.utc).isoformat(),
        full_original_roles=1620, original_successful_roles=1596, original_failed_roles=24,
        failed_input_group_counts=dict(inputs), representative_cases=cases, source_hashes=bindings,
        native_model_prior_input_seed_preserved=True, existing_jobs_restarted=False, posterior_qualified=False,
        full_twenty_iteration_correction_qualified=False, diagnostic_samples_excluded_from_posteriors=True,
        scope='Full24failure accounting with two causal initial-state representatives. Default8MiB stack reproduces SIGSEGV below mapping;64MiB exits normally with strict independently decoded622tip frames. Each native model,prior,input,seed and other caps unchanged between interventions. This does not prove20iteration success for all24 roles, full computational closure, long resources or adequate posteriors.')
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','representative_cases']},indent=2))


if __name__ == '__main__': main()
