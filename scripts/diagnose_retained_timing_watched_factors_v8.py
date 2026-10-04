#!/usr/bin/env python3
"""Replay original probes with two unprotected source-factor hardware watchpoints.

Calls the preserved original probe, including its two guards, likelihood
evaluations and independent replays. Never optimizes or accepts a biological
model. Original failed stages remain immutable and stopped.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
import os
import signal
import psutil

import numpy as np
from scipy import sparse

from ancestral_chain_attempt import sha, write_json
from full_expanded_model_design_sources import digest
from full_retained_shared_entity_timing import probe
from prepare_baliphy_scalar_v6_preflight_v1 import project_sources
from reference_measurement_union_sources import bind, verify


def array_hash(value):
    h = hashlib.sha256()
    h.update(json.dumps([str(value.dtype), list(value.shape)]).encode())
    h.update(np.ascontiguousarray(value).tobytes())
    return h.hexdigest()


def input_hashes(source, operators, x, y):
    values = dict(labels=array_hash(source['labels']), design=array_hash(x), response=array_hash(y),
        factors={name:array_hash(value) for name,value in source['factors'].items()}, operators={})
    for mode, bank in operators.items():
        values['operators'][mode] = {name:dict(shape=list(z.shape), data=array_hash(z.data),
            indices=array_hash(z.indices), indptr=array_hash(z.indptr)) for name,z in bank.items()}
    return values


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    a = p.parse_args(); assert not a.receipt.exists()
    plan = json.loads(a.plan.read_text()); verify(plan['pins'])
    pins = dict(plan['pins']); bind(pins, a.plan)
    modules = project_sources(pins, [Path(__file__)])
    diagnostic_path = Path(plan['diagnostic']); diagnostic = json.loads(diagnostic_path.read_text())
    assert diagnostic['status'] == 'completed_exact_failed_cohort_covariance_diagnostic_v1'
    assert (diagnostic['records'], diagnostic['selected_groups'], diagnostic['candidate_census_rows']) == (22881,40,1200)
    assert diagnostic['guard_status_counts'] == {'passed_unchanged_original_guard':40}
    for name in ['diagnostic_transport', 'reference_transport']:
        path = Path(plan[name]); transport = json.loads(path.read_text())
        assert transport['status'] == 'verified_original_retained_covariance_diagnostic_wait_zero'
        assert transport['actual_tool_terminal_exit_code'] == 0 and transport['entire_terminal_payload_matched']
        verify(transport['source_hashes'])
        bind(pins, path)
    reference_path = Path(plan['reference']); reference = json.loads(reference_path.read_text())
    assert reference['status'] == 'completed_extended_precision_failed_cohort_raw_reference_v1'
    assert reference['producer_sha256'] == sha(diagnostic_path)
    assert reference['selected_groups'] == 40 and reference['raw_reference_contexts'] == 10
    assert all(r[k]['passes_original_allclose'] for r in reference['groups']
               for k in ['reference_vs_fresh','reference_vs_inherited','reference_vs_latent'])
    for path,h in diagnostic['artifacts'].items(): bind(pins,path,h)
    verify(pins)
    original_plan_path = Path(plan['original_timing_plan'])
    timing = json.loads(original_plan_path.read_text())
    assert timing['scaled_variance_points'] == [0.,1.]
    fit_path = Path(timing['fit_plan']); assert sha(fit_path) == timing['fit_plan_sha256']
    fit = json.loads(fit_path.read_text())
    assert fit['launch_state'] == 'not_launched_or_queued' and not Path(fit['output']).exists()
    assert fit['independent_backend'] == 'component_spectral_v1'
    for path in [diagnostic_path, reference_path, original_plan_path, fit_path]: bind(pins,path)
    source_root = Path(json.loads(Path(diagnostic['plan']).read_text())['output'])
    source = dict(labels=np.load(source_root/'component_labels.npy'), factors={})
    rows = np.arange(diagnostic['records'])
    assert source['labels'].shape == rows.shape
    details = {}; operators = {}
    for item in diagnostic['groups']:
        group = item['group_id']; details[group] = json.loads((source_root/'groups'/(group+'.json')).read_text())
        identity = item['representative']; mode = identity['loading_mode']; tree = identity['tree']
        names = details[group]['retained_source_audit']['retained_kernel_names']
        if tree not in source['factors']:
            source['factors'][tree] = np.load(source_root/(tree+'-factor.npy'))
            assert source['factors'][tree].shape == (22881,301) and source['factors'][tree].dtype == np.float64
        if mode not in operators: operators[mode] = {name:sparse.load_npz(source_root/(mode+'-'+name+'.npz')) for name in names[1:-1]}
        assert list(operators[mode]) == names[1:-1]
    assert len(source['factors']) == 5 and len(operators) == 2 and len(details) == 40
    root = Path(plan['output']); root.mkdir(exist_ok=False); (root/'groups').mkdir()
    layout = {name:dict(shape=list(value.shape),dtype=str(value.dtype),strides=list(value.strides),
        pointer=int(value.__array_interface__['data'][0]),nbytes=value.nbytes,
        c_contiguous=bool(value.flags.c_contiguous),f_contiguous=bool(value.flags.f_contiguous))
        for name,value in source['factors'].items()}
    assert all(not np.shares_memory(a,b) for i,a in enumerate(source['factors'].values())
               for j,b in enumerate(source['factors'].values()) if i<j)
    write_json(root/'factor_buffer_layout.json',layout)
    inferior = psutil.Process()
    write_json(root/'watchpoint_inferior_identity.json', dict(pid=inferior.pid,
        created=inferior.create_time(), cmdline=inferior.cmdline()))
    targets = {}
    for tree in plan['watch_trees']:
        value = source['factors'][tree]
        address = layout[tree]['pointer'] + 8513*value.strides[0] + 116*value.strides[1]
        assert address % 8 == 0
        targets[tree] = dict(address=address, bytes=8, row=8513, column=116,
            initial_uint64=str(int(value.view(np.uint64)[8513,116])),
            factor_base=layout[tree]['pointer'], factor_bytes=value.nbytes)
    write_json(Path(plan['watchpoint_metadata']), dict(targets=targets,
        inferior_pid=inferior.pid, inferior_created=inferior.create_time(),
        process_local_read_only_protection=False, artificial_control=False))
    print('WATCHPOINT_TARGETS_READY', flush=True)
    os.kill(os.getpid(), signal.SIGSTOP)
    started = time.monotonic(); results = []; mutation = None
    ordered = sorted(diagnostic['groups'],key=lambda r:r['group_id'])
    write_json(root/'ordered_group_ids.json', [r['group_id'] for r in ordered])
    for number,item in enumerate(ordered,1):
        group = item['group_id']; detail = details[group]; identity = item['representative']
        with np.load(source_root/'groups'/(group+'.npz'),allow_pickle=False) as saved:
            x = saved['design']; y = saved['response']
        representative = dict(identity=identity,matrix=x,response=y,rank=item['selection_rank'],
            source_audit=detail['retained_source_audit'], original_audit=detail['original_source_audit'],
            certificate=detail['certificate'],key=(identity['loading_mode'],identity['tree'],identity['method'],identity['outcome']))
        before = input_hashes(source,operators,x,y)
        write_json(root/'current_group.json',dict(group_id=group,sequence=number,before=before,scientific_eligibility=False))
        try:
            result = probe(source,fit,timing,rows,operators[identity['loading_mode']],representative,item['eligible_candidates'],group)
        except (AssertionError,ValueError,ArithmeticError,np.linalg.LinAlgError) as error:
            result = dict(group_id=group,representative=identity,eligible_candidates=item['eligible_candidates'],
                status='original_probe_rejected_in_ordered_diagnostic',error_type=type(error).__name__,
                error_message=str(error),scientific_eligibility=False)
        after = input_hashes(source,operators,x,y)
        capture = {}
        second_after = input_hashes(source,operators,x,y)
        if before != after:
            for tree,changed in source['factors'].items():
                if before['factors'][tree] == after['factors'][tree]: continue
                baseline = np.load(source_root/(tree+'-factor.npy'))
                observed = changed.copy()
                assert baseline.dtype == observed.dtype == np.float64 and baseline.shape == observed.shape
                old_bits = np.ascontiguousarray(baseline).view(np.uint64)
                new_bits = np.ascontiguousarray(observed).view(np.uint64)
                offsets = np.argwhere(old_bits != new_bits)
                snapshots = root/'captured_factors'; snapshots.mkdir(exist_ok=True)
                snapshot = snapshots/(group+'-'+tree+'-observed.npy'); np.save(snapshot,observed)
                with np.errstate(over='ignore',invalid='ignore'):
                    delta = abs(observed-baseline)
                finite = delta[np.isfinite(delta)]
                samples = []
                for row,column in offsets[:1000]:
                    old = int(old_bits[row,column]); new = int(new_bits[row,column]); xor = old^new
                    samples.append(dict(row=int(row),column=int(column),baseline_hex=float(baseline[row,column]).hex(),
                        observed_hex=float(observed[row,column]).hex(),baseline_bits=str(old),observed_bits=str(new),
                        xor_hex=hex(xor),changed_bits=xor.bit_count()))
                capture[tree] = dict(changed_entries=len(offsets),sampled_changed_entries=len(samples),samples=samples,
                    current_probe_tree=identity['tree'],current_probe_uses_changed_tree=identity['tree']==tree,
                    baseline_array_hash=array_hash(baseline),observed_array_hash=array_hash(observed),
                    repeated_current_array_hash=array_hash(changed),
                    baseline_matches_before_hash=array_hash(baseline)==before['factors'][tree],
                    observed_values_all_finite=bool(np.isfinite(observed).all()),
                    maximum_finite_absolute_difference=float(finite.max()) if finite.size else None,
                    observed_snapshot=str(snapshot),observed_snapshot_sha256=sha(snapshot))
        record = dict(sequence=number,group_id=group,original_probe_result=result,
            input_hashes_before=before,input_hashes_after=after,second_input_hashes_after=second_after,
            factor_value_change_capture=capture,inputs_preserved=before==after,
            scientific_eligibility=False)
        write_json(root/'groups'/(group+'.json'),record); results.append(record)
        print('ordered_original_timing_probe',number,'/40',group,result['status'],'inputs_preserved',before==after,flush=True)
        if before != after:
            mutation = dict(sequence=number,group_id=group,factor_value_change_capture=capture)
            break
    pending = [r['group_id'] for r in ordered[len(results):]]
    verify(pins)
    artifacts = {str(path):sha(path) for path in root.rglob('*') if path.is_file()}
    result = dict(status='completed_watched_factor_timing_probe_diagnostic_v8',
        checked_utc=datetime.now(timezone.utc).isoformat(),plan=str(a.plan),plan_sha256=sha(a.plan),
        cohort_id=diagnostic['cohort_id'],records=22881,planned_groups=40,executed_groups=len(results),
        unattempted_group_ids=pending,detected_input_mutation=mutation,
        probe_status_counts=dict(Counter(r['original_probe_result']['status'] for r in results)),
        inputs_preserved_for_all_attempted_groups=all(r['inputs_preserved'] for r in results),
        elapsed_seconds=time.monotonic()-started,transitive_project_source_modules=len(modules),
        source_hashes=pins,artifacts=artifacts,
        full_original_source_verification_inherited=True,full_original_source_archive_rehashed=False,
        original_prior_ten_cohort_probe_sequence_replayed=False,
        scientific_eligibility=False,production_tolerance_changed=False,biological_fits=0,gpu=False,
        original_jobs_restarted=False, process_local_read_only_protection=False, hardware_watchpoints_requested=2,
        scope='Ordered40original selected timing probes on closed exact cohort11 exports, including '
              'preserved original guard/primary likelihood/independent qualification/likelihood points. '
              'Shared cohort factors/operators and designs/responses hashed before/after each probe; repeated after hashes, nonoverlapping factor-buffer layouts and affected factor snapshots/bitwise differences retained; '
              'mutation stops further attempts and retains all unattempted identities. Diagnostics '
              'only, not whole prior10cohort execution-state replay, fit/source repair, tolerance '
              'relaxation or biological acceptance. Full source verification is inherited from '
              'the closed exact-cohort producer; consumed exported artifact bindings rehash.')
    with a.receipt.open('x') as f: json.dump(result,f,indent=2); f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']},indent=2))


if __name__ == '__main__':
    main()
