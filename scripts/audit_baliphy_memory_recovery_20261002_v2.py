#!/usr/bin/env python3
"""Verify completed recovery attempts and retain the complete original chain grid.

The original native readback is replayed for new exit-zero attempts. This is a
separate invocation of the same integrity parser, not an independent parser or
a posterior convergence assessment. Original closed outputs are hash-checked.
"""
import argparse
from collections import Counter, defaultdict
import copy
import fcntl
import hashlib
import json
from pathlib import Path
import shutil
import time

from ancestral_chain_attempt import live_group
from background_measurement_union_sources import closed_source
from readback_independent_baliphy_chain import readback
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha

CHECKED = 'all_saved_alignments_and_candidate_nodes_checked'


def encoded_sha(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def changed_configuration(original, address_space):
    config = copy.deepcopy(original)
    assert config['command'][0] == '/usr/bin/prlimit'
    flags = [i for i, part in enumerate(config['command']) if part.startswith('--as=')]
    assert len(flags) == 1 and config['command'][flags[0]] == '--as=12884901888'
    assert address_space == 48 * 2**30
    config['command'][flags[0]] = '--as=' + str(address_space)
    return config


def build_overlay(jobs, originals, recovery_rows, failed_ids):
    """Select each whole attempt once; never join chains or drop failed groups."""
    index = {j['chain']['chain_id']: j for j in jobs}
    original = {r['chain_id']: r for r in originals}
    recovery = {r['chain_id']: r for r in recovery_rows}
    assert len(jobs) == len(index) == len(originals) == len(original) == 1620
    assert set(index) == set(original)
    failed = {cid for cid, row in original.items() if row['status'] != CHECKED}
    assert failed == set(failed_ids) == set(recovery) and len(recovery_rows) == len(recovery) == 3
    assert all(original[cid]['status'] == 'failed' for cid in failed)
    rows = []; groups = defaultdict(list)
    for cid, job in sorted(index.items()):
        chain, config = job['chain'], job['config']
        group = config['model_input_identity']
        assert group == chain['effective_input_group'] + '-' + chain['prior_label']
        assert cid == group + '-chain' + str(chain['chain'])
        assert chain['prior_label'] in {'package', 'centered', 'broad'}
        assert config['seed'] == chain['seed']
        command = config['command']
        assert command[command.index('--seed') + 1] == str(chain['seed'])
        assert command[command.index('--iterations') + 1] == '1000'
        assert chain['original_configuration_ids']
        original_row = original[cid]
        selected = copy.deepcopy(original_row); origin = 'original_checked_attempt'
        retry = recovery.get(cid)
        if retry is not None:
            assert retry['original_failed_attempt'] == original_row['receipt']
            assert retry['original_failed_attempt_sha256'] == original_row['receipt_sha256']
            assert retry['status'] in {CHECKED, 'failed'}
            assert retry['native_bad_alloc'] is (retry['status'] == 'failed')
            selected = dict(chain_id=cid, status=retry['status'], receipt=retry['new_attempt'],
                            receipt_sha256=retry['new_attempt_sha256'])
            if retry['status'] == CHECKED:
                assert retry['sample_audit'] and retry['sample_audit_sha256']
                selected.update(sample_audit=retry['sample_audit'],
                                sample_audit_sha256=retry['sample_audit_sha256'])
                origin = 'new_whole_same_seed_integrity_checked_attempt'
            else:
                assert 'sample_audit' not in retry and 'sample_audit_sha256' not in retry
                origin = 'new_failed_attempt_original_failure_retained'
        assert selected['status'] in {CHECKED, 'failed'}
        if selected['status'] == CHECKED:
            assert selected['sample_audit'] and selected['sample_audit_sha256']
        row = dict(chain=copy.deepcopy(chain), model_input_identity=group,
                   original_disposition=copy.deepcopy(original_row),
                   recovery_disposition=copy.deepcopy(retry), selected_disposition=selected,
                   selection_origin=origin, scientific_eligibility=False)
        rows.append(row); groups[group].append(row)
    assert len(groups) == 405
    quartets = []; effective = defaultdict(set)
    for group, members in sorted(groups.items()):
        assert len(members) == 4
        assert {m['chain']['chain'] for m in members} == {1, 2, 3, 4}
        assert len({m['chain']['seed'] for m in members}) == 4
        for field in ['effective_input_group', 'prior_label', 'family', 'alignment',
                      'alignment_sha256', 'tree', 'tree_sha256', 'program', 'program_sha256',
                      'original_configuration_ids']:
            assert all(m['chain'][field] == members[0]['chain'][field] for m in members)
        effective[members[0]['chain']['effective_input_group']].add(members[0]['chain']['prior_label'])
        missing = [m['chain']['chain_id'] for m in members if m['selected_disposition']['status'] != CHECKED]
        quartets.append(dict(model_input_identity=group,
            chain_ids=[m['chain']['chain_id'] for m in members],
            status='complete_integrity_checked_quartet' if not missing else 'unresolved_failed_chain',
            unresolved_chain_ids=missing, scientific_eligibility=False))
    assert len(effective) == 135 and all(priors == {'package', 'centered', 'broad'} for priors in effective.values())
    counts = Counter(r['selected_disposition']['status'] for r in rows)
    complete = sum(q['status'] == 'complete_integrity_checked_quartet' for q in quartets)
    summary = dict(full_native_chains=len(rows), full_quartets=len(quartets),
        original_checked_chains=sum(r['status'] == CHECKED for r in originals),
        original_failed_chains=len(failed), recovery_attempts=len(recovery_rows),
        recovered_integrity_checked_chains=sum(r['status'] == CHECKED for r in recovery_rows),
        unresolved_failed_chains=counts['failed'], selected_checked_chains=counts[CHECKED],
        complete_quartets=complete, unresolved_quartets=len(quartets) - complete)
    return rows, quartets, summary


def launch_plan_path(launch):
    command = launch['cmdline']
    assert command.count('--plan') == 1
    position = command.index('--plan')
    assert position + 1 < len(command)
    command_plan = command[position + 1]
    if 'plan' in launch:
        assert launch['plan'] == command_plan
    return command_plan


def audit_attempt(job, record, original_root, recovery_root, iterations, mapping, address_space, bindings):
    cid = job['chain']['chain_id']; config = changed_configuration(job['config'], address_space)
    rp = Path(record['new_attempt']); bind(bindings, rp, record['new_attempt_sha256'])
    assert rp.parent.parent == (recovery_root / cid).resolve() and rp.parent.name == 'attempt-0001'
    assert record['original_failed_attempt'] == str((original_root / cid / 'attempt-0001/receipt.json').resolve())
    old = Path(record['original_failed_attempt']); bind(bindings, old, record['original_failed_attempt_sha256'])
    original = json.loads(old.read_text()); assert original['status'] == 'failed' and original['exit_code'] != 0
    for name, digest in original['artifacts'].items(): bind(bindings, old.parent / name, digest)
    assert 'std::bad_alloc' in (old.parent / 'stderr.log').read_text()
    expected_hash = encoded_sha(config); receipt = json.loads(rp.read_text())
    assert receipt['configuration_sha256'] == expected_hash
    for path, expected in [(rp.parent.parent / 'configuration.json', config),
                           (rp.parent / 'command.json', config['command'])]:
        assert json.loads(path.read_text()) == expected; bind(bindings, path)
    assert not any('{attempt}' in part for part in config['command'])
    for name, digest in receipt['artifacts'].items():
        assert Path(name).is_relative_to('.') and not Path(name).is_absolute() and '..' not in Path(name).parts
        bind(bindings, rp.parent / name, digest)
    for directory in [old.parent, rp.parent]:
        process = json.loads((directory / 'process.json').read_text())
        assert not live_group(process['pgid'])
        assert fingerprint(dict(pid=process['pid'], created=process['created'], cmdline=process['command'])) is None
    process = json.loads((rp.parent / 'process.json').read_text())
    assert process['command'] == config['command'] and process['pgid'] == process['pid']
    stderr = rp.parent / 'stderr.log'; text = stderr.read_text(); bad_alloc = 'std::bad_alloc' in text
    assert record['native_bad_alloc'] is bad_alloc
    audit = dict(chain_id=cid, configuration_sha256=expected_hash, original_seed_preserved=True,
                 native_process=process, native_process_group_terminal=True,
                 allocation_warnings=text.count('Allocation failed in sample_tri_multi'),
                 native_bad_alloc=bad_alloc, exit_code=receipt['exit_code'], status=record['status'],
                 stderr=str(stderr), stderr_sha256=sha(stderr))
    if record['status'] == CHECKED:
        assert receipt['exit_code'] == 0 and receipt['status'] == 'exited_zero_pending_scientific_validation' and not bad_alloc
        ap = Path(record['sample_audit']); bind(bindings, ap, record['sample_audit_sha256'])
        saved = json.loads(ap.read_text())
        repeated = readback(job['chain'], rp, iterations, mapping)
        assert repeated == saved
        assert len(saved['candidate_samples']) == 4 * (iterations // 10 + 1)
        audit.update(replayed_full_saved_alignment_integrity=True,
                     saved_alignment_samples=iterations // 10 + 1,
                     candidate_node_samples=len(saved['candidate_samples']))
    else:
        assert record['status'] == receipt['status'] == 'failed' and receipt['exit_code'] != 0 and bad_alloc
        assert not list(rp.parent.parent.glob('attempt-*-sample-audit.json'))
        audit['replayed_full_saved_alignment_integrity'] = False
    for path, digest in config['pins'].items(): bind(bindings, path, digest)
    return audit


def run(path):
    started = time.monotonic(); plan = json.loads(path.read_text()); bindings = dict(plan['pins']); bind(bindings, path)
    verify(bindings)
    root = Path(plan['output']); root.mkdir(exist_ok=False)
    lock = (root / 'stage.lock').open('x'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    assert shutil.disk_usage(root).free >= plan['resources']['minimum_free_disk_gib'] * 2**30
    recovery_plan = json.loads(Path(plan['recovery_plan']).read_text())
    source_plan = json.loads(Path(recovery_plan['producer_plan']).read_text())
    recovery_launch = json.loads(Path(plan['recovery_launch']).read_text())
    assert recovery_launch['plan'] == plan['recovery_plan']
    assert recovery_launch['actual_cgroup_limits'] == {
        'cpu.max': '100000 100000', 'memory.max': str(64 * 2**30), 'memory.swap.max': '0'}
    assert source_plan['iterations'] == 1000 and source_plan['diagnostics']['burn_in_fractions'] == [0.25, 0.5]
    horizon = closed_source(plan['horizon_completion'], 'complete_verified_baliphy_initial_horizon_accounting',
        'complete_verified_baliphy_initial_horizon_accounting_archive', 2, bindings)
    assert horizon['chains'] == 1620 and horizon['checked_chains'] == 1617 and horizon['failed_chains'] == 3
    for p, d in recovery_plan['pins'].items(): bind(bindings, p, d)
    journals = []
    for launch_path in [*recovery_plan['original_launches'], plan['recovery_launch']]:
        launch = json.loads(Path(launch_path).read_text()); launch['launch'] = launch_path
        launch_plan = launch_plan_path(launch)
        assert sha(launch_plan) == launch['plan_sha256']; bind(bindings, launch_path); bind(bindings, launch_plan)
        journal = journal_terminal(launch); journal['status'] = 'verified_original_terminal_resource_journal'
        journals.append(journal)
    original_root = Path(source_plan['output']); recovery_root = Path(recovery_plan['output'])
    original_path = original_root / 'receipt.json'; original = json.loads(original_path.read_text())
    retry_path = recovery_root / 'receipt.json'; retry = json.loads(retry_path.read_text())
    bind(bindings, original_path); bind(bindings, retry_path)
    assert original['status'] == 'all_initial_horizon_dispositions_recorded' and original['plan_sha256'] == sha(recovery_plan['producer_plan'])
    assert retry['status'] == 'complete_baliphy_memory_recovery_pending_full_readback' and retry['plan_sha256'] == sha(plan['recovery_plan'])
    assert retry['original_failed_samples_concatenated'] is False and retry['original_attempts_overwritten'] is False and retry['scientific_eligibility'] is False
    assert json.loads((recovery_root / 'stage_plan.json').read_text()) == dict(plan_sha256=sha(plan['recovery_plan']))
    bind(bindings, recovery_root / 'stage_plan.json')
    jobs = json.loads(Path(source_plan['jobs']).read_text()); index = {j['chain']['chain_id']: j for j in jobs}
    rows, quartets, summary = build_overlay(jobs, original['chains'], retry['recovery_chains'], recovery_plan['failed_chain_ids'])
    assert retry['full_original_grid_chains'] == summary['full_native_chains'] and retry['original_checked_chains'] == summary['original_checked_chains']
    assert retry['recovery_status_counts'] == dict(Counter(r['status'] for r in retry['recovery_chains']))
    attempts = []
    for record in retry['recovery_chains']:
        attempts.append(audit_attempt(index[record['chain_id']], record, original_root, recovery_root,
            source_plan['iterations'], source_plan['mapping'], recovery_plan['native_address_space_bytes'], bindings))
        print('recovery_attempt_checked', record['chain_id'], record['status'], flush=True)
    # All original checked sample-audit metadata stay linked to their original
    # closed attempts. Their native alignments are hash-checked, not reparsed.
    for row in rows:
        selected = row['selected_disposition']; bind(bindings, selected['receipt'], selected['receipt_sha256'])
        if selected['status'] == CHECKED:
            bind(bindings, selected['sample_audit'], selected['sample_audit_sha256'])
            audit = json.loads(Path(selected['sample_audit']).read_text())
            assert audit['status'] == CHECKED and audit['attempt_receipt_sha256'] == selected['receipt_sha256']
            assert Path(audit['attempt_receipt']).resolve() == Path(selected['receipt']).resolve()
            assert audit['iterations'] == 1000 and audit['mapping_sha256'] == sha(source_plan['mapping'])
            assert len(audit['candidate_samples']) == 404
    verify(bindings)
    overlay = root / 'full_chain_overlay.json'
    with overlay.open('x') as f:
        json.dump(dict(status='verified_full_original_grid_with_whole_recovery_attempt_overlay',
            rows=rows, quartets=quartets, summary=summary,
            original_samples_concatenated=False, scientific_eligibility=False), f, indent=2); f.write('\n')
    receipt = dict(status='passed_full_baliphy_recovery_integrity_and_original_grid_readback',
        plan_sha256=sha(path), recovery_plan_sha256=sha(plan['recovery_plan']),
        recovery_receipt=str(retry_path), recovery_receipt_sha256=sha(retry_path),
        **summary, attempts=attempts, original_terminal_journals=journals,
        overlay=str(overlay), overlay_sha256=sha(overlay), source_hashes=bindings,
        elapsed_seconds=time.monotonic()-started, scientific_eligibility=False,
        scope=plan['scope'])
    with (root / 'receipt.json').open('x') as f: json.dump(receipt, f, indent=2); f.write('\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in {'source_hashes','attempts','original_terminal_journals'}}, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True)
    run(parser.parse_args().plan)
