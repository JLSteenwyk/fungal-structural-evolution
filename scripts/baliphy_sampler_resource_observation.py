"""Read-only identity-bound native memory observations; no peak guarantee."""
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

import psutil

from ancestral_chain_attempt import sha


MEMORY_FIELDS = ['VmPeak', 'VmSize', 'VmHWM', 'VmRSS', 'VmSwap']
SUMMARY_FIELDS = ['expected_roles', 'roles_with_native_attempt', 'roles_with_live_observation',
    'roles_without_live_observation', 'roles_without_native_attempt', 'native_outcome_status_counts',
    'native_observations', 'unavailable_observations', 'sampler_controller_terminal_status',
    'maximum_observed_cgroup_memory_bytes', 'maximum_observed_cgroup_reported_peak_bytes',
    'maximum_observed_poll_gap_seconds', 'posterior_qualified']


def parse_status(text):
    result = {}
    for line in text.splitlines():
        name = line.split(':', 1)[0]
        if name not in MEMORY_FIELDS: continue
        assert name not in result, 'Duplicate process memory field'
        match = re.fullmatch(re.escape(name) + r':\s+(\d+)\s+kB\s*', line)
        assert match, ('Malformed process memory field', name)
        result[name] = int(match[1]) * 1024
    assert set(result) == set(MEMORY_FIELDS), 'Incomplete process memory observation'
    return result


def parse_limits(text, job):
    result = {}
    expected = [('Max address space', job['memory_reservation_bytes'], 'bytes'),
                ('Max cpu time', job['config']['timeout_seconds'], 'seconds'),
                ('Max file size', 2 * 2**30, 'bytes')]
    for name, value, unit in expected:
        matches = re.findall(r'^' + re.escape(name) + r'\s+(\d+)\s+(\d+)\s+(\S+)\s*$', text, re.M)
        assert len(matches) == 1
        soft, hard, actual_unit = matches[0]
        assert int(soft) == int(hard) == value and actual_unit == unit
        result[name] = value
    return result


def describe(job, identity_path):
    identity_path = Path(identity_path); identity = json.loads(identity_path.read_text())
    assert identity['command'] == job['config']['command'] and identity['pgid'] == identity['pid']
    config = identity_path.parent.parent / 'configuration.json'
    assert json.loads(config.read_text()) == job['config']
    command = identity_path.parent / 'command.json'
    assert json.loads(command.read_text()) == identity['command']
    return dict(chain_id=job['chain']['chain_id'], pid=identity['pid'], created=identity['created'],
        command=identity['command'], identity_path=str(identity_path), identity_sha256=sha(identity_path),
        configuration_path=str(config), configuration_sha256=sha(config),
        command_path=str(command), command_sha256=sha(command))


def probe_native(job, descriptor, expected_cgroup, proc_root=Path('/proc')):
    """Check the original process twice around read-only procfs snapshots."""
    base = dict(chain_id=descriptor['chain_id'], pid=descriptor['pid'], created=descriptor['created'],
        identity_sha256=descriptor['identity_sha256'], checked_utc=datetime.now(timezone.utc).isoformat())
    try:
        process = psutil.Process(descriptor['pid'])
        if process.create_time() != descriptor['created'] or process.status() == psutil.STATUS_ZOMBIE:
            return dict(**base, status='original_process_unavailable_at_observation')
        command = job['config']['command']; observed = process.cmdline()
        allowed = [command, command[command.index('--') + 1:]]
        if observed not in allowed:
            # cmdline can become empty while the same process exits or execs.
            # No memory sample is accepted during uncertain identity. A stable
            # unrelated nonempty command still fails the identity contract.
            confirm = psutil.Process(descriptor['pid'])
            if confirm.create_time() != descriptor['created'] or confirm.status() == psutil.STATUS_ZOMBIE:
                return dict(**base, status='original_process_unavailable_at_observation')
            confirmed = confirm.cmdline()
            if not observed or not confirmed or confirmed in allowed:
                return dict(**base, status='native_command_unavailable_or_transitioning_at_observation')
            raise AssertionError('Stable native command mismatch: ' + repr(observed) + ' -> ' + repr(confirmed))
        if observed == command:
            # process.json can appear before prlimit applies limits and execs
            # the native binary. Do not mistake wrapper memory for native use.
            return dict(**base, status='limit_wrapper_before_native_exec_at_observation')
        folder = proc_root / str(descriptor['pid'])
        cgroup = (folder / 'cgroup').read_text()
        assert [line[3:] for line in cgroup.splitlines() if line.startswith('0::')] == [expected_cgroup]
        status = (folder / 'status').read_text(); limits = (folder / 'limits').read_text()
        memory = parse_status(status); caps = parse_limits(limits, job)
        cpu_seconds = sum(process.cpu_times()[:2])
        after = psutil.Process(descriptor['pid'])
        if after.create_time() != descriptor['created'] or after.status() == psutil.STATUS_ZOMBIE:
            return dict(**base, status='original_process_unavailable_after_observation')
        later = after.cmdline()
        if not later or later == command:
            return dict(**base, status='native_command_unavailable_or_transitioning_at_observation')
        assert later == allowed[1], 'Native command changed after resource read: ' + repr(later)
        return dict(**base, status='verified_live_native_resource_observation', command=observed,
            cgroup=expected_cgroup, actual_limits=caps, reported_memory_bytes=memory, cpu_seconds=cpu_seconds,
            raw_proc_status_sha256=hashlib.sha256(status.encode()).hexdigest(),
            raw_proc_limits_sha256=hashlib.sha256(limits.encode()).hexdigest(),
            scope='Read-only Linux-reported virtual/resident values. RSS accounting is approximate; a sampled high-water value can miss a later peak. Not a measured final per-attempt peak or long-chain memory bound.')
    except (psutil.NoSuchProcess, psutil.ZombieProcess, ProcessLookupError, FileNotFoundError) as error:
        return dict(**base, status='process_unavailable_during_observation', error_type=type(error).__name__)


def read_cgroup(record):
    group = next(x[3:] for x in (Path('/proc') / str(record['pid']) / 'cgroup').read_text().splitlines() if x.startswith('0::'))
    root = Path('/sys/fs/cgroup') / group.lstrip('/')
    limits = {k: (root / k).read_text().strip() for k in ['cpu.max', 'memory.max', 'memory.swap.max']}
    assert limits == record['actual_cgroup_limits']
    counters = {k: int((root / k).read_text().strip()) for k in ['memory.current', 'memory.peak', 'memory.swap.current']}
    events = dict(line.split() for line in (root / 'memory.events').read_text().splitlines())
    assert {'low', 'high', 'max', 'oom', 'oom_kill'} <= set(events)
    events = {k: int(v) for k, v in events.items()}; assert all(v >= 0 for v in events.values())
    return dict(path=group, limits=limits, memory_bytes=counters, memory_events=events,
        scope='Cgroup counters include descendants and kernel/cache charges; not the sum of native worker RSS. No counter resets or configuration writes.')


def accumulate(accumulator, sample):
    cid = sample['chain_id']; row = accumulator[cid]
    assert sample['pid'] == row['pid'] and sample['created'] == row['created']
    assert sample['identity_sha256'] == row['identity_sha256']
    if sample['status'] != 'verified_live_native_resource_observation':
        assert sample['status'] in ['original_process_unavailable_at_observation',
            'original_process_unavailable_after_observation', 'process_unavailable_during_observation',
            'limit_wrapper_before_native_exec_at_observation',
            'native_command_unavailable_or_transitioning_at_observation']
        row['unavailable_observations'] += 1; return
    row['live_observations'] += 1
    if row['first_observed_utc'] is None: row['first_observed_utc'] = sample['checked_utc']
    row['last_observed_utc'] = sample['checked_utc']
    for name, number in sample['reported_memory_bytes'].items():
        assert name in MEMORY_FIELDS and isinstance(number, int) and number >= 0
        previous = row['maximum_reported_memory_bytes'][name]
        row['maximum_reported_memory_bytes'][name] = number if previous is None else max(previous, number)
    assert sample['cpu_seconds'] >= 0
    row['maximum_observed_cpu_seconds'] = max(row['maximum_observed_cpu_seconds'] or 0, sample['cpu_seconds'])


def empty_role(job, descriptor=None):
    c = job['chain']; role = dict(chain_id=c['chain_id'], family=c['family'], prior_label=c['prior_label'],
        chain_role=c['chain'], effective_input_group=c['effective_input_group'], seed=c['seed'],
        source_seed=job['source_seed'], original_configuration_ids=c['original_configuration_ids'],
        memory_reservation_bytes=job['memory_reservation_bytes'], pid=None, created=None, identity_sha256=None,
        live_observations=0, unavailable_observations=0, first_observed_utc=None, last_observed_utc=None,
        maximum_reported_memory_bytes={k: None for k in MEMORY_FIELDS}, maximum_observed_cpu_seconds=None,
        native_outcome_status='no_native_attempt_before_controller_termination', exit_code=None,
        native_receipt=None, native_receipt_sha256=None, scientific_eligibility=False, posterior_qualified=False)
    if descriptor is not None:
        assert descriptor['chain_id'] == c['chain_id']
        role.update({k: descriptor[k] for k in ['pid', 'created', 'identity_sha256']})
    return role


def finish_roles(jobs, descriptors, samples, receipts):
    rows = {j['chain']['chain_id']: empty_role(j, descriptors.get(j['chain']['chain_id'])) for j in jobs}
    assert len(rows) == len(jobs)
    for sample in samples: accumulate(rows, sample)
    for cid, entry in receipts.items():
        row = rows[cid]; descriptor = descriptors[cid]; path = Path(entry['path'])
        receipt = json.loads(path.read_text()); assert sha(path) == entry['sha256']
        assert receipt['artifacts']['process.json'] == descriptor['identity_sha256']
        assert receipt['artifacts']['command.json'] == descriptor['command_sha256']
        assert sha(descriptor['identity_path']) == descriptor['identity_sha256']
        assert sha(descriptor['configuration_path']) == descriptor['configuration_sha256']
        assert sha(descriptor['command_path']) == descriptor['command_sha256']
        expected = hashlib.sha256(json.dumps(next(j['config'] for j in jobs if j['chain']['chain_id'] == cid),
            sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
        assert receipt['configuration_sha256'] == expected
        row.update(native_outcome_status=receipt['status'], exit_code=receipt['exit_code'],
            native_receipt=str(path), native_receipt_sha256=entry['sha256'])
    for cid in set(descriptors) - set(receipts): rows[cid]['native_outcome_status'] = 'native_attempt_without_receipt_at_controller_termination'
    return sorted(rows.values(), key=lambda x: x['chain_id'])


def summarize(rows, snapshots, terminal_status):
    assert len(rows) == len({r['chain_id'] for r in rows}) == 1620
    maximum_current = maximum_peak = maximum_gap = None
    for snapshot in snapshots:
        current = snapshot['cgroup']['memory_bytes']['memory.current']
        peak = snapshot['cgroup']['memory_bytes']['memory.peak']
        gap = snapshot['elapsed_since_previous_snapshot_seconds']
        maximum_current = current if maximum_current is None else max(maximum_current, current)
        maximum_peak = peak if maximum_peak is None else max(maximum_peak, peak)
        maximum_gap = gap if maximum_gap is None else max(maximum_gap, gap)
    return dict(expected_roles=len(rows), roles_with_native_attempt=sum(r['pid'] is not None for r in rows),
        roles_with_live_observation=sum(r['live_observations'] > 0 for r in rows),
        roles_without_live_observation=sum(r['live_observations'] == 0 for r in rows),
        roles_without_native_attempt=sum(r['pid'] is None for r in rows),
        native_outcome_status_counts=dict(Counter(r['native_outcome_status'] for r in rows)),
        native_observations=sum(r['live_observations'] for r in rows),
        unavailable_observations=sum(r['unavailable_observations'] for r in rows),
        sampler_controller_terminal_status=terminal_status,
        maximum_observed_cgroup_memory_bytes=maximum_current,
        maximum_observed_cgroup_reported_peak_bytes=maximum_peak,
        maximum_observed_poll_gap_seconds=maximum_gap,
        posterior_qualified=False)
