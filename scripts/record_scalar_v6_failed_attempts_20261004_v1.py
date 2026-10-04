#!/usr/bin/env python3
"""Read preserved unsuccessful V6 attempts and extract their exact telemetry."""
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    output = Path('metadata/baliphy_scalar_v6_failed_attempts_recheck_20261004_v1.json')
    assert not output.exists()
    root = Path('results/ancestral/full-baliphy-scalar-v6-short-sampler-20261004-v1')
    pins = {}
    rows = []
    identities = {}
    for path in sorted((root/'chains').glob('*.json')):
        row = json.loads(path.read_text())
        if row['status'] != 'unsuccessful_sampler_qualification_attempt_retained':
            continue
        native = Path(row['native_receipt'])
        assert sha(native) == row['native_receipt_sha256']
        receipt = json.loads(native.read_text())
        assert receipt['exit_code'] == row['exit_code'] and receipt['status'] == row['native_status']
        artifacts = {str(native.parent/name): digest for name, digest in receipt['artifacts'].items()}
        verify(artifacts)
        for name, digest in artifacts.items():
            bind(pins, Path(name), digest)
        process_path = native.parent/'process.json'
        identity = json.loads(process_path.read_text())
        identities[row['chain_id']] = identity
        stderr = (native.parent/'stderr.log').read_text()
        stdout = (native.parent/'stdout.log').read_text()
        rows.append(dict(chain_id=row['chain_id'], family=row['family'], prior=row['prior_label'],
            role=row['chain_role'], effective_input_group=row['effective_input_group'],
            configuration_aliases=row['original_configuration_ids'], exit_code=row['exit_code'],
            original_identity=identity, elapsed_worker_seconds=row['elapsed_worker_seconds'],
            stderr_bytes=len(stderr.encode()), stdout_has_mcmc_start='Beginning MCMC computations.' in stdout,
            bad_alloc_reported=row['bad_alloc'], allocation_warning_lines=row['allocation_warning_lines'],
            scalar_integrity_accepted=row['scalar_integrity_accepted'],
            native_receipt=str(native), native_receipt_sha256=sha(native),
            sampled_memory_highwater_bytes={}, matching_telemetry_rows=0))
        assert not row['scalar_integrity_accepted'] and not row['scientific_eligibility']
        for saved in [path, native, process_path]:
            bind(pins, saved)
    assert rows, 'No preserved unsuccessful attempts found'
    selected = {row['chain_id']: row for row in rows}
    telemetry = Path('results/ancestral/full-baliphy-scalar-v6-resource-observation-20261004-v1/observations.jsonl')
    # The observer is live. Hash only a bounded, consumed prefix; never claim
    # that a growing full file has a stable digest or a final maximum.
    byte_limit = telemetry.stat().st_size
    prefix_hash = hashlib.sha256()
    prefix_bytes = 0
    matches = []
    cgroup_event_highwater = Counter()
    with telemetry.open('rb') as handle:
        while prefix_bytes < byte_limit:
            line = handle.readline(byte_limit-prefix_bytes)
            if not line.endswith(b'\n'):
                break
            prefix_hash.update(line)
            prefix_bytes += len(line)
            snapshot = json.loads(line)
            for observation in snapshot['native_observations']:
                chain = observation.get('chain_id')
                if chain not in identities or observation.get('status') != 'verified_live_native_resource_observation':
                    continue
                expected = identities[chain]
                assert observation['pid'] == expected['pid'] and observation['created'] == expected['created']
                command = expected['command']
                assert observation['command'] == command[command.index('--')+1:]
                row = selected[chain]
                row['matching_telemetry_rows'] += 1
                for key, value in observation['reported_memory_bytes'].items():
                    prior = row['sampled_memory_highwater_bytes'].get(key, 0)
                    row['sampled_memory_highwater_bytes'][key] = max(prior, value)
                cgroup_event_highwater |= Counter(snapshot['cgroup']['memory_events'])
                matches.append(dict(sequence=snapshot['sequence'], checked_utc=snapshot['checked_utc'],
                    native_observation=observation, contemporaneous_cgroup=snapshot['cgroup']))
    extracted = Path('metadata/baliphy_scalar_v6_failed_attempts_telemetry_extract_20261004_v1.jsonl')
    with extracted.open('x') as handle:
        for match in matches:
            handle.write(json.dumps(match, sort_keys=True)+'\n')
    for path in [Path(__file__), extracted]:
        bind(pins, path)
    verify(pins)
    result = dict(status='read_only_preserved_v6_failed_attempts_and_exact_telemetry_rechecked',
        checked_utc=datetime.now(timezone.utc).isoformat(), failed_attempts=len(rows),
        family_counts=dict(Counter(r['family'] for r in rows)),
        exit_code_counts=dict(Counter(str(r['exit_code']) for r in rows)),
        prior_counts=dict(Counter(r['prior'] for r in rows)), attempts=rows,
        observed_cgroup_event_highwater=dict(cgroup_event_highwater),
        original_telemetry_prefix=dict(path=str(telemetry), consumed_bytes=prefix_bytes,
            sha256=prefix_hash.hexdigest(), full_live_file_hashed=False),
        telemetry_extract=str(extracted), matching_telemetry_rows=len(matches), source_hashes=pins,
        new_native_runs=0, original_attempts_restarted=False, root_cause_established=False,
        posterior_qualified=False, scientific_eligibility=False, gpu=False,
        scope='Every presently saved unsuccessful role is matched to its original native '
              'receipt, identity and artifacts. Exact live identity telemetry is extracted '
              'from a bounded original prefix. Sampled memory values can miss later peaks; '
              'no final-memory or causal diagnosis follows. No priors, caps, seeds, sources '
              'or failed attempts are modified or retried. All roles remain in accounting.')
    with output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items()
                     if k not in ['source_hashes', 'attempts']}, indent=2))
    print('sampled_peak_vm_gib', max(r['sampled_memory_highwater_bytes'].get('VmPeak', 0) for r in rows)/2**30)
    print('sampled_peak_rss_gib', max(r['sampled_memory_highwater_bytes'].get('VmHWM', 0) for r in rows)/2**30)


if __name__ == '__main__':
    main()
