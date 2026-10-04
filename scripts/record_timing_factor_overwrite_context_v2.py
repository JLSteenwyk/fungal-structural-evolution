#!/usr/bin/env python3
"""Retain byte-span evidence and the exact live diagnostic's native library maps."""
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import psutil

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    output = Path('metadata/retained_timing_factor_overwrite_context_20261004_v2.json')
    assert not output.exists()
    prior_path = Path('metadata/retained_timing_ordered_probe_diagnostic_20261004_v2.json')
    prior = json.loads(prior_path.read_text())
    capture = prior['detected_input_mutation']['factor_value_change_capture']['pmsf_profile_profile']
    baseline_path = Path('results/phylogeny/retained-timing-covariance-diagnostic-20261004-v1/pmsf_profile_profile-factor.npy')
    snapshot_path = Path(capture['observed_snapshot'])
    assert sha(snapshot_path) == capture['observed_snapshot_sha256']
    baseline = np.load(baseline_path, mmap_mode='r', allow_pickle=False)
    observed = np.load(snapshot_path, mmap_mode='r', allow_pickle=False)
    assert baseline.shape == observed.shape == (22881, 301)
    assert baseline.dtype == observed.dtype == np.dtype('float64')
    old = baseline.view(np.uint8).reshape(-1)
    new = observed.view(np.uint8).reshape(-1)
    offsets = np.flatnonzero(old != new)
    assert len(offsets) == 56
    assert offsets.tolist() == list(range(20500226, 20500282))
    assert np.all(new[offsets] == 0)
    layout_path = Path('results/phylogeny/retained-timing-ordered-probe-diagnostic-20261004-v2/factor_buffer_layout.json')
    layout = json.loads(layout_path.read_text())
    address = layout['pmsf_profile_profile']['pointer'] + int(offsets[0])
    config_path = Path('metadata/retained_timing_ordered_probe_diagnostic_execution_20261004_v6/configuration.json')
    process_path = config_path.parent / 'process.json'
    config = json.loads(config_path.read_text())
    identity = json.loads(process_path.read_text())
    native = psutil.Process(identity['pid'])
    assert native.create_time() == identity['created']
    assert native.cmdline() == identity['command'][identity['command'].index('--')+1:]
    assert native.status() != psutil.STATUS_ZOMBIE
    maps = Path(f'/proc/{native.pid}/maps').read_text()
    libraries = sorted({line.split(None, 5)[-1] for line in maps.splitlines()
                        if len(line.split(None, 5)) == 6
                        and '.so' in line.split(None, 5)[-1]
                        and line.split(None, 5)[-1].startswith('/')})
    assert libraries and all(Path(name).is_file() for name in libraries)
    pins = {}
    for path in [Path(__file__), prior_path, baseline_path, snapshot_path, layout_path,
                 config_path, process_path, *map(Path, libraries)]:
        bind(pins, path)
    verify(pins)
    assert native.create_time() == identity['created'] and native.status() != psutil.STATUS_ZOMBIE
    result = dict(status='verified_saved_factor_byte_span_and_live_native_library_context',
        checked_utc=datetime.now(timezone.utc).isoformat(),
        changed_factor_entries=8, changed_bytes=56, changed_bytes_contiguous=True,
        changed_byte_offsets_in_array=[20500226, 20500281],
        changed_bytes_all_observed_zero=True, first_changed_address_modulo_64=address % 64,
        live_original_tool_session_id=82836, live_invocation_id=config['invocation_id'],
        live_native_pid=native.pid, live_native_created=native.create_time(),
        live_native_command=native.cmdline(), library_count=len(libraries),
        mapped_library_paths=libraries, source_hashes=pins,
        root_cause_established=False, original_jobs_restarted=False,
        scientific_eligibility=False, gpu=False,
        scope='Fresh byte comparison of unchanged baseline and actual retained V2 snapshot; '
              '56 contiguous changed bytes are zero in observed state. Current V6 native identity '
              'and mapped-library file hashes are recorded for further isolation. Maps belong '
              'to V6, not an observation of terminated V2 or the original failed timing process. '
              'No attribution to a library, native instruction or hardware; no repair or fit acceptance.')
    with output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ['source_hashes', 'mapped_library_paths', 'live_native_command']}, indent=2))


if __name__ == '__main__':
    main()
