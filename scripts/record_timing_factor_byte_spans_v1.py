#!/usr/bin/env python3
"""Independently compare retained factor snapshots at the byte level."""
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    output = Path('metadata/retained_timing_factor_byte_spans_20261004_v1.json')
    assert not output.exists()
    pins = {}; records = []
    for version in [2, 4]:
        receipt_path = Path(f'metadata/retained_timing_ordered_probe_diagnostic_20261004_v{version}.json')
        transport_path = Path(f'metadata/retained_timing_ordered_probe_diagnostic_transport_20261004_v{version}.json')
        receipt = json.loads(receipt_path.read_text())
        transport = json.loads(transport_path.read_text())
        assert transport['validation_sha256'] == sha(receipt_path)
        assert transport['actual_tool_terminal_exit_code'] == 0
        verify(transport['source_hashes'])
        for tree, capture in receipt['detected_input_mutation']['factor_value_change_capture'].items():
            baseline_path = Path(f'results/phylogeny/retained-timing-covariance-diagnostic-20261004-v1/{tree}-factor.npy')
            observed_path = Path(capture['observed_snapshot'])
            assert sha(observed_path) == capture['observed_snapshot_sha256']
            baseline = np.load(baseline_path, mmap_mode='r', allow_pickle=False)
            observed = np.load(observed_path, mmap_mode='r', allow_pickle=False)
            assert baseline.dtype == observed.dtype == np.dtype('float64')
            assert baseline.shape == observed.shape == (22881, 301)
            old = baseline.view(np.uint8).reshape(-1)
            new = observed.view(np.uint8).reshape(-1)
            changed = np.flatnonzero(old != new)
            assert changed.tolist() == list(range(20500226, 20500282))
            assert np.all(new[changed] == 0)
            records.append(dict(version=version, tree=tree, probe_tree=capture['current_probe_tree'],
                changed_entries=8, changed_bytes=56, first_array_byte_offset=20500226,
                last_array_byte_offset=20500281, contiguous=True, observed_changed_bytes_all_zero=True,
                snapshot=str(observed_path), snapshot_sha256=sha(observed_path)))
            for path in [baseline_path, observed_path]: bind(pins, path)
        for path in [receipt_path, transport_path]: bind(pins, path)
    bind(pins, Path(__file__)); verify(pins)
    result = dict(status='verified_two_saved_factor_overwrite_byte_spans',
        checked_utc=datetime.now(timezone.utc).isoformat(), records=records, source_hashes=pins,
        original_jobs_restarted=False, root_cause_established=False, gpu=False,
        scope='Fresh independent byte comparisons and snapshot checksums for two closed diagnostics. '
              'Both captured changes occupy exactly the same56bytearrayoffsetspan in different trees. '
              'No inference of responsible function, native instruction, hardware cause or original failure repair.')
    with output.open('x') as handle:
        json.dump(result, handle, indent=2); handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
