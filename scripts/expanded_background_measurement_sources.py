"""Closed input and complete native checkpoint contracts for new background work."""
import csv
import json
import math
from pathlib import Path
from expanded_background_native_handoff import load_handoff
from reference_measurement_union_sources import bind, verify
from run_cross_clan_alignments import parse_output
from run_ortholog_pair_guide_comparison import sha


def load_native(plan, include_positions=True):
    source_path = Path(plan['native_plan']); source = json.loads(source_path.read_text())
    inputs, pairs, bindings, bundle, partition, input_root = load_handoff(source, include_positions)
    bind(bindings, source_path)
    root = Path(source['output']); rp = root / 'receipt.json'; receipt = json.loads(rp.read_text())
    assert receipt['status'] == 'complete_expanded_background_alignment_dispositions_pending_numeric_geometry'
    assert receipt['plan_sha256'] == sha(source_path) and receipt['input_bundle_sha256'] == bundle
    assert receipt['upstream_bindings'] == bindings
    assert receipt['scientific_eligibility'] is False
    assert receipt['full_background_pairs'] == len(partition) == source['expected']['full_pairs']
    assert receipt['distinct_model_pairs'] == len(pairs) == source['expected']['new_pairs']
    assert receipt['directed_dispositions'] == 4 * len(pairs)
    assert receipt['existing_catalog_pairs_pending_reuse'] == len(partition) - len(pairs)
    bind(bindings, rp)
    for name, digest in receipt['artifacts'].items(): bind(bindings, root / name, digest)
    assert (root / 'full_background_work_partition.tsv').read_bytes() == (input_root / 'full_background_work_partition.tsv').read_bytes()
    new = {row['pair_key']: [(row['model_a'], int(row['version_a'])), (row['model_b'], int(row['version_b']))] for row in pairs}
    expected = {(pair, mask, order) for pair in new for mask in ['full', 'plddt70'] for order in [0, 1]}
    checkpoints = {}
    with (root / 'checkpoint_manifest.tsv').open() as handle:
        for item in csv.DictReader(handle, delimiter='\t'):
            pair, mask, order = Path(item['path']).stem.rsplit('-', 2); key = pair, mask, int(order)
            assert key in expected and key not in checkpoints
            assert item['path'] == f'pairs/{pair[:2]}/{pair}-{mask}-{order}.json'
            checkpoints[key] = item; bind(bindings, root / item['path'], item['sha256'])
    assert set(checkpoints) == expected
    for path, digest in plan['pins'].items(): bind(bindings, path, digest)
    verify(bindings)
    return source, root, receipt, inputs, new, checkpoints, bindings, bundle


def inspect_native(source, root, item, key, ends, rows, bundle):
    """Check immutable byte/provenance bindings and each retained disposition."""
    pair, mask, order = key; path = root / item['path']; assert sha(path) == item['sha256']
    record = json.loads(path.read_text()); ready = all(row['status'] == 'ready' for row in rows)
    assert (record['pair_key'], record['mask'], record['order']) == key
    assert record['status'] == item['status'] and record['input_manifest_sha256'] == bundle
    identities = [dict(model_id=end[0], version=end[1], status=row['status'], path=row.get('path'), sha256=row.get('sha256')) for end, row in zip(ends, rows)]
    assert record['inputs'] == identities
    command = [source['usalign'], *[row['path'] for row in rows], *source['options']] if ready else None
    assert record['command'] == command
    assert math.isfinite(record['elapsed_seconds']) and record['elapsed_seconds'] >= 0
    status = record['status']
    if status == 'aligned':
        assert ready and record['returncode'] == 0
        assert parse_output(record['stdout'], [row['sequence'] for row in rows]) == record['metrics']
    elif status == 'input_unavailable':
        assert not ready and record['input_statuses'] == [row['status'] for row in rows]
        assert all(name not in record for name in ['metrics', 'stdout', 'returncode'])
    elif status == 'native_error': assert ready and record['returncode'] != 0 and 'metrics' not in record
    elif status == 'parse_error':
        assert ready and record['returncode'] == 0 and record.get('error') and 'metrics' not in record
        try: parse_output(record['stdout'], [row['sequence'] for row in rows])
        except (ValueError, AttributeError, StopIteration, IndexError): pass
        else: raise AssertionError('Parse-error checkpoint is actually parseable')
    elif status == 'timeout': assert ready and record['timeout_seconds'] == source['per_pair_timeout_seconds'] and 'metrics' not in record
    else: raise AssertionError('Unknown native disposition: ' + status)
    return record
