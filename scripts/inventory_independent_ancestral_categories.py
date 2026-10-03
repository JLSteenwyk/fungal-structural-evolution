#!/usr/bin/env python3
"""Inventory the complete closed categorical grid before a numerical replay.

Rechecks small source JSON and the archive hash. Large NPZ/JSONL data hashes
are recorded from the closed archive, not rehashed or numerically replayed here.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path

from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def run(completion, output):
    bindings = {str(Path(__file__)): sha(__file__)}
    closed = json.loads(completion.read_text()); bind(bindings, completion)
    assert closed['status'] == 'complete_verified_full_baliphy_recovery_diagnostics'
    assert closed['full_quartets'] == 405 and closed['complete_quartets'] == 403
    assert closed['unresolved_quartets'] == 2 and closed['scientific_eligibility'] is False
    archive_path = Path(closed['full_hash_archive'])
    bind(bindings, archive_path, closed['full_hash_archive_sha256'])
    assert sha(archive_path) == closed['full_hash_archive_sha256']
    archive = json.loads(archive_path.read_text())
    assert archive['status'] == 'complete_verified_full_baliphy_recovery_diagnostics_archive'
    assert len(archive['services']) == closed['exact_process_journals_checked'] == 2
    assert len(archive['source_hashes']) == closed['bound_source_hashes'] == 50019
    lookup = {str(Path(p).resolve()): d for p, d in archive['source_hashes'].items()}
    data_pins = {}

    def document(path, digest=None):
        path = Path(path); expected = lookup[str(path.resolve())]
        assert digest is None or digest == expected
        bind(bindings, path, expected); assert sha(path) == expected
        return json.loads(path.read_text())

    def data(path, digest):
        path = Path(path); assert lookup[str(path.resolve())] == digest
        assert str(path) not in data_pins or data_pins[str(path)] == digest
        data_pins[str(path)] = digest
        return dict(path=str(path), sha256=digest, bytes_at_observation=path.stat().st_size)

    producer = document(closed['producer_receipt'], closed['producer_receipt_sha256'])
    assert len(producer['groups']) == 405
    groups = {}; seen = set(); totals = Counter(); patterns = Counter(); coordinates = Counter()
    status = {key:Counter() for key in ['pattern_status_counts', 'coordinate_status_counts']}
    largest = dict(full_values_bytes=0, raw_trace_coordinate_count=0, patterns_per_cutoff=0)
    for group, info in sorted(producer['groups'].items()):
        ids = info['chain_ids']
        assert len(ids) == len(set(ids)) == 4 and not seen.intersection(ids)
        seen.update(ids)
        entry = dict(chain_ids=ids, source_status=info['status'], scientific_eligibility=False)
        if info['status'] == 'unresolved_failed_native_chain_retained':
            assert 'categorical' not in info
            totals['unresolved_quartets'] += 1; groups[group] = entry; continue
        assert info['status'] == 'complete_scalar_length_category_screens_not_posterior_qualification'
        category = document(info['categorical']['receipt'], info['categorical']['receipt_sha256'])
        assert category['group'] == group
        assert category['status'] == 'verified_quartet_categorical_reports_complete_not_posterior_qualification'
        manifest = document(category['manifest'], category['manifest_sha256'])
        report = document(category['report'], category['report_sha256'])
        assert manifest['status'] == 'provenance_checked_state_quartet'
        assert report['status'] == 'both_cutoff_categorical_reports_complete_not_posterior_qualification'
        assert report['arviz_version'] == '0.22.0' and report['numpy_version'] == '2.2.6'
        assert {c['chain_id'] for c in manifest['chains']} == set(ids)
        assert len(manifest['chains']) == len({c['seed'] for c in manifest['chains']}) == 4
        assert all(c['model_input_identity'] == group for c in manifest['chains'])
        assert manifest['expected_iterations'] == list(range(0, 1001, 10))
        coords = manifest['coordinates']; alphabet = list(coords['alphabet'])
        assert alphabet == list('ACDEFGHIKLMNPQRSTVWYX-')
        assert len(coords['nodes']) == len(set(coords['nodes'])) == 4
        count = 4 * sum(t['length'] for t in coords['tips'])
        full_bytes = 4 * 101 * count
        largest['full_values_bytes'] = max(largest['full_values_bytes'], full_bytes)
        largest['raw_trace_coordinate_count'] = max(largest['raw_trace_coordinate_count'], count)
        totals['full_values_bytes'] += full_bytes
        entry.update(category_receipt=info['categorical'], manifest=category['manifest'],
            manifest_sha256=category['manifest_sha256'], report=category['report'],
            report_sha256=category['report_sha256'], alphabet=alphabet,
            full_array=data(manifest['arrays'], manifest['arrays_sha256']), cutoffs={})
        assert set(report['outputs']) == {'250', '500'}
        for cutoff in ['250', '500']:
            folder = Path(category['report']).parent / ('discard-' + cutoff)
            summary = document(folder / 'summary.json')
            expected = report['outputs'][cutoff]
            assert summary['discard_through'] == int(cutoff)
            assert summary['coordinates'] == expected['coordinates'] == count
            assert summary['patterns'] == expected['patterns']
            assert summary['retained_samples_per_chain'] == expected['retained_samples_per_chain'] == (75 if cutoff == '250' else 50)
            assert sum(summary['coordinate_status_counts'].values()) == count
            assert sum(summary['pattern_status_counts'].values()) == summary['patterns']
            patterns[cutoff] += summary['patterns']; coordinates[cutoff] += count
            largest['patterns_per_cutoff'] = max(largest['patterns_per_cutoff'], summary['patterns'])
            for key in status:
                status[key].update({cutoff + ':' + k:v for k,v in summary[key].items()})
            files = {}
            for name in ['patterns.npz', 'diagnostics.jsonl.gz']:
                path = folder / name; relative = str(path.relative_to(Path(category['report']).parent))
                files[name] = data(path, report['artifacts'][relative])
                totals[name + '_compressed_bytes'] += files[name]['bytes_at_observation']
            entry['cutoffs'][cutoff] = dict(**expected, files=files)
        totals['complete_quartets'] += 1; groups[group] = entry
    assert len(seen) == 1620 and totals['complete_quartets'] == 403 and totals['unresolved_quartets'] == 2
    assert dict(patterns) == closed['categorical_pattern_counts']
    assert dict(coordinates) == closed['categorical_coordinate_counts']
    assert all(dict(counts) == closed[key] for key,counts in status.items())
    verify(bindings)
    result = dict(status='complete_full_categorical_replay_inventory_not_numerical_validation',
        checked_utc=datetime.now(timezone.utc).isoformat(), full_quartets=405,
        complete_quartets=403, unresolved_quartets=2, original_chain_ids=1620,
        original_cutoffs=[250, 500], alphabet=list('ACDEFGHIKLMNPQRSTVWYX-'),
        pattern_counts=dict(patterns), coordinate_counts=dict(coordinates),
        total_pattern_cutoff_rows=sum(patterns.values()),
        declared_indicator_rows=22 * sum(patterns.values()),
        source_status_counts={k:dict(v) for k,v in status.items()}, totals=dict(totals), largest=largest,
        groups=groups, small_sources_rehashed=bindings, large_data_pins_not_rehashed=data_pins,
        proposed_resources=dict(cpu_equivalents=2, memory_gib=32, swap_bytes=0,
            blas_threads=1, output_allowance_gib=128, minimum_free_disk_gib=228,
            estimated_array_workspace_bytes=16 * largest['full_values_bytes'],
            array_workspace_estimate_is_hard_bound=False,
            uncalibrated_wall_hours_per_stage=[4, 96], production_runtime_measured=False,
            finish_eta=None, fitting_or_sampling_launched=False, gpu=False, new_cost_usd=0),
        scientific_eligibility=False, native_parser_independently_verified=False,
        scope='All405originalquartets/bothcutoffs/22declaredstates, whole selected attempts '
              'and failed groups retained. Small source metadata and archive SHA rechecked. '
              'Compressed-data hashes come from the closed50019-binding archive and have not '
              'been rehashed/decompressed/numerically replayed here. Resource planning only; '
              'full data hashes, complete pattern/coordinate reconstruction and every indicator '
              'numeric comparison are required before any production numerical-validation claim.')
    with output.open('x') as handle: handle.write(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['groups', 'small_sources_rehashed', 'large_data_pins_not_rehashed']}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--completion', type=Path, default=Path('metadata/baliphy_recovery_full_diagnostics_completed_20261002.json'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); run(args.completion, args.output)
