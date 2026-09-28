"""Freeze all selected-fit/cache bindings for candidate interval evaluation.

This prepares the full grid; it neither runs intervals nor qualifies coverage.
The full cache replay must finish successfully before the evaluation launch.
"""
import argparse
import json
from pathlib import Path
import shutil
import pandas as pd
from ancestral_chain_attempt import sha, write_json


def read_receipt(root, status):
    receipt = json.loads((root / 'receipt.json').read_text())
    assert receipt['status'] == status
    for name, digest in receipt['artifacts'].items():
        assert sha(root / name) == digest, name
    return receipt


def indexed_manifest(path):
    result = {}
    for line in path.open():
        row = json.loads(line)
        key = (row['fit_input_id'], row['tree'])
        assert key not in result, key
        result[key] = row
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    cache = Path('results/model_validation/matched-simulation-input-cache-20260928-v1')
    selected = Path('results/structural_comparisons/refined-working-model-grid-export-20260928-v2')
    original = Path('results/structural_comparisons/full-matched-working-models-20260927-v1')
    refined = Path('results/structural_comparisons/matched-reml-analytic-refinement-20260928-v1')
    sr = read_receipt(selected, 'complete_refined_estimate_overlay')
    cr = read_receipt(cache, 'all_exact_simulation_inputs_cached_pending_likelihood_replay')
    for name, digest in sr['source_bindings'].items():
        assert sha(name) == digest, name
    old_receipt = json.loads((original / 'receipt.json').read_text())
    assert sha(original / 'fit_manifest.jsonl') == old_receipt['artifacts']['fit_manifest.jsonl']
    originals = indexed_manifest(original / 'fit_manifest.jsonl')
    refinements = indexed_manifest(refined / 'manifest.jsonl')
    frame = pd.read_parquet(selected / 'unique_fits.parquet')
    settings = pd.read_parquet(selected / 'full_settings.parquet')
    keys = ['fit_input_id', 'tree']
    assert len(frame) == 144040 and not frame.duplicated(keys).any()
    assert len(settings) == 414720
    assert set(map(tuple, settings[keys].to_numpy())) == set(map(tuple, frame[keys].to_numpy()))
    # Check the expansion preserves every selected estimate and disposition.
    columns = [c for c in frame if c.startswith('selected_')] + ['selection']
    joined = settings.merge(frame[keys + columns], on=keys, validate='many_to_one', suffixes=('_setting', '_unique'))
    for col in columns:
        a, b = joined[col + '_setting'], joined[col + '_unique']
        assert ((a == b) | (a.isna() & b.isna())).all(), col
    factors = {p.stem: {'path': str(p), 'sha256': sha(p)} for p in sorted((cache / 'factors').glob('*.npz'))}
    assert len(factors) == 5 and set(frame.tree) == set(factors)
    entries = {}
    for line in (cache / 'manifest.jsonl').open():
        entry = json.loads(line)
        assert entry['fit_input_id'] not in entries
        entries[entry['fit_input_id']] = entry
    assert len(entries) == 28808 and set(frame.fit_input_id) == set(entries)
    tasks = []
    for identifier, group in frame.groupby('fit_input_id', sort=True):
        assert len(group) == 5 and set(group.tree) == set(factors)
        entry = entries[identifier]
        assert set(group.records) == {entry['records']}
        sources = {}
        for row in group.to_dict('records'):
            key = (identifier, row['tree'])
            assert originals[key]['sha256'] == row['source_fit_sha256']
            if row['selection'] == 'refined_candidate':
                source = refinements[key]
            else:
                assert row['selection'] == 'original_not_targeted'
                source = originals[key]
            assert source['sha256'] == row['selected_source_sha256']
            sources[row['tree']] = {'path': source['path'], 'sha256': source['sha256']}
        tasks.append({'fit_input_id': identifier, 'cache_entry': entry, 'selected_sources': sources})
    assert sum(len(t['selected_sources']) for t in tasks) == 144040
    assert shutil.disk_usage('.').free > 100 * 2**30
    args.output.mkdir(parents=True, exist_ok=False)
    manifest = args.output / 'tasks.jsonl'
    with manifest.open('w') as handle:
        for task in tasks:
            handle.write(json.dumps(task, sort_keys=True, allow_nan=False) + '\n')
    bindings = {str(p): sha(p) for p in [selected/'receipt.json', selected/'unique_fits.parquet',
        selected/'full_settings.parquet', cache/'receipt.json', cache/'manifest.jsonl',
        original/'receipt.json', original/'fit_manifest.jsonl', refined/'manifest.jsonl']}
    receipt = dict(status='prepared_complete_selected_interval_task_grid_pending_cache_replay_and_runner',
        inputs=len(tasks), unique_fits=len(frame), settings=len(settings), factors=factors,
        selected_review_fits=int(frame.selected_review_required.sum()),
        selections=frame.selection.value_counts().to_dict(), source_bindings=bindings,
        artifacts={'tasks.jsonl': sha(manifest)}, script_sha256=sha(__file__),
        required_prerequisite=dict(launch='metadata/matched_simulation_cache_replay_launch_20260928.json',
            receipt='results/model_validation/matched-simulation-cache-replay-20260928-v1/receipt.json',
            status='all_144040_original_fits_replayed_from_simulation_cache',
            terminal_state='inactive', result='success', exit_status=0),
        proposed_resources=dict(cpus=16, workers=16, maximum_inflight_tasks=32, memory_gib=32,
            swap_gib=0, output_allowance_gib=64, free_disk_reserve_gib=100,
            planning_active_hours=[1,120], paid_cost=0, gpu=False,
            basis='15 existing real-design timings: approximately .02–15.6 seconds per fit; heterogeneous sizes. '
                  '144040 fits all at 15.6 seconds would take 39 ideal hours on16 workers. '
                  'Planning range is not a guaranteed bound; record measured throughput during full run.'),
        scope='All selected estimates and full-setting mappings bound. Source payload and cache content hashes '
              'must be rechecked by the runner, with all numerical failures retained as explicit review outcomes. '
              'No intervals computed here. Marginal candidate intervals require coverage, model adequacy and multiplicity work.')
    write_json(args.output/'receipt.json', receipt)
    write_json(Path('metadata/selected_matched_interval_grid_preparation_20260928.json'), receipt)
    print(json.dumps({k:receipt[k] for k in ['status','inputs','unique_fits','settings','selections']}))


if __name__ == '__main__':
    main()
