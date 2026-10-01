#!/usr/bin/env python3
"""Run every new full-scope background pair under both masks and both orders."""
import argparse
import csv
import fcntl
import itertools
import json
import shutil
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from expanded_background_native_handoff import load_handoff
from reference_measurement_union_sources import verify
from run_duplication_alignments import run_job
from run_ortholog_pair_guide_comparison import sha


def run(plan_path):
    plan_path = Path(plan_path); plan = json.loads(plan_path.read_text()); ph = sha(plan_path)
    initial = {str(plan_path): ph, **plan['pins']}; verify(initial)
    inputs, pairs, bindings, bundle, partition, source_root = load_handoff(plan)
    for path, digest in initial.items(): assert path not in bindings or bindings[path] == digest; bindings[path] = digest
    assert sha(plan['usalign']) == plan['pins'][plan['usalign']]
    out = Path(plan['output']); out.mkdir(parents=True, exist_ok=True); lock = (out / 'run.lock').open('w'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    assert not (out / 'receipt.json').exists(), 'Completed run requires review before restart'
    ledger = out / 'full_background_work_partition.tsv'
    if ledger.exists(): assert ledger.read_bytes() == (source_root / ledger.name).read_bytes()
    else: ledger.write_bytes((source_root / ledger.name).read_bytes())
    def jobs():
        for pair in pairs:
            a, b = [(pair['model_' + side], int(pair['version_' + side])) for side in ['a', 'b']]
            for mask in ['full', 'plddt70']:
                for order, (left, right) in enumerate([(a, b), (b, a)]): yield pair['pair_key'], left, right, mask, order
    stream = iter(jobs()); totals = Counter(); completed = 0
    with ThreadPoolExecutor(max_workers=plan['workers']) as pool, (out / 'checkpoint_manifest.tsv').open('w') as manifest:
        writer = csv.writer(manifest, delimiter='\t', lineterminator='\n'); writer.writerow(['path', 'sha256', 'status'])
        while batch := list(itertools.islice(stream, 64)):
            assert shutil.disk_usage(out).free >= plan['resources']['minimum_free_disk_gib'] * 2 ** 30
            futures = [pool.submit(run_job, job, inputs, plan, ph, bundle) for job in batch]
            for job, future in zip(batch, futures):
                path, status = future.result(); totals[job[3] + ':' + status] += 1; completed += 1
                writer.writerow([str(path.relative_to(out)), sha(path), status])
            manifest.flush(); (out / 'state.json').write_text(json.dumps(dict(stage='aligning_full_new_background_partition', completed=completed, total=4 * len(pairs), counts=dict(totals))) + '\n')
            print('Full new-background directed dispositions', completed, '/', 4 * len(pairs), flush=True)
    assert completed == 4 * plan['expected']['new_pairs']; verify(bindings)
    result = dict(status='complete_expanded_background_alignment_dispositions_pending_numeric_geometry', plan_sha256=ph, input_bundle_sha256=bundle, upstream_bindings=bindings,
                  full_background_pairs=len(partition), distinct_model_pairs=len(pairs), directed_dispositions=completed, existing_catalog_pairs_pending_reuse=len(partition) - len(pairs),
                  counts=dict(totals), artifacts={name: sha(out / name) for name in ['checkpoint_manifest.tsv', 'full_background_work_partition.tsv']}, scientific_eligibility=False,
                  scope=plan['scope'])
    with (out / 'receipt.json').open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'upstream_bindings'}, indent=2), flush=True); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True); run(parser.parse_args().plan)
