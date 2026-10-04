#!/usr/bin/env python3
"""Observe exact original timing controllers and their numerical admission gate."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from record_process_covariance_pipeline_checkpoint import observe
from reference_measurement_union_sources import bind, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    pp = Path('metadata/full_weighted_timing_parallel_plan_20261004_v1.json')
    ip = Path('metadata/full_weighted_timing_parallel_launches_20261004_v1.json')
    plan, inventory = [json.loads(p.read_text()) for p in [pp, ip]]
    assert inventory['source_plan_sha256'] == sha(pp)
    pins = dict(plan['pins'])
    verify(pins)
    handles = []
    for index, path in enumerate(inventory['launches']):
        limits = {'cpu.max': str((16 if index < 2 else 2)*100000) + ' 100000',
                  'memory.max': str((200 if index < 2 else 32)*2**30), 'memory.swap.max': '0'}
        handles.append(observe(path, limits))
        bind(pins, Path(path))
    dependency = observe(plan['numerical_closure_launch'])
    fit = json.loads(Path(plan['fit_plan']).read_text())
    completion = Path(fit['qualification_completion'])
    root = Path(plan['output'])
    native = [child for handle in handles for child in handle.get('child_resource_observations', [])
              if 'scripts/full_weighted_timing_parallel_v1.py' in child.get('cmdline', [])]
    if not completion.exists():
        assert not native and not root.exists(), 'Timing was admitted before full numerical closure'
    else:
        closed = json.loads(completion.read_text())
        assert sha(closed['full_hash_archive']) == closed['full_hash_archive_sha256']
        bind(pins, completion)
    assert fit['launch_state'] == 'not_launched_or_queued' and not Path(fit['output']).exists()
    for path in [Path(__file__), pp, ip, Path(plan['numerical_closure_launch'])]:
        bind(pins, path)
    verify(pins)
    result = dict(status='verified_original_gated_full_parallel_timing_runtime',
        checked_utc=datetime.now(timezone.utc).isoformat(), source_plan_sha256=sha(pp),
        original_handles=handles, exact_original_numerical_dependency=dependency,
        numerical_completion_present=completion.exists(), native_timing_process_count=len(native),
        output_root_present=root.exists(),
        producer_checkpoints=len(list((root/'checkpoints').glob('*.producer.json'))),
        reader_checkpoints=len(list((root/'checkpoints').glob('*.reader.json'))),
        failure_files=[str(p) for p in (root/'failures').rglob('*') if p.is_file()],
        source_hashes=pins, candidate_rows=20832000, maximum_eligible_groups=694400,
        working_model_fits_computed=0, production_finish_eta=None, scientific_eligibility=False,
        gpu=False, all_eight_aims_incomplete=True,
        scope='Exact original controller identities/cgroups or invocation-linked terminal evidence '
              'are observed, including the original numerical dependency. Absence of numerical '
              'closure requires no native timing child or output root. Artifact census and compact '
              'closure hashes are status evidence only; they do not independently qualify full '
              'arithmetic, timings, optimization, ancestral uncertainty or biological findings.')
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items()
                     if k not in ['source_hashes', 'original_handles', 'exact_original_numerical_dependency']}, indent=2))


if __name__ == '__main__':
    main()
