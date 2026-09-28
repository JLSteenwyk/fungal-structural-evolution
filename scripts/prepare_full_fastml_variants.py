#!/usr/bin/env python3
"""Bind all 156 inputs to both qualified software-regression variants."""
import argparse
import json
from pathlib import Path
import shutil
from ancestral_chain_attempt import sha, write_json


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    ap.add_argument('--inputs', type=Path, required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    assert not args.plan.exists() and not Path(args.output).exists()
    regression = Path('results/ancestral/fastml-cache-regression-20260928-v1/receipt.json')
    assert json.loads(regression.read_text())['status'] == 'passed_paired_cache_regression'
    bp = Path('results/software/fastml-precision-cache-build-20260928-v1/receipt.json')
    builds = json.loads(bp.read_text())['builds']
    sp = Path('metadata/ancestral_fastml_indel_plan_20260927.json')
    source = json.loads(sp.read_text())
    terminal = Path('metadata/fastml_indel_terminal_readback_20260928.json')
    proof = {r['job_id']: r for r in json.loads(terminal.read_text())['jobs']}
    pins = dict(source['pins'])
    for path in [str(sp), str(bp), str(regression), str(terminal), __file__,
                 'scripts/run_full_fastml_variants.py', 'scripts/readback_fastml_variant.py',
                 'scripts/ancestral_chain_attempt.py', 'scripts/replay_fastml_indel_probabilities.py']:
        pins[path] = sha(path)
    args.inputs.mkdir(parents=True, exist_ok=False)
    limiter = str(Path(shutil.which('prlimit')).resolve())
    jobs = []
    for job in source['jobs']:
        original = Path(source['output']) / job['job_id']
        assert sha(original / 'receipt.json') == proof[job['job_id']]['source_receipt_sha256']
        tree = (original / 'input_tree.nwk').resolve()
        seq = Path(job['characters']).resolve()
        for variant in builds:
            binary = str(Path(variant['binary']).resolve())
            assert sha(binary) == variant['binary_sha256']
            identity = job['job_id'] + '__' + variant['label']
            params = (args.inputs / (identity + '.txt')).resolve()
            settings = dict(source['settings'], _seqFile=str(seq), _treeFile=str(tree), _outDir='RESULTS')
            params.write_text(''.join(f'{k} {v}\n' for k, v in settings.items()))
            config = dict(command=[limiter, '--as=4294967296', binary, str(params)],
                timeout_seconds=43200, pins={p: sha(p) for p in [limiter, binary, str(params), str(tree), str(seq)]},
                model_input_identity=identity)
            pins.update(config['pins'])
            jobs.append(dict(id=identity, variant=variant['label'], job=job, config=config))
    assert len(jobs) == 312
    write_json(args.plan, dict(output=args.output, workers=4, jobs=jobs, pins=pins,
        minimum_free_disk_bytes=64 * 2**30,
        resources=dict(cpus=4, memory_gib=24, swap_gib=0, per_fit_address_space_gib=4,
            output_allowance_gib=32, planning_hours=[2, 24],
            basis='Original full batch used 5.2 CPU hours; two variants imply 10.4 CPU hours '
                  'before replay or changed optimizer behavior. Four workers; broad planning '
                  'allowance rather than an ETA. 12-hour per-fit timeout retains failures.',
            gpu=False, paid_resources=False),
        scope='Full paired inference, including empty inputs, plus independent replay. '
              'No globally optimized parameter or qualified ancestor claim from execution alone.'))


if __name__ == '__main__':
    main()
