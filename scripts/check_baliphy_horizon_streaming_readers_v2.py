#!/usr/bin/env python3
"""Qualify dynamic readers on synthetic horizons and the complete current corpus."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import copy
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import resource
import shutil

from ancestral_chain_attempt import sha
from baliphy_horizon_streaming_readers_v2 import (horizon, scalar_trace, latent_trace, frame_context,
    joint_frames, fasta_trace, verify_existing_exports, verify_staged_exports, stage_exports)
from readback_baliphy_log_alpha_v10 import HEADER
from reference_measurement_union_sources import bind


SUCCESS = 'full_joint_short_sampler_output_integrity_checked_not_posterior'
REVIEW = 'explicit_nonfinite_or_literal_null_scalar_output_retained_for_review'


def reject(action):
    try:
        action()
    except (AssertionError, ValueError, KeyError, TypeError, ArithmeticError, RuntimeError, OSError):
        return
    raise AssertionError('Invalid horizon/trace/export was accepted')


def controls(root):
    root.mkdir()
    context = dict(observed={k: 'A' for k in 'abcde'}, runtime_labels=set('abcde') | {'n1', 'n2', 'n3', 'n4'},
        candidates={k: k for k in ['n1', 'n2', 'n3', 'n4']}, negative_branch_pairs_require_review=0)
    stats = dict(prior=-2, likelihood=-3, posterior=-5,
        __project_scalar_v6_quality__=dict(numericLeafCount=3, nonfinite=[], literalNullPaths=[]))
    scalar = {'iter': 0, 'statistics//': stats, 'parameters//': {'S1/': {'ASRV.Gamma:alpha': 1}},
              'numericParameterQuality//': dict(numericLeafCount=1, nonfinite=[], literalNullPaths=[])}
    latent = {'iter': 0, 'statistics//': stats, 'parameters//': {'S1/': dict(latentLogAlpha=0, derivedAlpha=1,
        categoryRates=[1, 1, 1, 1], laplaceLocation=0, laplaceScale=2, latentLogDensity=-math.log(4))},
        'numericParameterQuality//': dict(numericLeafCount=9, nonfinite=[], literalNullPaths=[])}
    frame = dict(iter=0, catStates={k: dict(states=[0], categories=[0]) for k in sorted(context['runtime_labels'])},
        alignmentLines=[line for k in sorted(context['runtime_labels']) for line in ['>' + k, 'A']],
        properties={'rate': [[1] * 20 for _ in range(4)]}, conditions={})
    mapping = {'iter': 'iter', 'prior': 'prior', 'likelihood': 'likelihood', 'posterior': 'posterior',
               'S1/ASRV.Gamma:alpha': 'alpha'}
    positives = []

    def fixture(name, iterations):
        directory = root / name
        directory.mkdir()
        (directory / 'runtime-tree.nwk').write_text('((((a:0.1,b:0.1)n1:0.1,c:0.1)n2:0.1,d:0.1)n3:0.1,e:0.1)n4;\n')
        (directory / 'C1.log.column-map.json').write_text(json.dumps(mapping) + '\n')
        with (directory / 'C1.log.json').open('x') as out, (directory / 'C1.P1.log-alpha-samples.jsonl').open('x') as diagnostic, (directory / 'C1.log').open('x') as tsv:
            out.write(json.dumps(HEADER) + '\n')
            diagnostic.write(json.dumps(HEADER) + '\n')
            tsv.write('iter\tprior\tlikelihood\tposterior\talpha\n')
            for iteration in range(iterations + 1):
                out.write(json.dumps(dict(scalar, iter=iteration)) + '\n')
                diagnostic.write(json.dumps(dict(latent, iter=iteration)) + '\n')
                tsv.write(str(iteration) + '\t-2\t-3\t-5\t1\n')
        with (directory / 'C1.P1.site-property-samples.jsonl').open('x') as out, (directory / 'C1.P1.fastas').open('x') as fasta:
            for iteration in range(0, iterations + 1, 10):
                out.write(json.dumps(dict(frame, iter=iteration)) + '\n')
                fasta.write('iterations = ' + str(iteration) + '\n' + '\n'.join(frame['alignmentLines']) + '\n')
        return directory

    short = fixture('synthetic-horizon17', 17)
    long = fixture('synthetic-horizon10000', 10000)
    for directory, iterations, states, frames in [(short, 17, 18, 2), (long, 10000, 10001, 1001)]:
        result = stage_exports(directory, context, iterations, root / ('sealed-' + str(iterations)), 'broad')
        assert result['status'] == 'complete_horizon_output_integrity_not_posterior'
        assert result['saved_joint_frames'] == frames and result['scalar']['rows'] == result['latent']['rows'] == states
        assert result['fasta']['saved_alignments'] == frames
        checked = verify_staged_exports(directory, context, iterations, root / ('sealed-' + str(iterations)), 'broad')
        assert checked['saved_joint_frames'] == frames and checked['reference_ledger_streamed']
        assert not result['scientific_eligibility'] and not result['posterior_qualified'] and not result['native_custody_qualified_here']
        positives.append(dict(kind='synthetic_parser_control_not_native_sampling', iterations=iterations, scalar_states=states,
            saved_joint_frames=frames, transactional_final_directory_created=True, scientific_eligibility=False))
    failures = []
    for name, action in [('boolean_horizon', lambda: horizon(True)), ('negative_horizon', lambda: horizon(-1)),
                         ('zero_interval', lambda: horizon(17, 0)), ('boolean_interval', lambda: horizon(17, True))]:
        reject(action)
        failures.append(name)

    def altered(name, filename, modify, action):
        directory = root / name
        shutil.copytree(short, directory)
        p = directory / filename
        p.write_text(modify(p.read_text()))
        reject(lambda: action(directory))
        failures.append(name)
        return directory

    def last_row(text, mutate):
        lines = text.splitlines()
        row = json.loads(lines[-1])
        mutate(row)
        lines[-1] = json.dumps(row)
        return '\n'.join(lines) + '\n'

    scalar_action = lambda d: scalar_trace(d, 17)
    altered('missing_scalar_last', 'C1.log.json', lambda t: '\n'.join(t.splitlines()[:-1]) + '\n', scalar_action)
    altered('extra_scalar_row', 'C1.log.json', lambda t: t + t.splitlines()[-1] + '\n', scalar_action)
    altered('unordered_scalar_last', 'C1.log.json', lambda t: last_row(t, lambda r: r.update(iter=16)), scalar_action)
    altered('boolean_scalar_last', 'C1.log.json', lambda t: last_row(t, lambda r: r.update(iter=True)), scalar_action)
    altered('duplicate_scalar_key', 'C1.log.json', lambda t: t.replace('"iter": 17', '"iter": 17, "iter": 17'), scalar_action)
    altered('untagged_nonfinite_scalar', 'C1.log.json', lambda t: t.replace('"ASRV.Gamma:alpha": 1', '"ASRV.Gamma:alpha": NaN'), scalar_action)
    altered('changed_finite_scalar', 'C1.log.json', lambda t: last_row(t, lambda r: r['parameters//']['S1/'].update({'ASRV.Gamma:alpha': 2})), scalar_action)
    altered('missing_tsv_row', 'C1.log', lambda t: '\n'.join(t.splitlines()[:-1]) + '\n', scalar_action)
    altered('duplicate_tsv_column', 'C1.log', lambda t: t.replace('posterior\talpha', 'posterior\tprior', 1), scalar_action)
    altered('invalid_scalar_header', 'C1.log.json', lambda t: t.replace('"version": "0.2"', '"version": "0.3"', 1), scalar_action)
    joint_action = lambda d: sum(1 for _ in joint_frames(d, context, 17))
    altered('missing_joint_last', 'C1.P1.site-property-samples.jsonl', lambda t: t.splitlines()[0] + '\n', joint_action)
    altered('extra_joint_frame', 'C1.P1.site-property-samples.jsonl', lambda t: t + t.splitlines()[-1] + '\n', joint_action)
    altered('unordered_joint_last', 'C1.P1.site-property-samples.jsonl', lambda t: last_row(t, lambda r: r.update(iter=0)), joint_action)
    altered('bad_late_category', 'C1.P1.site-property-samples.jsonl', lambda t: last_row(t, lambda r: r['catStates']['n4']['categories'].__setitem__(0, 4)), joint_action)
    altered('duplicate_joint_key', 'C1.P1.site-property-samples.jsonl', lambda t: t.replace('"iter": 10', '"iter": 10, "iter": 10'), joint_action)
    altered('missing_fasta_last', 'C1.P1.fastas', lambda t: t.split('iterations = 10')[0], lambda d: fasta_trace(d, context, 17))
    altered('changed_fasta_tip', 'C1.P1.fastas', lambda t: t.replace('>a\nA', '>a\nV'), lambda d: fasta_trace(d, context, 17))
    altered('missing_latent_last', 'C1.P1.log-alpha-samples.jsonl', lambda t: '\n'.join(t.splitlines()[:-1]) + '\n', lambda d: latent_trace(d, 17, 'broad'))
    altered('changed_latent_density', 'C1.P1.log-alpha-samples.jsonl', lambda t: last_row(t, lambda r: r['parameters//']['S1/'].update(latentLogDensity=0)), lambda d: latent_trace(d, 17, 'broad'))
    reject(lambda: latent_trace(short, 17, 'centered'))
    failures.append('wrong_latent_prior')
    bad = root / 'bad_late_category'
    final = root / 'late-failure-final'
    reject(lambda: stage_exports(bad, context, 17, final, 'broad'))
    pending = final.with_name(final.name + '.pending')
    assert pending.is_dir() and not final.exists() and not (pending / 'complete.json').exists()
    failure = json.loads((pending / 'failure.json').read_text())
    assert failure['complete_joint_frames_written'] == 1 and len(list(pending.glob('*.npz'))) == 1
    failures.append('late_failure_retains_pending_without_final_export')
    rejected_role = root / 'tagged_review'
    shutil.copytree(short, rejected_role)
    p = rejected_role / 'C1.log.json'
    lines = p.read_text().splitlines()
    row = json.loads(lines[2])
    row['parameters//']['S1/']['ASRV.Gamma:alpha'] = '__project_scalar_v6__:positive_infinity'
    row['numericParameterQuality//']['nonfinite'] = [dict(path=['S1/', 'ASRV.Gamma:alpha'], kind='positive_infinity')]
    lines[2] = json.dumps(row)
    p.write_text('\n'.join(lines) + '\n')
    result = stage_exports(rejected_role, context, 17, root / 'review-final', 'broad', require_latent=False)
    assert result['status'] == 'numeric_review_no_exports_created' and result['scalar']['nonfinite_reviews'] == 1
    assert not (root / 'review-final').exists() and not (root / 'review-final.pending').exists()
    positives.append(dict(kind='tagged_numeric_review_retained', nonfinite_rows=1, arrays_created=0))
    reject(lambda: stage_exports(short, context, 17, root / 'sealed-17', 'broad'))
    failures.append('existing_namespace_not_reused')
    for kind in ['extra_ledger_row', 'missing_array', 'foreign_file']:
        directory = root / ('altered-seal-' + kind)
        shutil.copytree(root / 'sealed-17', directory)
        if kind == 'extra_ledger_row':
            p = directory / 'frames.jsonl'
            p.write_text(p.read_text() + p.read_text().splitlines()[-1] + '\n')
            p = directory / 'complete.json'
            receipt = json.loads(p.read_text())
            receipt['frame_ledger_sha256'] = sha(directory / 'frames.jsonl')
            p.write_text(json.dumps(receipt) + '\n')
        elif kind == 'missing_array':
            (directory / 'frame-10.npz').unlink()
        else:
            (directory / 'invented.npz').write_bytes(b'foreign')
        reject(lambda: verify_staged_exports(short, context, 17, directory, 'broad'))
        failures.append('staged_readback_' + kind)
    return dict(synthetic_controls=positives, rejection_controls=failures,
        retained_pending_failure=failure, native_sampler_launched=False, scientific_eligibility=False)


def native_role(item):
    phase, job, row, mapping = item
    chain = job['chain']
    path = Path(row['native_receipt'])
    assert sha(path) == row['native_receipt_sha256']
    native = json.loads(path.read_text())
    encoded = json.dumps(job['config'], sort_keys=True, separators=(',', ':'), allow_nan=False).encode()
    assert native['configuration_sha256'] == hashlib.sha256(encoded).hexdigest()
    config = path.parent.parent / 'configuration.json'
    assert json.loads(config.read_text()) == job['config']
    assert json.loads((path.parent / 'command.json').read_text()) == job['config']['command']
    assert json.loads((path.parent / 'process.json').read_text())['command'] == job['config']['command']
    assert job['config']['command'][job['config']['command'].index('--iterations') + 1] == '20'
    pins = {str(path): row['native_receipt_sha256'], str(config): sha(config)}
    for name, h in native['artifacts'].items():
        p = path.parent / name
        assert sha(p) == h
        bind(pins, p, h)
    assert native['exit_code'] == row['exit_code']
    assert all(row[k] == chain[{'chain_role': 'chain'}.get(k, k)] for k in
               ['chain_id', 'effective_input_group', 'prior_label', 'chain_role', 'family', 'seed'])
    result = dict(phase=phase, chain_id=chain['chain_id'], family=chain['family'],
        effective_input_group=chain['effective_input_group'], prior_label=chain['prior_label'],
        chain_role=chain['chain'], inherited_integrity_disposition=row['status'], native_exit_code=native['exit_code'],
        scalar_rows=0, scalar_mapped_values=0, scalar_review_events=0, literal_null_rows=0,
        latent_rows=0, latent_overflow_rows=0, fasta_frames=0, joint_frames=0, joint_tip_pairs=0, joint_ancestral_pairs=0,
        original_failure_retained=native['exit_code'] != 0, original_numeric_review_retained=row['status'] == REVIEW,
        new_native_runs=0, scientific_eligibility=False, posterior_qualified=False)
    if native['exit_code'] != 0:
        assert not row['joint_frames']
        return result, pins
    directories = list(path.parent.glob('independent-chain-*'))
    assert len(directories) == 1
    directory = directories[0]
    scalar = scalar_trace(directory, 20)
    old = row['scalar_v6_audit']
    assert scalar['rows'] == old['rows'] == 21
    assert scalar['nonfinite_reviews'] == len(old['nonfinite_reviews'])
    assert scalar['literal_null_rows'] == len(old['literal_null_iterations'])
    assert scalar['mapped_values_compared'] == old['mapped_values_compared']
    assert scalar['fully_finite_mapped_integrity_checked'] == old['fully_finite_mapped_integrity_checked']
    result.update(scalar_rows=scalar['rows'], scalar_mapped_values=scalar['mapped_values_compared'],
        scalar_review_events=scalar['nonfinite_reviews'], literal_null_rows=scalar['literal_null_rows'])
    if phase == 'V10_comparison':
        latent = latent_trace(directory, 20, chain['prior_label'])
        result.update(latent_rows=latent['rows'], latent_overflow_rows=latent['overflow_rows'])
    context = frame_context(chain, directory, mapping)
    for name in ('alignment', 'tree'):
        bind(pins, chain[name], chain[name + '_sha256'])
    fasta = fasta_trace(directory, context, 20)
    assert fasta['saved_alignments'] == row['legacy_saved_alignments'] == 3
    assert fasta['candidate_frames'] == row['legacy_candidate_frames'] == 12
    result['fasta_frames'] = fasta['saved_alignments']
    if row['status'] == SUCCESS:
        assert scalar['fully_finite_mapped_integrity_checked']
        joint = verify_existing_exports(directory, context, 20, row['joint_frames'])
        result.update(joint_frames=joint['saved_joint_frames'], joint_tip_pairs=joint['joint_tip_pairs'],
                      joint_ancestral_pairs=joint['joint_ancestral_pairs'])
        for frame in row['joint_frames']:
            bind(pins, frame['projection_array'], frame['projection_array_sha256'])
    else:
        assert row['status'] == REVIEW and not scalar['fully_finite_mapped_integrity_checked'] and not row['joint_frames']
    return result, pins


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('output', 'receipt'):
        parser.add_argument('--' + name, required=True, type=Path)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    assert args.workers == 4 and not args.output.exists() and not args.receipt.exists()
    args.output.mkdir(parents=True)
    pins, expected = {}, {}
    for name in ['metadata/current_ancestral_horizon_cost_readback_transport_20261006_v1.json',
                 'metadata/baliphy_log_alpha_v10_full_grid_readback_transport_20261006_v1.json',
                 'metadata/ancestral_derived_horizon_cost_completed_20261006_v1.json']:
        p = Path(name)
        proof = json.loads(p.read_text())
        if 'original_tool_terminal_exit_code' in proof:
            assert proof['original_tool_terminal_exit_code'] == 0
        else:
            assert proof['status'] == 'complete_verified_full_ancestral_derived_storage_census'
        bind(pins, p)
        for name, h in proof['source_hashes'].items():
            key = str(Path(name).resolve())
            assert key not in expected or expected[key] == h
            expected[key] = h

    def closed(path):
        p = Path(path)
        h = expected[str(p.resolve())]
        assert sha(p) == h
        bind(pins, p, h)
        return json.loads(p.read_text())

    completion = closed('metadata/baliphy_scalar_v6_short_sampler_completed_20261004_v1.json')
    archive = closed(completion['full_hash_archive'])
    assert sha(completion['full_hash_archive']) == completion['full_hash_archive_sha256']
    assert len(archive['source_hashes']) == completion['bound_source_hashes']
    for name, h in archive['source_hashes'].items():
        key = str(Path(name).resolve())
        assert key not in expected or expected[key] == h
        expected[key] = h
    mapping = 'results/ancestral/case-local-trees-20260927-v1/ancestral_node_mapping.tsv'
    assert sha(mapping) == expected[str(Path(mapping).resolve())]
    bind(pins, mapping)
    work = []
    for phase, name in [('V6_original', 'metadata/baliphy_scalar_v6_sampler_execution_plan_20261004_v1.json'),
                        ('V10_comparison', 'metadata/baliphy_log_alpha_v10_full_grid_plan_20261005_v1.json')]:
        plan = closed(name)
        jobs = {j['chain']['chain_id']: j for j in closed(plan['jobs'])}
        rows = closed(Path(plan['output']) / 'dispositions.json')
        assert len(rows) == len(jobs) == (1620 if phase == 'V6_original' else 24)
        assert {p.name for p in (Path(plan['output']) / 'frames').iterdir()} == {r['chain_id'] for r in rows if r['status'] == SUCCESS}
        for row in rows:
            work.append((phase, jobs[row['chain_id']], row, mapping))
    assert len(work) == 1644
    # Each distinct native/input/source pin is checked once, rather than hashing
    # the same installed binary and libraries in every worker.
    for _, job, _, _ in work:
        for name, h in job['config']['pins'].items():
            bind(pins, name, h)
    for name, h in pins.items():
        assert sha(name) == h
    synthetic = controls(args.output / 'controls')
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for row, checked in pool.map(native_role, work, chunksize=1):
            for name, h in checked.items():
                assert expected[str(Path(name).resolve())] == h
                bind(pins, name, h)
            rows.append(row)
    assert len(rows) == 1644 and len({(r['phase'], r['chain_id']) for r in rows}) == 1644
    assert sum(r['original_failure_retained'] for r in rows) == 24
    assert sum(r['original_numeric_review_retained'] for r in rows) == 14
    assert sum(r['scalar_rows'] for r in rows) == 34020
    assert sum(r['latent_rows'] for r in rows) == 504 and sum(r['latent_overflow_rows'] for r in rows) == 8
    assert sum(r['fasta_frames'] for r in rows) == 4860 and sum(r['joint_frames'] for r in rows) == 4818
    path = args.output / 'role_readback.tsv'
    with path.open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    bind(pins, path)
    control_path = args.output / 'software_controls.json'
    with control_path.open('x') as handle:
        json.dump(synthetic, handle, indent=2)
        handle.write('\n')
    for p in (args.output / 'controls').rglob('*'):
        if p.is_file():
            bind(pins, p)
    for name in [Path(__file__), Path('scripts/baliphy_horizon_streaming_readers_v2.py'),
                 Path('scripts/baliphy_scalar_json_logger_v6c.py'), Path('scripts/read_baliphy_scalar_json_v6b.py'),
                 Path('scripts/readback_baliphy_log_alpha_v10.py'), Path('scripts/independent_joint_ancestral_frames.py'),
                 Path('scripts/independent_native_ancestral_alignment.py'), Path('scripts/independent_native_ancestral_topology.py'),
                 Path('scripts/independent_short_sampler_outputs_v2.py'), Path('scripts/reference_measurement_union_sources.py'),
                 Path('scripts/ancestral_chain_attempt.py'), control_path]:
        bind(pins, name)
    result = dict(status='passed_horizon_aware_streaming_readers_full_current_corpus_and_software_controls',
        checked_utc=datetime.now(timezone.utc).isoformat(), full_attempts_checked=1644, original_v6_roles=1620,
        separate_v10_comparisons=24, original_native_failures_retained=24, original_numeric_reviews_retained=14,
        finite_integrity_roles=1606, scalar_rows_checked=34020, mapped_scalar_values_checked=sum(r['scalar_mapped_values'] for r in rows),
        diagnostic_rows_checked=504, overflow_rows_retained=8, native_fasta_frames_checked=4860,
        existing_joint_exports_checked=4818, joint_tip_pairs_checked=sum(r['joint_tip_pairs'] for r in rows),
        joint_ancestral_pairs_checked=sum(r['joint_ancestral_pairs'] for r in rows),
        inherited_status_counts=dict(Counter(r['inherited_integrity_disposition'] for r in rows)),
        synthetic_positive_controls=synthetic['synthetic_controls'], rejection_controls=synthetic['rejection_controls'],
        synthetic_long_iterations=10000, synthetic_long_scalar_states=10001, synthetic_long_joint_frames=1001,
        maximum_retained_joint_frames_per_worker=1, reference_ledger_streaming_supported=True,
        transient_json_parser_allocation_not_a_long_chain_peak_bound=True, complete_trace_needed_before_final_export=True,
        original_native_arrays_unchanged=True, prior_short_reader_sources_unchanged=True,
        role_table=str(path), control_receipt=str(control_path), source_hashes=pins,
        producer_self_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        long_native_horizon_validated=False, native_custody_for_future_runs_qualified=False,
        long_memory_or_disk_bound_qualified=False, production_launch_allowed=False,
        scientific_eligibility=False, posterior_qualified=False, new_native_runs=0, new_predictions=0, gpu=False,
        scope='All current1644native attempts retain exact custody/outcomes. Dynamic scalar/FASTA/joint readers '
              'replay all native-zero scalar rows,24latent traces and every existing admitted export. '
              'Long horizons are software controls only, not MCMC or a pilot. Original qualified decoders '
              'and sources are shared dependencies; no independent biological truth or long-run '
              'peak/convergence/disk guarantee is established. Reviewed arrays remain excluded.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
