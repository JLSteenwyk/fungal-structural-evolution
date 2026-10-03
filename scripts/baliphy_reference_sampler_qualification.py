"""Build and audit the complete short sampler/resource qualification grid."""
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path

from ancestral_chain_attempt import sha
from baliphy_horizon_resource_inventory import trace_summary
from baliphy_reference_initialization import transform
from readback_independent_baliphy_chain import readback


ITERATIONS = 20
SUMMARY_FIELDS = ['full_chains', 'full_quartets', 'effective_inputs', 'original_configuration_aliases',
                  'checked_sampler_attempts', 'unsuccessful_sampler_attempts', 'complete_quartets',
                  'unresolved_quartets', 'saved_alignments', 'candidate_frames', 'status_counts',
                  'allocation_warning_attempts', 'bad_alloc_attempts', 'gamma_equal_rate_limit_events',
                  'other_nonfinite_events', 'nonpositive_parameter_review_events', 'posterior_qualified']


def build_jobs(startups, census_chains, census_attempts, binary, prlimit, api_paths):
    assert len(startups) == len({x['chain']['chain_id'] for x in startups}) == 1620
    census = {x['chain_id']: x for x in census_chains}; assert len(census) == 1620
    risk_families = {census[x['chain_id']]['family'] for x in census_attempts
                     if x['allocation_warning_lines'] or x['bad_alloc']}
    assert risk_families
    jobs = []
    for startup in startups:
        source = startup['chain']; cid = source['chain_id']; old = census[cid]
        for key in ['seed', 'effective_input_group', 'prior_label', 'family', 'alignment_sha256',
                    'tree_sha256', 'original_configuration_ids', 'proteins']:
            assert source[key] == old[key], (cid, key)
        program = Path(startup['program']); assert program.read_text() == transform(Path(source['program']).read_text())
        memory = (48 if source['family'] in risk_families else 12) * 2**30
        baseline = old['selected_attempt']['elapsed_worker_seconds'] if old['selected_status'] == 'all_saved_alignments_and_candidate_nodes_checked' else None
        timeout = 7200 if baseline is None else max(900, math.ceil(baseline * ITERATIONS / 1000 * 4 + 60))
        assert timeout <= 7200
        chain = dict(source, seed=startup['fresh_seed'], program=str(program), program_sha256=sha(program))
        paths = [binary, prlimit, program, Path(source['tree']).resolve(), Path(source['alignment']).resolve(), *api_paths]
        config = dict(command=[str(prlimit), '--as=' + str(memory), '--cpu=' + str(timeout),
            '--fsize=' + str(2 * 2**30), '--', str(binary), '--seed', str(chain['seed']),
            'run', str(program), '--iterations', str(ITERATIONS), '--log-format', 'json,tsv',
            '--name', 'independent-chain'], timeout_seconds=timeout,
            pins={str(p.resolve()): sha(p) for p in paths})
        jobs.append(dict(chain=chain, source_seed=source['seed'], memory_reservation_bytes=memory,
            resource_class='48GiB_family_with_historical_allocation_warning' if memory > 12 * 2**30 else '12GiB_other_family',
            previous_successful_elapsed_worker_seconds=baseline,
            linear_twenty_iteration_worker_seconds=None if baseline is None else baseline * ITERATIONS / 1000,
            config=config))
    assert len({j['chain']['seed'] for j in jobs}) == 1620
    assert {j['chain']['seed'] for j in jobs}.isdisjoint({j['source_seed'] for j in jobs})
    assert Counter(j['chain']['prior_label'] for j in jobs) == {'broad': 540, 'centered': 540, 'package': 540}
    groups = defaultdict(list)
    for job in jobs:
        c = job['chain']; groups[c['effective_input_group'] + '-' + c['prior_label']].append(c)
        assert isinstance(c['seed'], int) and 1 <= c['seed'] <= 2**31 - 1
    assert len(groups) == 405 and len({j['chain']['effective_input_group'] for j in jobs}) == 135
    assert all(len(v) == 4 and {c['chain'] for c in v} == {1, 2, 3, 4} for v in groups.values())
    assert len({a for j in jobs for a in j['chain']['original_configuration_ids']}) == 324
    return jobs, sorted(risk_families)


def inspect(job, receipt_path, plan_digest, mapping):
    receipt_path = Path(receipt_path); native = json.loads(receipt_path.read_text()); chain = job['chain']
    encoded = json.dumps(job['config'], sort_keys=True, separators=(',', ':'), allow_nan=False)
    assert native['configuration_sha256'] == hashlib.sha256(encoded.encode()).hexdigest()
    assert json.loads((receipt_path.parent.parent / 'configuration.json').read_text()) == job['config']
    assert json.loads((receipt_path.parent / 'command.json').read_text()) == job['config']['command']
    assert json.loads((receipt_path.parent / 'process.json').read_text())['command'] == job['config']['command']
    for name, h in native['artifacts'].items(): assert sha(receipt_path.parent / name) == h
    row = dict(chain_id=chain['chain_id'], effective_input_group=chain['effective_input_group'],
        model_input_identity=chain['effective_input_group'] + '-' + chain['prior_label'],
        prior_label=chain['prior_label'], chain_role=chain['chain'], family=chain['family'],
        original_configuration_ids=chain['original_configuration_ids'], seed=chain['seed'],
        source_seed=job['source_seed'], native_receipt=str(receipt_path), native_receipt_sha256=sha(receipt_path),
        native_status=native['status'], exit_code=native['exit_code'], elapsed_worker_seconds=native['elapsed_seconds'],
        memory_reservation_bytes=job['memory_reservation_bytes'], plan_sha256=plan_digest,
        scientific_eligibility=False, posterior_qualified=False)
    stderr = (receipt_path.parent / 'stderr.log').read_text()
    row.update(allocation_warning_lines=stderr.count('Not enough memory to allocate a DP matrix'),
               bad_alloc='bad_alloc' in stderr)
    directories = list(receipt_path.parent.glob('independent-chain-*'))
    row['partial_scalar_trace'] = None
    row['partial_scalar_trace_issue'] = None
    if len(directories) == 1 and (directories[0] / 'C1.log').is_file():
        try:
            row['partial_scalar_trace'] = trace_summary(directories[0] / 'C1.log', False, ITERATIONS)
        except (AssertionError, KeyError, ValueError, ArithmeticError) as error:
            # Header-only or truncated failure logs must not erase the native
            # failure disposition or stop accounting for the remaining roles.
            row['partial_scalar_trace_issue'] = dict(error_type=type(error).__name__, error=str(error))
    if native['exit_code'] != 0:
        return dict(**row, status='unsuccessful_sampler_qualification_attempt_retained', saved_alignments=0, candidate_frames=0)
    try:
        assert len(directories) == 1
        complete = trace_summary(directories[0] / 'C1.log', True, ITERATIONS)
        row['partial_scalar_trace'] = complete
        audit = readback(chain, receipt_path, ITERATIONS, mapping)
        assert len(audit['candidate_samples']) == 12
        return dict(**row, status='full_short_sampler_output_integrity_checked_not_posterior',
                    saved_alignments=3, candidate_frames=12, sample_audit=audit)
    except (AssertionError, KeyError, ValueError, ArithmeticError, OSError) as error:
        return dict(**row, status='invalid_sampler_qualification_output_retained', saved_alignments=0,
                    candidate_frames=0, error_type=type(error).__name__, error=str(error))


def summarize(rows):
    from run_baliphy_reference_preflight import summarize as grid_summary
    grid = grid_summary([dict(r, status='reference_startup_homology_density_and_representation_checked'
        if r['status'] == 'full_short_sampler_output_integrity_checked_not_posterior' else 'unsuccessful_native_startup_retained') for r in rows])
    limit = other = support = 0
    for row in rows:
        trace = row['partial_scalar_trace']
        if trace is None: continue
        for name, param in trace['parameters'].items():
            for event in param['nonfinite']:
                if name == 'ASRV.Gamma:alpha' and event['token'] in ['Infinity', '+Infinity']:
                    limit += 1
                else: other += 1
            support += len(param['nonpositive_iterations'])
    return dict(full_chains=1620, full_quartets=405, effective_inputs=grid['effective_inputs'],
        original_configuration_aliases=grid['original_configuration_aliases'],
        checked_sampler_attempts=grid['validated_startups'], unsuccessful_sampler_attempts=grid['unsuccessful_startups'],
        complete_quartets=grid['complete_startup_quartets'], unresolved_quartets=grid['unresolved_startup_quartets'],
        saved_alignments=sum(r['saved_alignments'] for r in rows), candidate_frames=sum(r['candidate_frames'] for r in rows),
        status_counts=dict(Counter(r['status'] for r in rows)),
        allocation_warning_attempts=sum(bool(r['allocation_warning_lines']) for r in rows),
        bad_alloc_attempts=sum(r['bad_alloc'] for r in rows), gamma_equal_rate_limit_events=limit,
        other_nonfinite_events=other, nonpositive_parameter_review_events=support, posterior_qualified=False)
