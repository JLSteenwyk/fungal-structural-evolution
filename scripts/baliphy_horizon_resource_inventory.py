"""Full initial/selected-attempt resource and numerical survey; no new sampling."""
from collections import Counter, defaultdict
import csv
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import statistics

from run_ortholog_pair_guide_comparison import sha

PARAMETERS = ['RS07:rate', 'RS07:meanLength', 'ASRV.Gamma:alpha', '|A|', '|indels|']


def trace_summary(path, complete, horizon):
    """Stream every row, retaining invalid numeric tokens rather than imputing."""
    ranges = {k:dict(finite_min=None, finite_max=None, nonfinite=[], nonpositive_iterations=[]) for k in PARAMETERS}
    first = last = None; count = 0
    with Path(path).open() as f:
        reader = csv.DictReader(f, delimiter='\t')
        assert reader.fieldnames and len(reader.fieldnames) == len(set(reader.fieldnames))
        assert {'iter', *PARAMETERS} <= set(reader.fieldnames)
        for row in reader:
            assert None not in row and all(v is not None for v in row.values())
            iteration = int(row['iter']); assert iteration == count
            selected = dict(iteration=iteration, values={k:row[k] for k in PARAMETERS})
            if first is None: first = selected
            last = selected; count += 1
            for k in PARAMETERS:
                value = Decimal(row[k]); entry = ranges[k]
                if not value.is_finite() or not math.isfinite(float(value)):
                    entry['nonfinite'].append(dict(iteration=iteration, token=row[k])); continue
                number = float(value)
                entry['finite_min'] = number if entry['finite_min'] is None else min(entry['finite_min'], number)
                entry['finite_max'] = number if entry['finite_max'] is None else max(entry['finite_max'], number)
                if k != '|indels|' and value <= 0: entry['nonpositive_iterations'].append(iteration)
                if k == '|indels|' and value < 0: entry['nonpositive_iterations'].append(iteration)
    assert count > 0 and count <= horizon + 1
    if complete: assert count == horizon + 1
    alignment_ratio = None
    if first is not None:
        initial = Decimal(first['values']['|A|']); maximum = ranges['|A|']['finite_max']
        if initial.is_finite() and initial > 0 and maximum is not None:
            alignment_ratio = float(Decimal(str(maximum)) / initial)
    return dict(rows=count, first=first, last=last, parameters=ranges,
        maximum_alignment_to_initial_ratio=alignment_ratio,
        numerical_parameter_review=any(v['nonfinite'] or v['nonpositive_iterations'] for v in ranges.values()),
        complete_horizon_logged=count == horizon + 1)


def fresh_seeds(chains, horizon, namespace):
    used = {r['seed'] for r in chains}; result = {}
    assert len(used) == len(chains)
    for chain in sorted(chains, key=lambda r:r['chain_id']):
        salt = 0
        while True:
            value = f"{namespace}:{horizon}:{chain['chain_id']}:{salt}"
            seed = int.from_bytes(hashlib.sha256(value.encode()).digest()[:8], 'big') % (2**31 - 1) + 1
            if seed not in used: break
            salt += 1
        used.add(seed); result[chain['chain_id']] = seed
    return result


def aggregate(chains, attempts, horizon, namespace):
    assert len(chains) == len({r['chain_id'] for r in chains}) == 1620
    groups = defaultdict(list)
    for row in chains: groups[row['model_input_identity']].append(row)
    assert len(groups) == 405 and len({r['effective_input_group'] for r in chains}) == 135
    assert Counter(r['prior_label'] for r in chains) == {'broad':540, 'centered':540, 'package':540}
    seeds = fresh_seeds(chains, horizon, namespace)
    grid = []; proposed = []
    for group, rows in sorted(groups.items()):
        rows = sorted(rows, key=lambda r:r['chain_id'])
        assert len(rows) == len({r['seed'] for r in rows}) == 4
        assert {r['chain'] for r in rows} == {1,2,3,4}
        success = [r for r in rows if r['selected_status'] == 'all_saved_alignments_and_candidate_nodes_checked']
        review = len(success) != 4 or any(r['selected_attempt']['trace']['numerical_parameter_review'] or
            r['selected_attempt']['allocation_warning_lines'] for r in rows)
        measured = [r['selected_attempt']['elapsed_worker_seconds'] for r in success]
        grid.append(dict(group=group, chain_ids=[r['chain_id'] for r in rows], checked_chains=len(success),
            failed_chains=4-len(success), measured_successful_worker_seconds=sum(measured),
            median_successful_worker_seconds=statistics.median(measured) if measured else None,
            sampler_review_required=review, scientific_eligibility=False))
        for row in rows:
            baseline = row['selected_attempt']['elapsed_worker_seconds'] if row in success else None
            proposed.append(dict(source_chain_id=row['chain_id'], source_seed=row['seed'], fresh_seed=seeds[row['chain_id']],
                model_input_identity=group, original_configuration_ids=row['original_configuration_ids'],
                iterations=horizon, saved_alignment_interval=10, expected_saved_alignments=horizon//10+1,
                burn_in_cutoffs=[horizon//4,horizon//2], retained_alignment_samples=[horizon*3//40,horizon//20],
                linear_successful_worker_seconds=None if baseline is None else baseline*horizon/1000,
                failed_chain_memory_requirement_unresolved=row not in success,
                allocation_or_parameter_review_required=bool(row['selected_attempt']['allocation_warning_lines'] or
                    row['selected_attempt']['trace']['numerical_parameter_review']),
                production_launch_allowed=False, scientific_eligibility=False))
    assert horizon > 1000 and horizon % 40 == 0
    successful = [r for r in chains if r['selected_status'] == 'all_saved_alignments_and_candidate_nodes_checked']
    assert len(successful) == 1618 and len(attempts) == 1623
    assert sum(r['checked_chains']==4 for r in grid)==403
    assert sum(r['checked_chains']!=4 for r in grid)==2
    observed_seconds = sum(r['selected_attempt']['elapsed_worker_seconds'] for r in successful)
    observed_bytes = sum(r['selected_attempt']['native_alignment_bytes_size_only'] for r in successful)
    diagnostics = Counter()
    for row in attempts:
        diagnostics['allocation_warning_attempts'] += bool(row['allocation_warning_lines'])
        diagnostics['bad_alloc_attempts'] += row['bad_alloc']
        diagnostics['numerical_parameter_review_attempts'] += row['trace']['numerical_parameter_review']
        diagnostics['nonfinite_parameter_observations'] += sum(len(v['nonfinite']) for v in row['trace']['parameters'].values())
    summary = dict(full_chains=1620, effective_input_groups=135, full_quartets=405,
        selected_checked_chains=1618, selected_failed_chains=2, original_and_selected_recovery_attempts=1623,
        complete_quartets=sum(r['checked_chains']==4 for r in grid), unresolved_quartets=sum(r['checked_chains']!=4 for r in grid),
        proposed_fresh_chain_seeds=len(seeds), proposed_iterations=horizon,
        observed_successful_elapsed_worker_seconds=observed_seconds,
        observed_successful_native_alignment_bytes_size_only=observed_bytes,
        linear_successful_worker_seconds=observed_seconds*horizon/1000,
        linear_successful_native_alignment_bytes_size_only=math.ceil(observed_bytes*(horizon//10+1)/101),
        memory_requirements_for_two_failed_chains_unresolved=True,
        runtime_projection_uncalibrated=True, finish_eta=None, convergence_eta=None,
        sampler_review_quartets=sum(r['sampler_review_required'] for r in grid),
        selected_successful_allocation_warning_chains=sum(bool(r['selected_attempt']['allocation_warning_lines']) for r in successful),
        selected_successful_parameter_review_chains=sum(r['selected_attempt']['trace']['numerical_parameter_review'] for r in successful),
        attempt_diagnostics=dict(diagnostics), production_launch_allowed=False, scientific_eligibility=False)
    return grid, proposed, summary


def collect(plan):
    """Rehash used small inputs, logs and receipts; large sample files are sized only."""
    bindings = dict(plan['pins'])
    def doc(path, digest=None):
        path = str(path); actual = sha(path)
        if digest is not None: assert actual == digest, path
        if path in bindings: assert bindings[path] == actual, path
        bindings[path] = actual
        return json.loads(Path(path).read_text())
    for path, digest in bindings.items(): assert sha(path) == digest
    complete = doc(plan['diagnostic_completion'])
    assert complete['status'] == 'complete_verified_full_baliphy_recovery_diagnostics'
    proof = doc(complete['full_hash_archive'], complete['full_hash_archive_sha256'])
    assert len(proof['source_hashes']) == complete['bound_source_hashes'] and len(proof['services']) == 2
    expected = {str(Path(p).resolve()):d for p,d in proof['source_hashes'].items()}
    def bound_doc(path, digest=None):
        assert str(Path(path).resolve()) in expected, path
        if digest is not None: assert expected[str(Path(path).resolve())] == digest
        return doc(path, expected[str(Path(path).resolve())])
    producer = bound_doc(complete['producer_receipt'], complete['producer_receipt_sha256'])
    recovery = bound_doc(plan['recovery_completion'])
    overlay = bound_doc(recovery['overlay'], recovery['overlay_sha256'])
    assert overlay['original_samples_concatenated'] is False
    originals = doc(plan['chain_inputs']); assert len(originals) == 1620
    original_lookup = {r['chain_id']:r for r in originals}; assert len(original_lookup) == 1620
    inventory = doc(plan['native_inventory'])
    details = doc(inventory['full_chain_details'], inventory['full_chain_details_sha256'])
    scalar = doc(plan['scalar_completion'])
    assert scalar['status'] == 'complete_verified_full_independent_baliphy_scalar_comparison'
    scalar_proof = doc(scalar['full_hash_archive'], scalar['full_hash_archive_sha256'])
    assert len(scalar_proof['services']) == 2 and len(scalar_proof['source_hashes']) == scalar['bound_source_hashes']
    assert scalar['quartets_passing_every_scalar'] == scalar['quartets_passing_every_length_scalar'] == {'250':0,'500':0}
    chains = []; attempts = {}; seen = set()
    def measured(disposition, cid):
        p = Path(disposition['receipt']); receipt = bound_doc(p, disposition['receipt_sha256'])
        if str(p) in attempts: return attempts[str(p)]
        assert receipt['elapsed_seconds'] >= 0 and math.isfinite(receipt['elapsed_seconds'])
        configuration = bound_doc(p.parent.parent/'configuration.json')
        encoded = json.dumps(configuration,sort_keys=True,separators=(',',':'),allow_nan=False)
        assert hashlib.sha256(encoded.encode()).hexdigest() == receipt['configuration_sha256']
        chain = original_lookup[cid]
        assert configuration['seed'] == chain['seed'] and configuration['model_input_identity'] == chain['effective_input_group']+'-'+chain['prior_label']
        def artifact(name, optional=False):
            keys = [k for k in receipt['artifacts'] if Path(k).name == name]
            if optional and not keys: return None
            assert len(keys) == 1, (cid,name,keys)
            q = p.parent/keys[0]; digest = receipt['artifacts'][keys[0]]
            assert expected[str(q.resolve())] == digest and sha(q) == digest
            bindings[str(q)] = digest; return q
        command = json.loads(artifact('command.json').read_text())
        assert command == configuration['command'] and command[command.index('--seed')+1] == str(chain['seed'])
        assert command[command.index('--iterations')+1] == '1000'
        assert command[command.index('run')+1] == chain['program']
        stderr = artifact('stderr.log').read_text(); log = artifact('C1.log')
        success = disposition['status'] == 'all_saved_alignments_and_candidate_nodes_checked'
        assert success == (receipt['exit_code']==0 and receipt['status']=='exited_zero_pending_scientific_validation')
        native_keys = [k for k in receipt['artifacts'] if Path(k).name == 'C1.P1.fastas']; assert len(native_keys)==1
        native_path = p.parent/native_keys[0]
        row = dict(chain_id=cid, receipt=str(p),receipt_sha256=sha(p), disposition=disposition['status'],
            exit_code=receipt['exit_code'],elapsed_worker_seconds=receipt['elapsed_seconds'],
            address_space_limit_bytes=int(next(s.split('=',1)[1] for s in command if s.startswith('--as='))),
            allocation_warning_lines=sum('Allocation failed in sample_tri_multi!' in s for s in stderr.splitlines()),
            bad_alloc='std::bad_alloc' in stderr,
            scalar_log=str(log),scalar_log_sha256=sha(log),trace=trace_summary(log,success,1000),
            native_alignment_bytes_size_only=native_path.stat().st_size,
            native_alignment_content_hash_rechecked=False, memory_peak_unavailable=True)
        if success:
            entry=details[cid]; assert entry['selected_attempt_receipt_sha256']==sha(p)
            assert entry['native_files']['C1.P1.fastas']['bytes_at_observation']==native_path.stat().st_size
        attempts[str(p)]=row; return row
    for native in sorted(overlay['rows'],key=lambda r:r['chain']['chain_id']):
        chain=native['chain'];cid=chain['chain_id']; assert cid not in seen;seen.add(cid)
        assert chain==original_lookup[cid] and native['scientific_eligibility'] is False
        for field in ['alignment','tree','program']:
            path=chain[field];digest=chain[field+'_sha256'];assert sha(path)==digest;bindings[path]=digest
        original=measured(native['original_disposition'],cid)
        selected=measured(native['selected_disposition'],cid)
        chains.append(dict(**chain, model_input_identity=native['model_input_identity'],
            selected_status=native['selected_disposition']['status'],
            original_attempt=original,selected_attempt=selected,scientific_eligibility=False))
    assert seen==set(original_lookup)==set(details)
    assert {cid for info in producer['groups'].values() for cid in info['chain_ids']}==seen
    assert len(producer['groups'])==405
    for group,info in producer['groups'].items():
        assert len(info['chain_ids'])==len(set(info['chain_ids']))==4
        assert all(original_lookup[cid]['effective_input_group']+'-'+original_lookup[cid]['prior_label']==group for cid in info['chain_ids'])
    grid, proposed, summary=aggregate(chains,list(attempts.values()),plan['proposed_iterations'],plan['seed_namespace'])
    for path,digest in bindings.items(): assert sha(path)==digest,path
    return dict(chains=chains,attempts=[attempts[p] for p in sorted(attempts)],quartets=grid,proposed=proposed,summary=summary),bindings
