#!/usr/bin/env python3
"""Close every original state/length/category disposition and preserve failures."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import subprocess
import sys
from background_measurement_union_sources import closed_source
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True); a = p.parse_args()
    plan = json.loads(a.plan.read_text()); bindings = dict(plan['pins']); bind(bindings, a.plan)
    horizon = closed_source(plan['horizon_completion'], 'complete_verified_baliphy_initial_horizon_accounting',
        'complete_verified_baliphy_initial_horizon_accounting_archive', 2, bindings)
    pp = Path(plan['producer_plan']); source = json.loads(pp.read_text()); bind(bindings, pp)
    jobs = json.loads(Path(source['jobs']).read_text()); bind(bindings, source['jobs'])
    job_index = {j['chain']['chain_id']: j for j in jobs}; assert len(job_index) == len(jobs) == 1620
    native = json.loads((Path(source['output']) / 'receipt.json').read_text())
    native_status = {r['chain_id']: r for r in native['chains']}
    checked = {c for c, r in native_status.items() if r['status'] == 'all_saved_alignments_and_candidate_nodes_checked'}
    groups = defaultdict(set)
    for cid, job in job_index.items(): groups[job['config']['model_input_identity']].add(cid)
    complete_groups = {g for g, cids in groups.items() if cids <= checked}
    assert len(checked) == horizon['checked_chains'] == 1617 and len(complete_groups) == horizon['complete_quartets'] == 402
    def document(path, digest=None):
        path = Path(path); bind(bindings, path, digest); r = json.loads(path.read_text())
        for key in ['pins', 'evidence', 'source_hashes']:
            for q, d in r.get(key, {}).items(): bind(bindings, q, d)
        for n, d in r.get('artifacts', {}).items(): bind(bindings, path.parent / n, d)
        return r
    configs = {name: document(path) for name, path in plan['postprocessing_plans'].items()}
    reports = {name: document(Path(config['output']) / 'receipt.json') for name, config in configs.items()}
    for name, r in reports.items(): assert r['plan_sha256'] == sha(plan['postprocessing_plans'][name])
    traces = reports['states']; assert traces['status'] == 'terminal_full_state_trace_accounting_requires_review'
    assert set(traces['chains']) == set(job_index) and traces['processed_chains'] == 1617 and traces['unresolved_chains'] == 3
    state_observations = anchors = unanchored = 0
    for cid, x in traces['chains'].items():
        if cid not in checked:
            assert x['status'] == 'unresolved_no_verified_terminal_chain'; continue
        assert x['status'] == 'state_trace_complete_not_posterior_qualification'
        attempt = document(x['attempt_receipt'], x['attempt_receipt_sha256']); assert attempt['exit_code'] == 0
        r = document(x['receipt'], x['receipt_sha256'])
        assert r['status'] == 'complete_independently_checked_anchored_state_traces' and r['chains'] == 1
        assert len(r['summaries']) == 1; summary = r['summaries'][0]
        assert summary['chain'] == cid and summary['state_observations'] == x['state_observations'] == r['state_observations']
        assert summary['samples'] == 101 and summary['state_observations'] == summary['candidate_anchor_coordinates'] * 101
        assert summary['source_audit_sha256'] == native_status[cid]['sample_audit_sha256']
        bind(bindings, summary['source_audit'], summary['source_audit_sha256'])
        for q, d in summary['artifacts'].items(): bind(bindings, q, d)
        state_observations += summary['state_observations']; anchors += summary['candidate_anchor_coordinates']; unanchored += summary['unanchored_residue_observations']
    length = reports['length']; assert length['status'] == 'terminal_length_group_accounting_requires_review'
    assert set(length['groups']) == set(groups) and length['complete_quartets'] == 402 and length['unresolved_quartets'] == 3
    length_counts = Counter(); length_fully_passing = Counter()
    for group, x in length['groups'].items():
        if group not in complete_groups:
            assert x['status'] == 'unresolved_no_complete_scalar_quartet'; continue
        assert x['status'] == 'length_screens_complete_pending_review'
        r = document(x['receipt'], x['receipt_sha256']); document(x['scalar_receipt'], x['scalar_receipt_sha256'])
        assert r['status'] == 'candidate_length_screens_complete_not_posterior_qualification' and r['model_input_identity'] == group
        assert r['producer_plan_sha256'] == sha(pp) and set(r['outputs']) == {'250', '500'}
        for cutoff, info in r['outputs'].items():
            assert info['retained_samples_per_chain'] == (75 if cutoff == '250' else 50)
            d = document(info['path'], info['sha256']); assert d['status'] == 'scalar_diagnostics_complete_not_posterior_qualification'
            length_counts.update(cutoff + ':' + v['status'] for v in d['variables'].values())
            length_fully_passing[cutoff] += all(v['status'] == 'passes_scalar_screen_only' for v in d['variables'].values())
    categorical = reports['categorical']; assert categorical['status'] == 'terminal_full_categorical_report_accounting_requires_review'
    assert set(categorical['groups']) == set(groups) and categorical['completed_groups'] == 402 and categorical['unresolved_groups'] == 3
    coordinate_counts = Counter(); pattern_counts = Counter(); coordinates = Counter(); patterns = Counter()
    for group, x in categorical['groups'].items():
        if group not in complete_groups:
            assert x['status'] == 'unresolved_no_verified_extracted_quartet'; continue
        assert x['status'] == 'categorical_reports_complete_not_posterior_qualification'
        attempt = document(x['attempt_receipt'], x['attempt_receipt_sha256']); assert attempt['exit_code'] == 0
        r = document(x['receipt'], x['receipt_sha256']); assert r['status'] == 'verified_quartet_categorical_reports_complete_not_posterior_qualification' and r['group'] == group
        report = document(r['report'], r['report_sha256']); manifest = document(r['manifest'], r['manifest_sha256'])
        assert report['status'] == 'both_cutoff_categorical_reports_complete_not_posterior_qualification'
        assert manifest['status'] == 'provenance_checked_state_quartet'
        assert {c['chain_id'] for c in manifest['chains']} == groups[group] and len(manifest['chains']) == 4
        assert len({c['seed'] for c in manifest['chains']}) == 4 and all(c['model_input_identity'] == group for c in manifest['chains'])
        assert manifest['expected_iterations'] == list(range(0, 1001, 10))
        bind(bindings, manifest['arrays'], manifest['arrays_sha256'])
        for c in manifest['chains']: bind(bindings, c['log'], c['log_sha256'])
        assert set(report['outputs']) == {'250', '500'}
        for cutoff, info in report['outputs'].items():
            sp = Path(r['report']).parent / ('discard-' + cutoff) / 'summary.json'; summary = document(sp)
            assert summary['discard_through'] == int(cutoff)
            assert summary['retained_samples_per_chain'] == info['retained_samples_per_chain'] == (75 if cutoff == '250' else 50)
            assert summary['coordinates'] == info['coordinates'] == len(manifest['coordinates'])
            assert summary['patterns'] == info['patterns']
            assert sum(summary['coordinate_status_counts'].values()) == summary['coordinates'] and sum(summary['pattern_status_counts'].values()) == summary['patterns']
            coordinates[cutoff] += summary['coordinates']; patterns[cutoff] += summary['patterns']
            coordinate_counts.update({cutoff + ':' + k: v for k, v in summary['coordinate_status_counts'].items()})
            pattern_counts.update({cutoff + ':' + k: v for k, v in summary['pattern_status_counts'].items()})
        print('full_ancestral_postprocessing_group', group, flush=True)
    verify(bindings)
    summary = dict(full_native_chains=1620, full_quartets=405, state_chains=1617, unresolved_state_chains=3,
        state_observations=state_observations, per_chain_anchor_occurrences=anchors, unanchored_residue_observations=unanchored,
        length_quartets=402, categorical_quartets=402, unresolved_quartets=3, length_status_counts=dict(length_counts),
        quartets_passing_every_length_scalar=dict(length_fully_passing), categorical_coordinate_counts=dict(coordinates),
        categorical_pattern_counts=dict(patterns), coordinate_status_counts=dict(coordinate_counts), pattern_status_counts=dict(pattern_counts))
    output = Path(plan['output']); output.mkdir(exist_ok=False); snapshot = output / 'accounting.json'
    snapshot.write_text(json.dumps(dict(status='verified_full_ancestral_postprocessing_accounting_pending_journals',
        **summary, source_hashes=bindings, scientific_eligibility=False, scope=plan['scope']), indent=2) + '\n')
    inner = output / 'closure_plan.json'; archive = output / 'completion_archive.json'
    inner.write_text(json.dumps(dict(output=str(archive), completed_status='complete_verified_baliphy_postprocessing_accounting_archive',
        evidence=dict(accounting=dict(path=str(snapshot), expected=dict(status='verified_full_ancestral_postprocessing_accounting_pending_journals', **summary))),
        links=[], launches=plan['launches'], pins={**bindings, str(snapshot): sha(snapshot)}, summary=summary, scope=plan['scope']), indent=2) + '\n')
    subprocess.run([sys.executable, 'scripts/record_completed_process_handoffs_v2.py', '--plan', str(inner)], check=True)
    proof = json.loads(archive.read_text()); assert len(proof['services']) == 3
    result = dict(status='complete_verified_baliphy_postprocessing_accounting', **summary,
        full_hash_archive=str(archive), full_hash_archive_sha256=sha(archive), bound_source_hashes=len(proof['source_hashes']),
        exact_process_journals_checked=3, scientific_eligibility=False, scope=plan['scope'])
    with Path(plan['completion']).open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__': main()
