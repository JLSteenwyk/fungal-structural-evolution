#!/usr/bin/env python3
"""Close full initial-horizon accounting without qualifying posterior samples."""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True); a = p.parse_args()
    plan = json.loads(a.plan.read_text()); bindings = dict(plan['pins']); bind(bindings, a.plan)
    pp, dp = Path(plan['producer_plan']), Path(plan['diagnostic_plan'])
    producer_plan, diagnostic_plan = [json.loads(p.read_text()) for p in [pp, dp]]
    assert diagnostic_plan['producer_plan'] == str(pp)
    for config in [producer_plan, diagnostic_plan]:
        for path, digest in config['pins'].items(): bind(bindings, path, digest)
    root = Path(producer_plan['output']); diagnostics = Path(diagnostic_plan['output'])
    rp, drp = root / 'receipt.json', diagnostics / 'receipt.json'
    r, d = [json.loads(p.read_text()) for p in [rp, drp]]; bind(bindings, rp); bind(bindings, drp)
    assert r['status'] == 'all_initial_horizon_dispositions_recorded' and r['plan_sha256'] == sha(pp)
    assert d['status'] == 'terminal_group_accounting_requires_scientific_review' and d['plan_sha256'] == sha(dp)
    assert (root / 'run_plan.json').read_bytes() == pp.read_bytes(); bind(bindings, root / 'run_plan.json')
    jobs = json.loads(Path(producer_plan['jobs']).read_text()); bind(bindings, producer_plan['jobs'])
    job_index = {j['chain']['chain_id']: j for j in jobs}; assert len(jobs) == len(job_index) == 1620
    observed = {x['chain_id']: x for x in r['chains']}; assert set(observed) == set(job_index) and len(observed) == len(r['chains'])
    groups = defaultdict(set); counts = Counter(); failed = []; allocation_notices = []
    for cid, job in sorted(job_index.items()):
        chain, config = job['chain'], job['config']; groups[config['model_input_identity']].add(cid)
        for path, digest in config['pins'].items(): bind(bindings, path, digest)
        x = observed[cid]; attempt_path = Path(x['receipt']); bind(bindings, attempt_path, x['receipt_sha256'])
        attempt = json.loads(attempt_path.read_text())
        assert attempt_path.parent.parent == (root / cid).resolve()
        digest = hashlib.sha256(json.dumps(config, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
        assert attempt['configuration_sha256'] == digest
        for path, expected in [(root / cid / 'configuration.json', config), (attempt_path.parent / 'command.json', config['command'])]:
            assert json.loads(path.read_text()) == expected; bind(bindings, path)
        for name, digest in attempt['artifacts'].items(): bind(bindings, attempt_path.parent / name, digest)
        stderr = attempt_path.parent / 'stderr.log'; text = stderr.read_text(errors='replace')
        notices = text.count('Allocation failed in sample_tri_multi')
        if notices: allocation_notices.append(dict(chain_id=cid, status=x['status'], allocation_notices=notices, stderr=str(stderr), stderr_sha256=sha(stderr)))
        if x['status'] == 'all_saved_alignments_and_candidate_nodes_checked':
            assert attempt['exit_code'] == 0 and attempt['status'] == 'exited_zero_pending_scientific_validation'
            ap = Path(x['sample_audit']); bind(bindings, ap, x['sample_audit_sha256']); audit = json.loads(ap.read_text())
            assert audit['status'] == x['status'] and audit['attempt_receipt_sha256'] == sha(attempt_path)
            assert audit['iterations'] == producer_plan['iterations'] and audit['mapping_sha256'] == sha(producer_plan['mapping'])
            bind(bindings, audit['scalar_log'], audit['scalar_log_sha256'])
            assert len(audit['candidate_samples']) == 4 * (producer_plan['iterations'] // 10 + 1)
        else:
            assert x['status'] == attempt['status'] and attempt['exit_code'] != 0
            failed.append(dict(**x, exit_code=attempt['exit_code'], native_bad_alloc='std::bad_alloc' in text,
                allocation_notices=notices, address_space_limit=next(s for s in config['command'] if s.startswith('--as=')),
                stderr=str(stderr), stderr_sha256=sha(stderr)))
        counts[x['status']] += 1
    assert dict(counts) == r['counts'] and len(groups) == 405 and all(len(g) == 4 for g in groups.values())
    assert set(d['groups']) == set(groups)
    diagnostic_counts = Counter(); complete = unresolved = 0; fully_passing = Counter()
    for group, x in sorted(d['groups'].items()):
        ready = all(observed[c]['status'] == 'all_saved_alignments_and_candidate_nodes_checked' for c in groups[group])
        if not ready:
            assert x['status'] == 'unresolved_missing_checked_quartet' and set(x['chain_ids']) == groups[group]; unresolved += 1; continue
        assert x['status'] == 'both_burnin_scalar_screens_complete'; complete += 1
        gp = Path(x['receipt']); bind(bindings, gp, x['receipt_sha256']); gr = json.loads(gp.read_text())
        for path, digest in gr['evidence'].items(): bind(bindings, path, digest)
        assert set(gr['outputs']) == {str(int(producer_plan['iterations'] * b)) for b in producer_plan['diagnostics']['burn_in_fractions']}
        for burnin, info in gr['outputs'].items():
            path = Path(info['path']); bind(bindings, path, info['sha256']); report = json.loads(path.read_text())
            assert report['status'] == 'scalar_diagnostics_complete_not_posterior_qualification'
            actual = Counter(v['status'] for v in report['variables'].values()); assert dict(actual) == info['counts']
            diagnostic_counts.update({burnin + ':' + k: v for k, v in actual.items()})
            fully_passing[burnin] += all(v['status'] == 'passes_scalar_screen_only' for v in report['variables'].values())
            for path, digest in report['pins'].items(): bind(bindings, path, digest)
    assert complete == d['complete_quartets'] and unresolved == d['unresolved_quartets']
    verify(bindings)
    summary = dict(chains=1620, checked_chains=counts['all_saved_alignments_and_candidate_nodes_checked'],
        failed_chains=len(failed), complete_quartets=complete, unresolved_quartets=unresolved,
        scalar_variable_counts=dict(diagnostic_counts), quartets_passing_every_scalar=dict(fully_passing),
        chains_with_allocation_notices=len(allocation_notices))
    output = Path(plan['output']); output.mkdir(exist_ok=False)
    snapshot = output / 'accounting.json'
    snapshot.write_text(json.dumps(dict(status='verified_full_baliphy_initial_horizon_accounting_pending_journals',
        **summary, failed_chains_evidence=failed, allocation_notices=allocation_notices, source_hashes=bindings,
        scientific_eligibility=False, scope=plan['scope']), indent=2) + '\n')
    inner = output / 'closure_plan.json'; archive = output / 'completion_archive.json'
    spec = dict(output=str(archive), completed_status='complete_verified_baliphy_initial_horizon_accounting_archive',
        evidence=dict(accounting=dict(path=str(snapshot), expected=dict(status='verified_full_baliphy_initial_horizon_accounting_pending_journals', **summary))),
        links=[], launches=plan['launches'], pins={**bindings, str(snapshot): sha(snapshot)}, summary=summary, scope=plan['scope'])
    inner.write_text(json.dumps(spec, indent=2) + '\n')
    subprocess.run([sys.executable, 'scripts/record_completed_process_handoffs_v2.py', '--plan', str(inner)], check=True)
    proof = json.loads(archive.read_text()); assert len(proof['services']) == 2
    result = dict(status='complete_verified_baliphy_initial_horizon_accounting', **summary,
        accounting=str(snapshot), accounting_sha256=sha(snapshot), full_hash_archive=str(archive),
        full_hash_archive_sha256=sha(archive), bound_source_hashes=len(proof['source_hashes']),
        exact_process_journals_checked=2, scientific_eligibility=False, scope=plan['scope'])
    with Path(plan['completion']).open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__': main()
