#!/usr/bin/env python3
"""Audit completed case arithmetic and summarize start disagreement, retaining pending cases."""
import argparse
import csv
import json
import math
import re
from pathlib import Path
from ancestral_chain_attempt import sha, write_json


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    pp = Path('metadata/local_mg94_unconstrained_multistart_plan_20260927.json')
    plan = json.loads(pp.read_text())
    for path, digest in plan['pins'].items():
        assert sha(path) == digest
    manifest = Path(plan['starts']) / 'start_manifest.tsv'
    starts = {(r['case_id'], r['start_label']): r for r in csv.DictReader(manifest.open(), delimiter='\t')}
    cases = sorted({key[0] for key in starts})
    assert len(cases) == 1632 and len(starts) == 13056
    pins = {str(pp): sha(pp), str(manifest): sha(manifest), str(Path(__file__)): sha(__file__)}
    rows = []
    root = Path(plan['output'])
    for case in cases:
        folder = root / case
        rp = folder / 'receipt.json'
        if not rp.exists():
            rows.append(dict(case_id=case, status='pending_complete_case'))
            continue
        pins[str(rp)] = sha(rp)
        receipt = json.loads(rp.read_text())
        assert receipt['case_id'] == case and receipt['status'] == 'complete_all_eight_unconstrained_starts_with_fresh_readbacks'
        fits = receipt['rows']
        assert len(fits) == len(receipt['proofs']) == 8
        assert {r['start_label'] for r in fits} == {k[1] for k in starts if k[0] == case}
        for proof in receipt['proofs']:
            for name, digest in proof['artifacts'].items():
                assert sha(folder / proof['start_label'] / name) == digest
        best = max(r['log_likelihood'] for r in fits)
        worst = min(r['log_likelihood'] for r in fits)
        for fit in fits:
            start = starts[case, fit['start_label']]
            assert all(math.isfinite(fit[k]) for k in ['log_likelihood', 'starting_log_likelihood', 'omega', 'target_parameter'])
            assert abs(fit['starting_log_likelihood'] - float(start['starting_log_likelihood'])) <= 1e-6
            assert fit['log_likelihood'] >= fit['starting_log_likelihood'] - 1e-5
            assert abs(fit['deficit_from_best_start'] - (best - fit['log_likelihood'])) <= 1e-9
            work = folder / fit['start_label']
            replay = re.findall(r'^READBACK=(\S+)', (work / 'readback.log').read_text(), re.M)
            assert len(replay) == 1 and abs(float(replay[0]) - fit['log_likelihood']) <= 1e-6
        near = [r for r in fits if best - r['log_likelihood'] <= 1e-5]
        rows.append(dict(case_id=case, status='completed_case_hashes_and_score_arithmetic_checked',
            between_start_log_likelihood_spread=best - worst,
            starts_within_1e_5_of_best=len(near),
            omega_min=min(r['omega'] for r in fits), omega_max=max(r['omega'] for r in fits),
            near_best_omega_min=min(r['omega'] for r in near), near_best_omega_max=max(r['omega'] for r in near)))
    args.output.mkdir(parents=True, exist_ok=False)
    write_json(args.output / 'case_dispositions.json', rows)
    complete = [r for r in rows if r['status'].startswith('completed_')]
    summary = dict(status='partial_codon_start_disagreement_snapshot', cases=1632,
        completed_cases=len(complete), pending_cases=len(rows) - len(complete),
        spreads_above={str(cutoff): sum(r['between_start_log_likelihood_spread'] > cutoff for r in complete)
                      for cutoff in [1e-5, .01, 1, 10]},
        maximum_spread=max((r['between_start_log_likelihood_spread'] for r in complete), default=None),
        near_best_omega_range_above_0_1=sum(r['near_best_omega_max']-r['near_best_omega_min'] > .1 for r in complete),
        pins=pins, artifacts={'case_dispositions.json': sha(args.output / 'case_dispositions.json')},
        scope='Frozen partial inventory and completed artifact/score checks. No full free-parameter/source-model audit, '
              'global optimum, dS interval, selection eligibility, or biological selection claim.')
    write_json(args.output / 'receipt.json', summary)
    print(json.dumps({k: v for k, v in summary.items() if k not in ['pins', 'artifacts']}))


if __name__ == '__main__':
    main()
