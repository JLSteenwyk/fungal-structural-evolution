#!/usr/bin/env python3
"""Measure ancestral indel posterior sensitivity across all completed starts."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
import subprocess

from Bio import Phylo
import numpy as np
from screen_duplication_alignment_reuse import sha


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plan', required=True, type=Path)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    pins = {str(args.plan): sha(args.plan), **plan['pins']}
    for launch_path in plan['completed_launches']:
        launch = json.loads(Path(launch_path).read_text())
        state = dict(x.split('=', 1) for x in subprocess.check_output([
            'systemctl', '--user', 'show', launch['unit'], '-p', 'ActiveState',
            '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
        assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0')
    root, audit = Path(plan['producer']), Path(plan['audit'])
    receipt = json.loads((root / 'receipt.json').read_text())
    checked = json.loads((audit / 'receipt.json').read_text())
    assert checked['status'] == 'full_refinement_output_readback_and_start_comparison_complete'
    assert checked['source_receipt_sha256'] == sha(root / 'receipt.json')
    for name, h in checked['artifacts'].items():
        pins[str(audit / name)] = h
    groups = defaultdict(list)
    for entry in receipt['dispositions'].values():
        pins[entry['path']] = entry['sha256']
        assert sha(entry['path']) == entry['sha256']
        row = json.loads(Path(entry['path']).read_text())
        assert row['status'] == entry['status']
        groups[row['job']['job_id']].append(row)
    assert len(groups) == 156 and sum(map(len, groups.values())) == 780
    def verify():
        for p, h in pins.items():
            assert sha(p) == h, p
    verify()
    output = Path(plan['output'])
    output.mkdir(parents=True, exist_ok=False)
    summaries = []
    total_rows = 0
    for job, rows in sorted(groups.items()):
        assert len(rows) == 5 and len({r['start'] for r in rows}) == 5
        if rows[0]['job']['character_count'] == 0:
            assert all(r['status'] == 'no_coded_characters' for r in rows)
            summaries.append(dict(job_id=job, status='no_coded_characters'))
            continue
        assert all(r['status'] == 'refined_native_settings_and_replay_checked_requires_optimization_review' for r in rows)
        rows.sort(key=lambda r: (-r['native_rate_replay_log_likelihood'], r['start']))
        likelihoods = np.array([r['native_rate_replay_log_likelihood'] for r in rows])
        posteriors, keys, internal = [], None, None
        for row in rows:
            rp = Path(row['attempt_receipt'])
            pins[str(rp)] = row['attempt_receipt_sha256']
            assert sha(rp) == pins[str(rp)]
            attempt = json.loads(rp.read_text())
            paths = ['RESULTS/AncestralReconstructPosterior.txt', 'RESULTS/TheTree.INodes.ph']
            for name in paths:
                pins[str(rp.parent / name)] = attempt['artifacts'][name]
                assert sha(rp.parent / name) == pins[str(rp.parent / name)]
            tree = Phylo.read(rp.parent / paths[1], 'newick')
            tips = {n.name for n in tree.get_terminals()}
            nodes = {n.name or n.comment for n in tree.find_clades()}
            values = {}
            with (rp.parent / paths[0]).open() as handle:
                for r in csv.DictReader(handle, delimiter='\t'):
                    key = (int(r['POS']), r['Node'], r['State'])
                    assert key not in values and key[1] in nodes and key[2] == '1'
                    assert 1 <= key[0] <= row['job']['character_count']
                    values[key] = float(r['Prob'])
            assert len(values) == row['probability_rows'] == len(nodes)*row['job']['character_count']
            ordered = sorted(values)
            mask = np.array([k[1] not in tips for k in ordered])
            if keys is None:
                keys, internal = ordered, mask
            assert keys == ordered and np.array_equal(internal, mask)
            vector = np.array([values[k] for k in ordered])
            assert np.isfinite(vector).all()
            tolerance = row['probability_boundary_tolerance']
            assert vector.min() >= -tolerance and vector.max() <= 1+tolerance
            posteriors.append(vector)
            total_rows += len(vector)
        matrix = np.stack(posteriors)
        near = likelihoods[0]-likelihoods <= plan['near_best_log_likelihood_tolerance']
        result = dict(job_id=job, status='descriptive_start_sensitivity_not_qualified',
                      best_start=rows[0]['start'], start_log_likelihood_range=float(np.ptp(likelihoods)),
                      near_best_starts=int(near.sum()), posterior_rows_per_start=len(keys),
                      internal_rows_per_start=int(internal.sum()),
                      starts_with_iteration_limit=sum(r['optimizer_settings']['model_iteration_limit_messages']>0 for r in rows),
                      best_has_iteration_limit=rows[0]['optimizer_settings']['model_iteration_limit_messages']>0)
        for label, selected in [('all', np.ones(5, dtype=bool)), ('near_best', near)]:
            ranges = np.ptp(matrix[selected][:, internal], axis=0)
            result[label+'_maximum_internal_probability_range'] = float(ranges.max())
            result[label+'_internal_rows_range_above_0_01'] = int((ranges > .01).sum())
            result[label+'_internal_rows_range_above_0_1'] = int((ranges > .1).sum())
            for parameter in ['alpha', 'gain', 'loss']:
                values = [r['fitted_parameters'][parameter] for r, use in zip(rows, selected) if use]
                result[label+'_'+parameter+'_minimum'] = min(values)
                result[label+'_'+parameter+'_maximum'] = max(values)
        summaries.append(result)
        print('Compared start posteriors', len(summaries), '/156', flush=True)
    verify()
    target = output / 'start_sensitivity.jsonl'
    target.write_text(''.join(json.dumps(r, sort_keys=True, allow_nan=False)+'\n' for r in summaries))
    nonempty = [r for r in summaries if r['status'] != 'no_coded_characters']
    result = dict(status='complete_full_fastml_start_posterior_sensitivity', groups=len(summaries),
                  nonempty_groups=len(nonempty), posterior_rows_read=total_rows,
                  groups_with_likelihood_range_above_tolerance=sum(r['start_log_likelihood_range']>plan['near_best_log_likelihood_tolerance'] for r in nonempty),
                  groups_with_best_iteration_limit=sum(r['best_has_iteration_limit'] for r in nonempty),
                  near_best_tolerance=plan['near_best_log_likelihood_tolerance'],
                  source_hashes=pins, artifacts={target.name: sha(target)},
                  scope='Full five-start internal-node posterior and parameter ranges, preserving raw values and empty inputs. Near-best is a descriptive likelihood tolerance, not a confidence region. A single near-best start yields zero range without stability evidence. No global optimum, calibrated uncertainty, biological indel model or ancestral-sequence qualification.')
    for label in ['all', 'near_best']:
        result[label+'_groups_with_internal_range_above_0_01'] = sum(r[label+'_internal_rows_range_above_0_01']>0 for r in nonempty)
        result[label+'_maximum_internal_probability_range'] = max(r[label+'_maximum_internal_probability_range'] for r in nonempty)
    (output / 'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes', 'artifacts']}))


if __name__ == '__main__':
    main()
