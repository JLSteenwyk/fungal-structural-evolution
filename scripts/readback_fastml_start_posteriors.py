#!/usr/bin/env python3
"""Recompute every start-sensitivity summary with scalar posterior accumulators."""
import argparse
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path
import subprocess
from Bio import Phylo
from screen_duplication_alignment_reuse import sha


def update_extrema(accumulator, key, value):
    if key not in accumulator:
        accumulator[key] = [value, value]
    else:
        lo, hi = accumulator[key]
        accumulator[key] = [min(lo, value), max(hi, value)]


def range_summary(accumulator):
    assert accumulator
    ranges = [hi-lo for lo, hi in accumulator.values()]
    return dict(maximum_internal_probability_range=max(ranges),
                internal_rows_range_above_0_01=sum(v > .01 for v in ranges),
                internal_rows_range_above_0_1=sum(v > .1 for v in ranges))


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--plan', required=True, type=Path)
    args = parser.parse_args(); config = json.loads(args.plan.read_text())
    plan = json.loads(Path(config['source_plan']).read_text())
    root = Path(plan['output']); rp = root/'receipt.json'; receipt = json.loads(rp.read_text())
    assert receipt['status'] == 'complete_full_fastml_start_posterior_sensitivity'
    bindings = {str(args.plan):sha(args.plan), **config['pins'], **receipt['source_hashes'], str(rp):sha(rp)}
    bindings.update({str(root/name):h for name,h in receipt['artifacts'].items()})
    states = {}
    for path in [config['source_launch'], *plan['completed_launches']]:
        launch = json.loads(Path(path).read_text())
        if 'plan' in launch: assert sha(launch['plan']) == launch['plan_sha256']
        state = dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',launch['unit'],
            '-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
        assert state == dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
        states[launch['unit']] = state
    def verify():
        for path,h in bindings.items(): assert sha(path) == h,path
    verify()
    source = json.loads((Path(plan['producer'])/'receipt.json').read_text())
    groups = defaultdict(list)
    for entry in source['dispositions'].values():
        assert sha(entry['path']) == entry['sha256']
        row = json.loads(Path(entry['path']).read_text())
        assert row['status'] == entry['status']
        groups[row['job']['job_id']].append(row)
    assert len(groups) == 156 and sum(map(len,groups.values())) == 780
    expected = {}
    for line in (root/'start_sensitivity.jsonl').open():
        row = json.loads(line); assert row['job_id'] not in expected
        expected[row['job_id']] = row
    assert set(expected) == set(groups)
    summaries = []; total_rows = 0; probability_files = 0; near_counts = Counter(); reported_excursions = 0
    tolerance = plan['near_best_log_likelihood_tolerance']
    for job, starts in sorted(groups.items()):
        assert len(starts) == 5 and len({r['start'] for r in starts}) == 5
        assert all(r['job'] == starts[0]['job'] for r in starts)
        characters = starts[0]['job']['character_count']
        if characters == 0:
            assert all(r['status'] == 'no_coded_characters' for r in starts)
            result = dict(job_id=job,status='no_coded_characters')
            assert expected[job] == result; summaries.append(result); continue
        assert all(r['status'] == 'refined_native_settings_and_replay_checked_requires_optimization_review' for r in starts)
        assert all(r['raw_values_clipped'] is False for r in starts)
        assert all(math.isfinite(r['native_rate_replay_log_likelihood']) for r in starts)
        best = min(starts,key=lambda r:(-r['native_rate_replay_log_likelihood'],r['start']))
        best_ll = best['native_rate_replay_log_likelihood']
        near = [r for r in starts if best_ll-r['native_rate_replay_log_likelihood'] <= tolerance]
        near_labels = {r['start'] for r in near}; near_counts[str(len(near))] += 1
        extrema = {'all':{},'near_best':{}}; node_set = tip_set = None
        # Traverse by start label, rather than the producer's likelihood-sorted
        # matrix row order. Explicit keys align every observation across starts.
        for row in sorted(starts,key=lambda r:r['start']):
            attempt_path = Path(row['attempt_receipt'])
            assert sha(attempt_path) == row['attempt_receipt_sha256']
            attempt = json.loads(attempt_path.read_text())
            posterior = attempt_path.parent/'RESULTS/AncestralReconstructPosterior.txt'
            tree_path = attempt_path.parent/'RESULTS/TheTree.INodes.ph'
            for path,name in [(posterior,'RESULTS/AncestralReconstructPosterior.txt'),(tree_path,'RESULTS/TheTree.INodes.ph')]:
                assert sha(path) == attempt['artifacts'][name] == bindings[str(path)]
            tree = Phylo.read(tree_path,'newick'); clades = list(tree.find_clades())
            nodes = {n.name or n.comment for n in clades}; tips = {n.name for n in tree.get_terminals()}
            assert None not in nodes and len(nodes) == len(clades) and tips <= nodes
            if node_set is None: node_set,tip_set = nodes,tips
            else: assert nodes == node_set and tips == tip_set
            observed = set()
            with posterior.open() as f:
                reader = csv.reader(f,delimiter='\t'); assert next(reader) == ['POS','Node','State','Prob']
                for fields in reader:
                    assert len(fields) == 4
                    position,node,state,probability = fields
                    position = int(position); probability = float(probability)
                    key = position,node
                    assert key not in observed and 1 <= position <= characters and node in nodes and state == '1'
                    observed.add(key)
                    assert math.isfinite(probability)
                    boundary = row['probability_boundary_tolerance']
                    assert -boundary <= probability <= 1+boundary
                    reported_excursions += probability < 0 or probability > 1
                    if node not in tips:
                        update_extrema(extrema['all'],key,probability)
                        if row['start'] in near_labels: update_extrema(extrema['near_best'],key,probability)
            # Unique keys within valid position/node bounds and the complete
            # cardinality prove the full Cartesian grid without materializing it.
            assert len(observed) == characters*len(nodes) == row['probability_rows']
            total_rows += len(observed); probability_files += 1
        internal_count = characters*(len(node_set)-len(tip_set))
        result = dict(job_id=job,status='descriptive_start_sensitivity_not_qualified',best_start=best['start'],
                      start_log_likelihood_range=max(r['native_rate_replay_log_likelihood'] for r in starts)-min(r['native_rate_replay_log_likelihood'] for r in starts),
                      near_best_starts=len(near),posterior_rows_per_start=characters*len(node_set),
                      internal_rows_per_start=internal_count,
                      starts_with_iteration_limit=sum(r['optimizer_settings']['model_iteration_limit_messages']>0 for r in starts),
                      best_has_iteration_limit=best['optimizer_settings']['model_iteration_limit_messages']>0)
        for label, chosen in [('all',starts),('near_best',near)]:
            assert len(extrema[label]) == internal_count
            result.update({label+'_'+k:v for k,v in range_summary(extrema[label]).items()})
            for parameter in ['alpha','gain','loss']:
                values = [r['fitted_parameters'][parameter] for r in chosen]
                assert all(math.isfinite(v) for v in values)
                result[label+'_'+parameter+'_minimum'] = min(values)
                result[label+'_'+parameter+'_maximum'] = max(values)
        assert result == expected[job],job
        summaries.append(result)
        print('Independently checked posterior start ranges',len(summaries),'/156',flush=True)
    nonempty = [r for r in summaries if r['status'] != 'no_coded_characters']
    assert len(nonempty) == receipt['nonempty_groups'] == 153
    assert total_rows == receipt['posterior_rows_read'] == 40291700 and probability_files == 765
    for key, value in [('groups',len(summaries)),('near_best_tolerance',tolerance),
        ('groups_with_likelihood_range_above_tolerance',sum(r['start_log_likelihood_range']>tolerance for r in nonempty)),
        ('groups_with_best_iteration_limit',sum(r['best_has_iteration_limit'] for r in nonempty))]:
        assert receipt[key] == value
    for label in ['all','near_best']:
        assert receipt[label+'_groups_with_internal_range_above_0_01'] == sum(r[label+'_internal_rows_range_above_0_01']>0 for r in nonempty)
        assert receipt[label+'_maximum_internal_probability_range'] == max(r[label+'_maximum_internal_probability_range'] for r in nonempty)
    verify()
    result = dict(status='passed_full_fastml_start_posterior_sensitivity_readback',groups=len(summaries),nonempty_groups=len(nonempty),
                  posterior_files=probability_files,posterior_rows_checked=total_rows,near_best_start_counts=dict(near_counts),
                  raw_probability_excursions_retained=reported_excursions,source_receipt_sha256=sha(rp),
                  source_hashes={str(p):sha(p) for p in [args.plan,Path(config['source_plan']),Path(config['source_launch']),rp,
                      root/'start_sensitivity.jsonl',Path(__file__)]},terminal_states=states,
                  scope='Every raw posterior row and exact position/node grid checked. Every saved group range/count, parameter range, start choice, likelihood range and aggregate independently recomputed with scalar keyed min/max rather than array stacking/reductions. Tree parser and source files shared; no optimizer or pruning rerun. Singleton near-best spread is automatic. No global-optimum, calibrated indel model or ancestral-ensemble qualification.')
    with Path(config['output']).open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','terminal_states']}),flush=True)


if __name__ == '__main__': main()
