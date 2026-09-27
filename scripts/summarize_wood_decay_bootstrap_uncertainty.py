#!/usr/bin/env python3
"""Summarize fully audited ecological edge ambiguity without treating it as origins."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import time
import subprocess
import psutil


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--wait-for', type=Path)
    args = parser.parse_args()
    script_hash = sha(__file__)
    if args.wait_for:
        identity = json.loads(args.wait_for.read_text())
        while True:
            try:
                process = psutil.Process(identity['pid'])
                if process.create_time() != identity['created'] or process.status() == psutil.STATUS_ZOMBIE:
                    break
                if process.cmdline() != identity['cmdline']:
                    raise ValueError('Audit identity changed')
            except psutil.NoSuchProcess:
                break
            print('waiting_for_ecology_bootstrap_audit', identity['pid'], flush=True)
            time.sleep(30)
    if args.wait_for:
        state = dict(line.split('=', 1) for line in subprocess.check_output(['systemctl', '--user', 'show', identity['unit'], '-p', 'ActiveState', '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
        assert state == {'ActiveState': 'inactive', 'Result': 'success', 'ExecMainStatus': '0'}, state
    bootstrap = Path('results/ecology/wood-decay-bootstrap-edges-20260927-v1')
    audit = Path('results/ecology/wood-decay-bootstrap-readback-20260927-v1/receipt.json')
    ml = Path('results/ecology/wood-decay-optimal-edge-states-20260927-v1')
    ml_audit = Path('metadata/wood_decay_optimal_edges_readback_20260927.json')
    proof, source, mlproof, mlsource = [json.loads(p.read_text()) for p in
                                     [audit, bootstrap / 'receipt.json', ml_audit, ml / 'receipt.json']]
    if (proof['status'] != 'passed_full_bootstrap_wood_decay_edge_network_flow_readback'
            or proof['producer_receipt_sha256'] != sha(bootstrap / 'receipt.json')
            or proof['trees'] != 2000 or proof['constrained_costs'] != 41960000):
        raise ValueError('Complete bootstrap network-flow audit required')
    if (mlproof['status'] != 'passed_full_wood_decay_optimal_edge_network_flow_readback'
            or mlproof['producer_receipt_sha256'] != sha(ml / 'receipt.json')):
        raise ValueError('Complete ML network-flow audit required')
    paths = [audit, bootstrap / 'receipt.json', ml_audit, ml / 'receipt.json']
    source_plan = Path('metadata/wood_decay_bootstrap_plan_20260927.json')
    assert source['plan_sha256'] == sha(source_plan)
    paths.append(source_plan)
    for root, receipt, names in [(bootstrap, source, ['split_frequencies.tsv', 'tree_summaries.tsv']),
                                 (ml, mlsource, ['edge_states.tsv'])]:
        for name in names:
            path = root / name
            if sha(path) != receipt['artifacts'][name]:
                raise ValueError('Changed input: ' + str(path))
            paths.append(path)
    hashes = {str(p): sha(p) for p in paths}
    frequency = read(bootstrap / 'split_frequencies.tsv')
    tree_rows = read(bootstrap / 'tree_summaries.tsv')
    ml_rows = read(ml / 'edge_states.tsv')
    key = lambda r: (r['tree_source'], r['coding'], r['split_id'])
    frequencies = {key(r): r for r in frequency}
    ml_edges = {key(r): r for r in ml_rows}
    if len(frequencies) != len(frequency) or len(ml_edges) != len(ml_rows):
        raise ValueError('Duplicate split keys')
    groups = defaultdict(list)
    universes = defaultdict(set)
    for row in tree_rows:
        groups[row['tree_source'], row['coding']].append(row)
    for row in ml_rows:
        universes[row['tree_source'], row['coding']].update(row['child_taxa'].split(';'))
    if len(groups) != 10 or set(universes) != set(groups):
        raise ValueError('Condition grid differs')
    distributions = []
    for group, rows in sorted(groups.items()):
        if len(rows) != 1000 or {int(r['tree_index']) for r in rows} != set(range(1, 1001)):
            raise ValueError('Incomplete bootstrap ensemble')
        for metric in ['minimum_changes', 'required_change', 'optional_change', 'no_change_in_any_optimum']:
            for value, count in sorted(Counter(int(r[metric]) for r in rows).items()):
                distributions.append(dict(tree_source=group[0], coding=group[1], metric=metric,
                                          value=value, trees=count, ensemble_trees=1000, fraction=count / 1000))
    paired = defaultdict(dict)
    for row in tree_rows:
        k = (row['tree_source'], int(row['tree_index']))
        assert row['coding'] not in paired[k]
        paired[k][row['coding']] = row
    coding_rows = []
    for (label, index), scenarios in sorted(paired.items()):
        assert len(scenarios) == 5 and 'uncertain_unknown' in scenarios
        for metric in ['minimum_changes', 'required_change', 'optional_change', 'no_change_in_any_optimum']:
            primary = int(scenarios['uncertain_unknown'][metric])
            values = [int(r[metric]) for name, r in scenarios.items() if name != 'uncertain_unknown']
            coding_rows.append(dict(tree_source=label, tree_index=index, metric=metric,
                                    primary_unknown=primary, assignment_min=min(values), assignment_max=max(values),
                                    minimum_difference=min(values)-primary, maximum_difference=max(values)-primary))
    assert len(paired) == 2000 and len(coding_rows) == 8000
    output = []
    for split_key in sorted(set(frequencies) | set(ml_edges)):
        f, m = frequencies.get(split_key), ml_edges.get(split_key)
        group = split_key[:2]
        if m:
            side = set(m['child_taxa'].split(';'))
            canonical = min([sorted(side), sorted(universes[group] - side)], key=lambda s: (len(s), s))
            if hashlib.sha256(json.dumps(canonical, separators=(',', ':')).encode()).hexdigest() != split_key[2]:
                raise ValueError('ML split identity differs')
            canonical = ';'.join(canonical)
            if f and f['canonical_side'] != canonical:
                raise ValueError('ML/bootstrap split side differs')
        else:
            canonical = f['canonical_side']
        counts = {name: int(f[name]) if f else 0 for name in
                  ['present', 'required_change', 'optional_change', 'no_change_in_any_optimum']}
        if f and int(f['ensemble_trees']) != 1000:
            raise ValueError('Unexpected ensemble denominator')
        if sum(counts[k] for k in ['required_change', 'optional_change', 'no_change_in_any_optimum']) != counts['present']:
            raise ValueError('Edge status partition differs')
        row = dict(tree_source=group[0], coding=group[1], split_id=split_key[2], canonical_side=canonical,
                   ml_status=m['status'] if m else 'absent_from_ml_tree', ensemble_trees=1000,
                   absent_trees=1000-counts['present'], **counts)
        row['presence_fraction'] = counts['present'] / 1000
        for status in ['required_change', 'optional_change', 'no_change_in_any_optimum']:
            row[status + '_fraction_all_trees'] = counts[status] / 1000
            row[status + '_fraction_when_present'] = counts[status] / counts['present'] if counts['present'] else ''
        output.append(row)
    for path, expected in hashes.items():
        if sha(path) != expected:
            raise ValueError('Source changed during summary')
    if sha(__file__) != script_hash:
        raise ValueError('Summary script changed during run')
    args.output.mkdir(parents=True, exist_ok=False)
    for name, rows in [('split_uncertainty.tsv', output), ('tree_metric_distributions.tsv', distributions), ('paired_coding_sensitivity.tsv', coding_rows)]:
        with (args.output / name).open('x') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader(); writer.writerows(rows)
    receipt = dict(status='complete_audited_wood_decay_bootstrap_uncertainty_summary_pending_readback',
                   sources=hashes, script_sha256=script_hash, split_rows=len(output), distribution_rows=len(distributions), paired_coding_rows=len(coding_rows),
                   trees=2000, coding_tree_combinations=10000,
                   artifacts={p.name: sha(p) for p in args.output.iterdir()},
                   scope='Descriptive root-free equal-cost state-mapping sensitivity. Split absence is separate from unchanged states; conditional and whole-ensemble denominators are explicit. Required/optional statuses concern optimal mappings, not identified biological origins, independent contrasts, posterior probabilities or ecological effects. The primary unknown coding and four explicit uncertain-taxon assignments remain separate; paired ranges are sensitivity bounds, not confidence intervals.')
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: receipt[k] for k in ['status', 'split_rows', 'distribution_rows']}), flush=True)


if __name__ == '__main__':
    main()
