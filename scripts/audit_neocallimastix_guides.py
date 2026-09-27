#!/usr/bin/env python3
"""Wait for six representative-sensitivity guides, audit and compare common tips."""
import argparse
import copy
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess
import time

from Bio import Phylo, SeqIO
import psutil


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def splits(tree, common):
    result = set()
    for clade in tree.find_clades():
        side = {tip.name for tip in clade.get_terminals()} & common
        other = common - side
        if min(len(side), len(other)) > 1:
            a, b = tuple(sorted(side)), tuple(sorted(other))
            result.add(min((a, b), key=lambda x: (len(x), x)))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    digest = sha(args.plan)
    def verify():
        assert sha(args.plan) == digest
        for path, expected in plan['pins'].items():
            assert sha(path) == expected, path
    verify()
    output = Path(plan['output'])
    output.mkdir(parents=True, exist_ok=False)
    producer = plan['producer']
    while True:
        state = dict(line.split('=', 1) for line in subprocess.check_output(
            ['systemctl', '--user', 'show', producer['unit'], '--property=ActiveState,SubState,Result,ExecMainStatus,MainPID'], text=True).splitlines())
        if state['ActiveState'] == 'inactive':
            assert state['Result'] == 'success' and state['ExecMainStatus'] == '0', state
            break
        assert state['ActiveState'] in ['active', 'activating', 'deactivating'], state
        if state['SubState'] == 'running':
            assert int(state['MainPID']) == producer['pid']
            process = psutil.Process(producer['pid'])
            assert process.create_time() == producer['created'] and process.cmdline() == producer['cmdline']
        print('waiting_for_six_neocallimastix_guides', producer['pid'], flush=True)
        time.sleep(30)
    verify()
    inputs = Path(plan['inputs'])
    runs = Path(plan['runs'])
    source = json.loads((inputs / 'receipt.json').read_text())
    config = json.loads((runs / 'config.json').read_text())
    result = json.loads((runs / 'receipt.json').read_text())
    assert result['status'] == 'complete_neocallimastix_guide_sensitivities'
    assert result['config_sha256'] == sha(runs / 'config.json')
    assert config['source_receipt_sha256'] == sha(inputs / 'receipt.json')
    assert config['model'] == 'LG+F+G4' and config['seed'] == 20260913
    assert len(source['matrices']) == len(result['results']) == 6
    expected_names = {m['alignment'] + '-retain_' + m['representative'] for m in source['matrices']}
    assert {r['name'] for r in result['results']} == expected_names
    with Path('metadata/analysis_manifest.tsv').open() as handle:
        manifest = {r['taxon_id']: r for r in csv.DictReader(handle, delimiter='\t')}
    trees, reports, bindings = {}, [], {}
    for item in source['matrices']:
        name = item['alignment'] + '-retain_' + item['representative']
        matrix = inputs / item['path']
        assert sha(matrix) == item['sha256']
        sequences = list(SeqIO.parse(matrix, 'fasta'))
        expected = set(manifest) - (set(source['group_taxa']) - {item['representative']})
        assert len(sequences) == len(expected) == 524
        assert {s.id for s in sequences} == expected
        assert all(len(s.seq) == item['columns'] for s in sequences)
        assert sum(manifest[t]['study_role'] == 'outgroup' for t in expected) == 25
        folder = runs / name
        rp = folder / 'receipt.json'
        receipt = json.loads(rp.read_text())
        assert receipt == next(r for r in result['results'] if r['name'] == name)
        assert receipt['status'] == 'complete_sensitivity_guide'
        assert receipt['config_sha256'] == sha(runs / 'config.json')
        assert receipt['input_sha256'] == sha(matrix)
        command = receipt['command']
        assert sha(command[0]) == config['executable_sha256']
        for flag, value in [('-s', str(matrix.resolve())), ('-st', 'AA'), ('-m', 'LG+F+G4'), ('-T', '8'), ('--mem', '32G'), ('--seed', '20260913')]:
            assert command[command.index(flag) + 1] == value
        for filename, expected_hash in receipt['artifacts'].items():
            assert sha(folder / filename) == expected_hash
        tree = Phylo.read(folder / 'guide.treefile', 'newick')
        tips = [t.name for t in tree.get_terminals()]
        assert len(tips) == len(set(tips)) == 524 and set(tips) == expected
        branches = [c.branch_length for c in tree.find_clades() if c is not tree.root]
        assert all(x is not None and math.isfinite(x) and x >= 0 for x in branches)
        trees[name] = tree
        bindings[str(rp)] = sha(rp)
        reports.append(dict(name=name, taxa=len(tips), outgroups=25, columns=item['columns'], internal_splits=len(splits(tree, expected)), total_branch_length=sum(branches)))
    comparisons = []
    for left, right in itertools.combinations(sorted(trees), 2):
        common = {t.name for t in trees[left].get_terminals()} & {t.name for t in trees[right].get_terminals()}
        splitsets = []
        for name in [left, right]:
            projected = splits(trees[name], common)
            pruned = copy.deepcopy(trees[name])
            for tip in list(pruned.get_terminals()):
                if tip.name not in common:
                    pruned.prune(tip)
            assert splits(pruned, common) == projected
            splitsets.append(projected)
        a, b = splitsets
        rf = len(a ^ b)
        comparisons.append(dict(left=left, right=right, common_taxa=len(common), left_splits=len(a), right_splits=len(b), shared_splits=len(a & b), rf_distance=rf, rf_over_split_sum=rf/(len(a)+len(b)) if a or b else 0.))
    assert len(comparisons) == 15
    for filename, rows in [('tree_audit.tsv', reports), ('pairwise_common_tip_topology.tsv', comparisons)]:
        with (output / filename).open('w') as handle:
            writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)
    verify()
    proof = dict(status='complete_six_neocallimastix_guide_audits_and_topology_comparisons', plan_sha256=digest, terminal_state=state, producer_receipt_sha256=sha(runs/'receipt.json'), tree_receipts=bindings, trees=6, pairwise_comparisons=15, artifacts={p.name:sha(p) for p in output.iterdir()}, scope='All six exact-tip guide trees and artifacts verified; unrooted splits compared on common taxa, cross-checked by explicit pruning. No likelihood recomputation, bootstrap support, species delimitation or final species phylogeny claim. Different representative comparisons omit all three disputed taxa from the shared-tip topology.')
    (output / 'receipt.json').write_text(json.dumps(proof, indent=2) + '\n')
    print(json.dumps(proof), flush=True)


if __name__ == '__main__':
    main()
