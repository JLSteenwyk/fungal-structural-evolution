#!/usr/bin/env python3
"""Read completed genus trees and independently recount all bootstrap splits."""
import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import re
from Bio import Phylo, SeqIO


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def split_map(tree, taxa):
    tips = [tip.name for tip in tree.get_terminals()]
    if len(tips) != len(taxa) or set(tips) != taxa:
        raise ValueError('Missing, duplicate or extra tree tips')
    splits = {}
    for node in tree.find_clades():
        if node is tree.root:
            continue
        if node.branch_length is None or not math.isfinite(node.branch_length) or node.branch_length < 0:
            raise ValueError('Invalid tree edge')
        if node.is_terminal():
            continue
        side = frozenset(t.name for t in node.get_terminals())
        key = min((tuple(sorted(side)), tuple(sorted(taxa - side))), key=lambda x: (len(x), x))
        if len(key) < 2 or key in splits:
            raise ValueError('Trivial or duplicated internal split')
        splits[key] = node
    if len(splits) != len(taxa) - 3:
        raise ValueError('Expected resolved unrooted topology')
    return splits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trees', required=True, type=Path)
    parser.add_argument('--inputs', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--allow-incomplete', action='store_true')
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new audit directory')
    config_path = args.trees / 'config.json'
    config = json.loads(config_path.read_text())
    if sha(args.inputs / 'receipt.json') != config['input_receipt_sha256']:
        raise ValueError('Changed input provenance')
    info_path = args.trees / 'information.tsv'
    if sha(info_path) != config['information_sha256']:
        raise ValueError('Changed information screen')
    expected = {r['case_id']: r for r in table(info_path) if r['status'] == 'ready_for_supported_tree_diagnostic'}
    completed = {p.parent.name: p for p in args.trees.glob('*/receipt.json')}
    if not completed or set(completed) - set(expected):
        raise ValueError('No completed cases or unexpected cases')
    pending = sorted(set(expected) - set(completed))
    if pending and not args.allow_incomplete:
        raise ValueError('Full tree execution is not complete')
    if not pending:
        batch = json.loads((args.trees / 'receipt.json').read_text())
        if batch['status'] != 'complete_genus_nucleotide_tree_execution_pending_full_audit' or batch['config_sha256'] != sha(config_path):
            raise ValueError('Invalid full execution receipt')
        if {r['case_id']: r['receipt_sha256'] for r in batch['case_receipts']} != {case: sha(p) for case, p in completed.items()}:
            raise ValueError('Full receipt grid differs')
    summaries, branches, warnings, sources = [], [], [], []
    for case, receipt_path in sorted(completed.items()):
        folder = receipt_path.parent
        receipt = json.loads(receipt_path.read_text())
        rc_path = folder / 'config.json'
        rc = json.loads(rc_path.read_text())
        if receipt['status'] != 'complete_supported_nucleotide_tree_pending_full_audit' or receipt['case_id'] != case or receipt['config_sha256'] != sha(rc_path) or rc['parent_config_sha256'] != sha(config_path):
            raise ValueError('Changed case provenance')
        for name, digest in receipt['artifacts'].items():
            if sha(folder / name) != digest:
                raise ValueError('Changed tree artifact')
        alignment = args.inputs / case / 'codons.fna'
        if sha(alignment) != rc['alignment_sha256']:
            raise ValueError('Changed alignment')
        records = list(SeqIO.parse(alignment, 'fasta'))
        taxa = {r.id for r in records}
        ncols = int(expected[case]['nucleotide_columns'])
        if len(taxa) != len(records) or len(taxa) != int(expected[case]['taxa']) or any(len(r.seq) != ncols for r in records):
            raise ValueError('Alignment grid differs')
        command = rc['command']
        for option, value in {'-st':'DNA', '-m':'GTR+F+G4', '--alrt':'1000', '-B':'1000'}.items():
            if command[command.index(option)+1] != value:
                raise ValueError('Wrong inference settings')
        if not {'--bnni','--boot-trees','-keep-ident'} <= set(command):
            raise ValueError('Missing inference flags')
        report = (folder / 'tree.iqtree').read_text()
        required = [f'Input data: {len(taxa)} sequences with {ncols} nucleotide sites',
                    'Model of substitution: GTR+F+G4',
                    'tree reconstruction + ultrafast bootstrap (1000 replicates)',
                    'SH-aLRT support (%) / ultrafast bootstrap support (%)']
        if not all(text in report for text in required):
            raise ValueError('Report settings differ: '+case)
        likelihood = float(re.search(r'Log-likelihood of the tree: (\S+)', report)[1])
        alpha = float(re.search(r'Gamma shape alpha: (\S+)', report)[1])
        if not math.isfinite(likelihood) or not math.isfinite(alpha) or alpha <= 0:
            raise ValueError('Nonfinite fit or invalid gamma shape')
        tree = Phylo.read(folder / 'tree.treefile', 'newick')
        ml = split_map(tree, taxa)
        counts = Counter()
        boot_count = 0
        for boot in Phylo.parse(folder / 'tree.ufboot', 'newick'):
            counts.update(split_map(boot, taxa).keys())
            boot_count += 1
        if boot_count != 1000:
            raise ValueError('Wrong bootstrap replicate count')
        largest_difference = 0
        for side, node in ml.items():
            parts = str(node.name).split('/')
            if len(parts) != 2:
                raise ValueError('Missing paired branch supports')
            alrt, ufb = map(float, parts)
            if any(not math.isfinite(s) or not 0 <= s <= 100 for s in [alrt, ufb]):
                raise ValueError('Invalid support range')
            observed = counts[side] / 10
            difference = abs(observed - ufb)
            largest_difference = max(largest_difference, difference)
            branches.append({'case_id':case,'smaller_side_taxa':';'.join(side),'length':node.branch_length,
                             'sh_alrt_percent':alrt,'reported_ufb_percent':ufb,'recounted_ufb_percent':observed,
                             'ufb_absolute_difference':difference})
        # Integer support labels can differ from empirical percentages by rounding.
        if largest_difference > 0.500001:
            raise ValueError('Bootstrap support differs beyond rounding: '+case)
        log = (folder / 'tree.log').read_text()
        case_warnings = [line.strip() for line in log.splitlines() if 'WARNING' in line.upper()]
        warnings.extend({'case_id':case,'warning':line} for line in case_warnings)
        edge_lengths = [n.branch_length for n in tree.find_clades() if n is not tree.root]
        total = sum(edge_lengths)
        reported_total = float(re.search(r'Total tree length \(sum of branch lengths\): (\S+)', report)[1])
        if abs(total - reported_total) > 0.000051:
            raise ValueError('Reported total length differs')
        summaries.append({'case_id':case,'taxa':len(taxa),'nucleotide_columns':ncols,'log_likelihood':likelihood,
                          'gamma_alpha':alpha,'total_tree_length':total,'internal_edges':len(ml),
                          'bootstrap_trees':boot_count,'maximum_ufb_rounding_difference':largest_difference,
                          'near_zero_edges_le1e_5':sum(x <= 1e-5 for x in edge_lengths),
                          'long_edges_ge10':sum(x >= 10 for x in edge_lengths),'warning_lines':len(case_warnings),
                          'translation_table':rc['translation_table'],'marker_copy_caveat':rc['marker_copy_caveat']})
        sources.append({'case_id':case,'receipt_sha256':sha(receipt_path)})
    args.output.mkdir(parents=True)
    for name, rows, fields in [('case_audit.tsv',summaries,list(summaries[0])),('branch_support.tsv',branches,list(branches[0])),('warnings.tsv',warnings,['case_id','warning'])]:
        with (args.output / name).open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fields, delimiter='\t', lineterminator='\n')
            writer.writeheader(); writer.writerows(rows)
    result = {'status':'passed_completed_case_snapshot' if pending else 'passed_full_genus_tree_audit',
              'audited_cases':len(summaries),'planned_cases':len(expected),'pending_cases':pending,
              'bootstrap_trees_read':1000*len(summaries),'internal_edges':len(branches),'warning_lines':len(warnings),
              'source_config_sha256':sha(config_path),'source_case_receipts':sources,'script_sha256':sha(Path(__file__)),
              'artifacts':{p.name:sha(p) for p in args.output.iterdir()},
              'interpretation':'All audited trees, bootstrap tip/edge/split grids, model reports and rounded UFB frequencies checked. SH-aLRT ranges checked, not independently recomputed. Numerical convergence, model adequacy, orthology, recombination, synonymous saturation and selection eligibility remain separate requirements.'}
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['pending_cases','source_case_receipts','artifacts']}, indent=2))


if __name__ == '__main__':
    main()
