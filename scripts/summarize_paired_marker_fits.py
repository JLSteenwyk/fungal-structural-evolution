#!/usr/bin/env python3
"""Audit completed paired branch fits and retain numerical/model warnings."""
import argparse
import csv
import json
import math
import re
from pathlib import Path
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import sha
from run_paired_marker_fits import tree_edges
from prepare_paired_phylogenetic_inputs import write_table


def number(text, pattern):
    match = re.search(pattern, text)
    if not match:
        raise ValueError('Missing report field: ' + pattern)
    value = float(match.group(1))
    if not math.isfinite(value):
        raise ValueError('Nonfinite report value')
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['inputs', 'models', 'fits', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    inputs = checked_receipt(args.inputs)
    checked_receipt(args.models)
    receipt = json.loads((args.fits / 'receipt.json').read_text())
    config_path = args.fits / 'config.json'
    config = json.loads(config_path.read_text())
    if (receipt['status'] != 'complete_matched_topology_point_estimates'
            or receipt['config_sha256'] != sha(config_path)
            or config['input_receipt_sha256'] != sha(args.inputs / 'receipt.json')
            or config['model_receipt_sha256'] != sha(args.models / 'receipt.json')):
        raise ValueError('Run provenance mismatch')
    markers = {r['marker']: r for r in csv.DictReader((args.inputs / 'marker_summary.tsv').open(), delimiter='\t')
               if r['status'] == 'ready_for_inference'}
    if len(receipt['results']) != inputs['ready_markers'] or {r['marker'] for r in receipt['results']} != set(markers):
        raise ValueError('Incomplete marker fits')
    if args.output.exists():
        raise FileExistsError('Use a new immutable audit output')
    fits, branches, warnings = [], [], []
    for result in receipt['results']:
        marker = result['marker']
        folder = args.fits / marker
        branch_path = folder / 'paired_branches.tsv'
        if sha(branch_path) != result['paired_branches_sha256']:
            raise ValueError('Changed paired branch table')
        aa = {r.id: str(r.seq) for r in SeqIO.parse(args.inputs / marker / 'aa.faa', 'fasta')}
        state = {r.id: str(r.seq) for r in SeqIO.parse(args.inputs / marker / '3di.faa', 'fasta')}
        if set(aa) != set(state):
            raise ValueError('Paired taxa differ')
        for taxon in aa:
            if [c == '?' for c in aa[taxon]] != [c == '?' for c in state[taxon]]:
                raise ValueError('Observation masks differ')
        all_edges = {}
        for label in ['aa', '3di_af', '3di_af_empirical', '3di_llm']:
            fit = json.loads((folder / (label + '.receipt.json')).read_text())
            run_path = folder / (label + '.config.json')
            run = json.loads(run_path.read_text())
            alignment = args.inputs / marker / ('aa.faa' if label == 'aa' else '3di.faa')
            if (fit['config_sha256'] != sha(run_path) or run['parent_config_sha256'] != sha(config_path)
                    or run['alignment_sha256'] != sha(alignment)):
                raise ValueError('Fit configuration or alignment changed')
            for name, checksum in fit['artifacts'].items():
                if sha(folder / name) != checksum:
                    raise ValueError('Changed fit artifact')
            if label != 'aa' and run['topology_sha256'] != sha(folder / 'aa.treefile'):
                raise ValueError('Fixed topology source changed')
            report = (folder / (label + '.iqtree')).read_text()
            expected_model = run['command'][run['command'].index('-m') + 1]
            if 'Model of substitution: ' + expected_model + '\n' not in report:
                raise ValueError('Report substitution model differs from command')
            edges = tree_edges(folder / (label + '.treefile'), set(aa))
            all_edges[label] = edges
            total = sum(edges.values())
            reported_total = number(report, r'Total tree length \(sum of branch lengths\): ([\d.eE+-]+)')
            if not math.isclose(total, reported_total, abs_tol=.000051):
                raise ValueError('Reported and parsed total tree lengths differ')
            log = (folder / (label + '.log')).read_text()
            warning_lines = sorted(set(line.strip() for line in (report + '\n' + log).splitlines()
                                       if re.match(r'\s*(WARNING|ERROR):', line)))
            warnings.extend({'marker': marker, 'fit': label, 'warning': line} for line in warning_lines)
            fits.append({'marker': marker, 'fit': label, 'taxa': len(aa), 'columns': len(next(iter(aa.values()))),
                'log_likelihood': number(report, r'Log-likelihood of the tree: ([\d.eE+-]+)'),
                'gamma_alpha': number(report, r'Gamma shape alpha: ([\d.eE+-]+)'),
                'total_branch_length': total, 'branches': len(edges),
                'branches_le_1e_minus5': sum(v <= 1e-5 for v in edges.values()),
                'branches_ge_10': sum(v >= 10 for v in edges.values()), 'warning_count': len(warning_lines),
                'elapsed_seconds': fit['elapsed_seconds']})
        table = list(csv.DictReader(branch_path.open(), delimiter='\t'))
        if len(table) != len(all_edges['aa']):
            raise ValueError('Branch table length mismatch')
        for label, edges in all_edges.items():
            if set(edges) != set(all_edges['aa']):
                raise ValueError('Paired topology mismatch')
            for row in table:
                split = tuple(row['split_taxa'].split(','))
                if split not in edges or not math.isclose(float(row[label + '_branch_length']), edges[split], abs_tol=1e-12):
                    raise ValueError('Paired branch values differ from trees')
        branches.extend(table)
    args.output.mkdir(parents=True)
    write_table(args.output / 'fit_summary.tsv', fits)
    write_table(args.output / 'paired_branches.tsv', branches)
    if warnings:
        write_table(args.output / 'warnings.tsv', warnings)
    audit = {'status': 'complete_paired_fit_audit', 'markers': len(markers), 'fits': len(fits),
        'paired_branches': len(branches), 'fits_with_warnings': sum(r['warning_count'] > 0 for r in fits),
        'fits_with_near_zero_branches': sum(r['branches_le_1e_minus5'] > 0 for r in fits),
        'fits_with_branches_ge_10': sum(r['branches_ge_10'] > 0 for r in fits),
        'fit_receipt_sha256': sha(args.fits / 'receipt.json'), 'script_sha256': sha(Path(__file__)),
        'interpretation': 'Verified model identity, provenance, paired masks/topologies and branch tables. Point estimates are not confidence intervals or acceleration tests. Likelihoods must not be compared across AA and 3Di data; model comparisons within the same 3Di data require parameter penalties. Near-zero lengths warrant uncertainty analysis, not division into branch-length ratios.',
        'artifacts': {p.name: sha(p) for p in args.output.iterdir()}}
    (args.output / 'receipt.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
