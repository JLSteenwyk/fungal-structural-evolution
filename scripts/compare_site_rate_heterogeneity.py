#!/usr/bin/env python3
"""Compare audited Gamma4 and FreeRate4 estimates on identical observations and trees."""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from Bio import SeqIO
from scipy.stats import rankdata
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import read_table, sha
from prepare_paired_phylogenetic_inputs import write_table
from run_paired_marker_fits import tree_edges


def rank_correlation(x, y):
    if len(x) < 2 or len(set(x)) < 2 or len(set(y)) < 2:
        return ''
    return float(np.corrcoef(rankdata(x), rankdata(y))[0, 1])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['gamma', 'gamma-audit', 'free', 'free-audit', 'inputs', 'output']:
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use a new immutable output')
    checks = {}
    configs = {}
    rates = {}
    fit_summaries = {}
    for label, run, audit in [('G4', a.gamma, a.gamma_audit), ('R4', a.free, a.free_audit)]:
        checks[label] = checked_receipt(audit)
        receipt = json.loads((run / 'receipt.json').read_text())
        configs[label] = json.loads((run / 'config.json').read_text())
        if (checks[label]['status'] != 'passed_full_site_rate_output_audit'
                or checks[label]['rate_receipt_sha256'] != sha(run / 'receipt.json')
                or receipt['config_sha256'] != sha(run / 'config.json')
                or configs[label].get('heterogeneity', 'G4') != label):
            raise ValueError('Rate audit or model provenance differs')
        for key, h in receipt['fit_receipts'].items():
            marker, fit = key.split('/')
            rp = run / marker / (fit + '.receipt.json')
            if sha(rp) != h:
                raise ValueError('Changed fit receipt')
            for name, digest in json.loads(rp.read_text())['artifacts'].items():
                if sha(run / marker / name) != digest:
                    raise ValueError('Changed fit artifact')
        rate_rows = read_table(audit / 'site_rates.tsv')
        rates[label] = {(r['marker'], r['fit'], int(r['paired_column_1based'])):
                        float(r['posterior_mean_relative_rate']) for r in rate_rows}
        if len(rates[label]) != len(rate_rows):
            raise ValueError('Duplicate site identity')
        fit_rows = read_table(audit / 'fit_summary.tsv')
        fit_summaries[label] = {(r['marker'], r['fit']): r for r in fit_rows}
        if len(fit_summaries[label]) != len(fit_rows):
            raise ValueError('Duplicate fit identity')
    if (configs['G4']['source_receipts'] != configs['R4']['source_receipts']
            or configs['G4']['executable_sha256'] != configs['R4']['executable_sha256']
            or configs['G4']['source_receipts']['inputs'] != sha(a.inputs / 'receipt.json')
            or rates['G4'].keys() != rates['R4'].keys()
            or fit_summaries['G4'].keys() != fit_summaries['R4'].keys()):
        raise ValueError('Unmatched source, site or fit population')
    checked_receipt(a.inputs)
    sites, branches, summary = [], [], []
    for marker, fit in sorted(fit_summaries['G4']):
        n = int(fit_summaries['G4'][marker, fit]['sites'])
        g = [rates['G4'][marker, fit, i] for i in range(1, n + 1)]
        r = [rates['R4'][marker, fit, i] for i in range(1, n + 1)]
        sites.extend({'marker': marker, 'fit': fit, 'paired_column_1based': i,
                      'gamma_rate': x, 'freerate_rate': y, 'freerate_minus_gamma': y - x}
                     for i, (x, y) in enumerate(zip(g, r), 1))
        taxa = {x.id for x in SeqIO.parse(a.inputs / marker / 'aa.faa', 'fasta')}
        trees = {label: tree_edges(run / marker / (fit + '.treefile'), taxa)
                 for label, run in [('G4', a.gamma), ('R4', a.free)]}
        if trees['G4'].keys() != trees['R4'].keys():
            raise ValueError('Unmatched topology')
        edges = sorted(trees['G4'])
        bg = [trees['G4'][e] for e in edges]
        br = [trees['R4'][e] for e in edges]
        if any(not math.isfinite(v) or v < 0 for v in bg + br):
            raise ValueError('Invalid branch length')
        branches.extend({'marker': marker, 'fit': fit, 'split_taxa': ','.join(e),
                         'gamma_length': x, 'freerate_length': y,
                         'freerate_minus_gamma': y - x}
                        for e, x, y in zip(edges, bg, br))
        ll = {label: float(fit_summaries[label][marker, fit]['log_likelihood'])
              for label in ['G4', 'R4']}
        summary.append({'marker': marker, 'fit': fit, 'sites': n, 'branches': len(edges),
                        'site_rate_spearman': rank_correlation(g, r),
                        'median_absolute_site_rate_change': float(np.median(np.abs(np.array(r) - g))),
                        'branch_length_spearman': rank_correlation(bg, br),
                        'gamma_total_tree_length': math.fsum(bg),
                        'freerate_total_tree_length': math.fsum(br),
                        'gamma_log_likelihood': ll['G4'], 'freerate_log_likelihood': ll['R4'],
                        'freerate_minus_gamma_log_likelihood': ll['R4'] - ll['G4'],
                        'freerate_ll_lower_by_more_than_point_one': ll['R4'] < ll['G4'] - .1})
    if len(sites) != len(rates['G4']):
        raise ValueError('Incomplete output site grid')
    a.output.mkdir(parents=True)
    for name, rows in [('sites.tsv', sites), ('branches.tsv', branches), ('fit_summary.tsv', summary)]:
        write_table(a.output / name, rows)
    result = {'status': 'complete_matched_rate_heterogeneity_comparison',
              'fits': len(summary), 'site_comparisons': len(sites), 'branch_comparisons': len(branches),
              'source_receipts': {k: sha(getattr(a, k) / 'receipt.json')
                                  for k in ['gamma', 'gamma_audit', 'free', 'free_audit', 'inputs']},
              'script_sha256': sha(Path(__file__)),
              'interpretation': 'Conditional matched model sensitivity; average-tie rank correlations without p-values. No model selection, acceleration test or calibrated uncertainty. Site estimates are within-model relative multipliers; category labels are not homologous between models. FreeRate adds flexibility, so likelihood improvement alone does not establish adequacy. Lower likelihoods remain explicit optimization diagnostics.',
              'artifacts': {f.name: sha(f) for f in a.output.iterdir()}}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
