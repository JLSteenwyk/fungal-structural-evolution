#!/usr/bin/env python3
"""Check every comparison value against raw reports, rates and independently collected tree edges."""
import argparse
import csv
import json
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path
from Bio import Phylo
from scipy.stats import spearmanr
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def likelihood(path):
    value = float(re.search(r'Log-likelihood of the tree: ([\d.eE+-]+)', path.read_text()).group(1))
    if not math.isfinite(value):
        raise ValueError('Nonfinite likelihood')
    return value


def raw_rates(path):
    lines = [line.split() for line in path.read_text().splitlines()
             if line.strip() and not line.startswith('#')]
    site, rate = lines[0].index('Site'), lines[0].index('Rate')
    if [int(row[site]) for row in lines[1:]] != list(range(1, len(lines))):
        raise ValueError('Invalid raw site grid')
    values = [float(row[rate]) for row in lines[1:]]
    if any(not math.isfinite(v) or v < 0 for v in values):
        raise ValueError('Invalid raw rate')
    return values


def edges(path):
    tree = Phylo.read(path, 'newick')
    tips = [tip.name for tip in tree.get_terminals()]
    if len(set(tips)) != len(tips):
        raise ValueError('Repeated tree tips')
    universe = frozenset(tips)
    descendants, result = {}, defaultdict(float)
    for node in tree.find_clades(order='postorder'):
        descendants[node] = (frozenset([node.name]) if node.is_terminal() else
                             frozenset().union(*(descendants[c] for c in node.clades)))
        if node is tree.root:
            continue
        side, other = sorted(descendants[node]), sorted(universe - descendants[node])
        selected = side if (len(side), side) < (len(other), other) else other
        if node.branch_length is None or not math.isfinite(node.branch_length) or node.branch_length < 0:
            raise ValueError('Invalid raw branch')
        result[','.join(selected)] += node.branch_length
    return dict(result)


def close(value, expected):
    if not math.isclose(float(value), expected, rel_tol=1e-10, abs_tol=1e-10):
        raise ValueError(f'Mismatch {value} != {expected}')


def rank_check(value, x, y):
    if len(x) < 2 or len(set(x)) < 2 or len(set(y)) < 2:
        if value != '':
            raise ValueError('Unavailable correlation reported')
    else:
        close(value, float(spearmanr(x, y).statistic))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('comparison', 'gamma', 'free', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--optimization', type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    receipt = checked_receipt(args.comparison)
    if receipt['status'] != 'complete_matched_rate_heterogeneity_comparison':
        raise ValueError('Incomplete comparison')
    if ('optimization' in receipt['source_receipts']) != (args.optimization is not None):
        raise ValueError('Provide optimization source exactly when comparison used it')
    source_receipts = {}
    for source in ('gamma', 'free') + (('optimization',) if args.optimization else ()):
        folder = getattr(args, source)
        if sha(folder / 'receipt.json') != receipt['source_receipts'][source]:
            raise ValueError('Source receipt changed')
        source_receipts[source] = json.loads((folder / 'receipt.json').read_text())
        for identity, digest in source_receipts[source]['fit_receipts'].items():
            path = (folder / identity / 'receipt.json' if source == 'optimization'
                    else folder / (identity + '.receipt.json'))
            if sha(path) != digest:
                raise ValueError('Fit receipt changed')
            fit = json.loads(path.read_text())
            for name, expected in fit['artifacts'].items():
                if sha(path.parent / name) != expected:
                    raise ValueError('Fit artifact changed')
    summaries = table(args.comparison / 'fit_summary.tsv')
    keys = [(r['marker'], r['fit']) for r in summaries]
    expected = {tuple(k.split('/')) for k in source_receipts['gamma']['fit_receipts']}
    if len(keys) != len(set(keys)) or set(keys) != expected or set(source_receipts['gamma']['fit_receipts']) != set(source_receipts['free']['fit_receipts']):
        raise ValueError('Comparison fit grid differs')
    grouped = {}
    for kind in ('sites', 'branches'):
        grouped[kind] = defaultdict(list)
        for row in table(args.comparison / (kind + '.tsv')):
            grouped[kind][row['marker'], row['fit']].append(row)
        if set(grouped[kind]) != expected:
            raise ValueError('Comparison row grid differs')
    requests = defaultdict(list)
    if args.optimization:
        config = json.loads((args.optimization / 'config.json').read_text())
        if sha(args.optimization / 'config.json') != source_receipts['optimization']['config_sha256']:
            raise ValueError('Optimization configuration changed')
        for request in config['requests']:
            requests[request['marker'], request['fit']].append(request)
        if set(requests) != expected:
            raise ValueError('Optimization grid differs')
    site_count = branch_count = optimized = 0
    for row in summaries:
        marker, fit = identity = (row['marker'], row['fit'])
        gp, fp = args.gamma / marker / fit, args.free / marker / fit
        original_ll = likelihood(fp.with_suffix('.iqtree'))
        if args.optimization:
            group = requests[identity]
            if len(group) != 4 or {(q['start'], q['optimizer']) for q in group} != {(s, o) for s in ('gamma', 'freerate') for o in ('EM', '2-BFGS')}:
                raise ValueError('Missing refit combination')
            candidates = [(likelihood(args.optimization / q['name'] / 'fit.iqtree'), q['name']) for q in group]
            best_ll, best_name = max(candidates)
            use_optimized = best_ll > original_ll
            if row['selected_freerate_source'] != ('optimized' if use_optimized else 'original') or row['selected_diagnostic_name'] != (best_name if use_optimized else ''):
                raise ValueError('Incorrect selected fit')
            close(row['original_freerate_log_likelihood'], original_ll)
            close(row['best_diagnostic_log_likelihood'], best_ll)
            if use_optimized:
                fp = args.optimization / best_name / 'fit'
                optimized += 1
        g, f = raw_rates(gp.with_suffix('.rate')), raw_rates(fp.with_suffix('.rate'))
        sr = sorted(grouped['sites'][identity], key=lambda r: int(r['paired_column_1based']))
        if len(g) != len(f) or len(g) != len(sr) or len(g) != int(row['sites']) or [int(r['paired_column_1based']) for r in sr] != list(range(1, len(g) + 1)):
            raise ValueError('Site grid mismatch')
        for record, x, y in zip(sr, g, f):
            close(record['gamma_rate'], x); close(record['freerate_rate'], y)
            close(record['freerate_minus_gamma'], y - x)
        rank_check(row['site_rate_spearman'], g, f)
        close(row['median_absolute_site_rate_change'], statistics.median(abs(x - y) for x, y in zip(g, f)))
        ge, fe = edges(gp.with_suffix('.treefile')), edges(fp.with_suffix('.treefile'))
        br = grouped['branches'][identity]
        if ge.keys() != fe.keys() or {r['split_taxa'] for r in br} != ge.keys() or len(br) != len(ge) or len(br) != int(row['branches']):
            raise ValueError('Branch grid mismatch')
        for record in br:
            x, y = ge[record['split_taxa']], fe[record['split_taxa']]
            close(record['gamma_length'], x); close(record['freerate_length'], y)
            close(record['freerate_minus_gamma'], y - x)
        rank_check(row['branch_length_spearman'], [ge[k] for k in sorted(ge)], [fe[k] for k in sorted(ge)])
        close(row['gamma_total_tree_length'], math.fsum(ge.values()))
        close(row['freerate_total_tree_length'], math.fsum(fe.values()))
        gl, fl = likelihood(gp.with_suffix('.iqtree')), likelihood(fp.with_suffix('.iqtree'))
        close(row['gamma_log_likelihood'], gl); close(row['freerate_log_likelihood'], fl)
        close(row['freerate_minus_gamma_log_likelihood'], fl - gl)
        if row['freerate_ll_lower_by_more_than_point_one'] != str(fl < gl - .1):
            raise ValueError('Likelihood flag mismatch')
        site_count += len(sr); branch_count += len(br)
    if len(summaries) != receipt['fits'] or site_count != receipt['site_comparisons'] or branch_count != receipt['branch_comparisons'] or (args.optimization and optimized != receipt['optimized_fits_selected']):
        raise ValueError('Receipt accounting mismatch')
    result = {'status': 'passed_full_rate_comparison_readback', 'fits': len(summaries),
              'sites': site_count, 'branches': branch_count, 'optimized_fits_selected': optimized,
              'comparison_receipt_sha256': sha(args.comparison / 'receipt.json'),
              'script_sha256': sha(Path(__file__)),
              'scope': 'All raw rate values, tree edges, likelihoods, grid identities and summary statistics recomputed. When optimization is present, all four raw diagnostic likelihoods determine selection independently. Tree parsing uses Bio.Phylo with separate postorder edge collection; Spearman uses scipy.stats.spearmanr. No independent likelihood optimization, model adequacy or biological inference.'}
    args.output.mkdir(parents=True)
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
