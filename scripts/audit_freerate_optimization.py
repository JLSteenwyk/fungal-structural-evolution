#!/usr/bin/env python3
"""Check all diagnostic fit outputs and summarize best available likelihoods."""
import argparse
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha
from audit_paired_site_rates import rounding_bound
from estimate_paired_site_rates import rate_rows
from prepare_paired_phylogenetic_inputs import write_table
from run_paired_marker_fits import tree_edges


def number(report, pattern):
    return float(re.search(pattern, report).group(1))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['diagnostics', 'free', 'output']:
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use a new immutable output')
    receipt = checked_receipt(a.diagnostics)
    cp = a.diagnostics / 'config.json'; config = json.loads(cp.read_text())
    if receipt['config_sha256'] != sha(cp) or config['source_receipts']['free'] != sha(a.free / 'receipt.json'):
        raise ValueError('Diagnostic provenance changed')
    for path, digest in config['pinned_files'].items():
        if sha(Path(path)) != digest:
            raise ValueError('Pinned source changed')
    requests = config['requests']
    if set(receipt['fit_receipts']) != {r['name'] for r in requests} or len(requests) != receipt['diagnostic_fits']:
        raise ValueError('Incomplete diagnostic grid')
    groups = defaultdict(list); rows = []; site_count = 0
    for request in requests:
        folder = a.diagnostics / request['name']; rp = folder / 'receipt.json'
        if sha(rp) != receipt['fit_receipts'][request['name']]:
            raise ValueError('Changed diagnostic receipt')
        fit = checked_receipt(folder)
        if fit['config_sha256'] != sha(cp):
            raise ValueError('Fit configuration changed')
        report = (folder / 'fit.iqtree').read_text()
        source_report = (a.free / request['marker'] / (request['fit'] + '.iqtree')).read_text()
        command = request['command']; model = command[command.index('-m') + 1].split('{')[0]
        if ('Model of substitution: ' + model + '\n' not in report
                or 'Model of rate heterogeneity: FreeRate with 4 categories' not in report
                or '-optfromgiven' not in command):
            raise ValueError('Model identity or optimization setting differs')
        pattern = r'Number of free parameters \(#branches \+ #model parameters\): (\d+)'
        if number(report, pattern) != number(source_report, pattern):
            raise ValueError('Unexpected parameter count: categories may have been fixed')
        alignment = Path(command[command.index('-s') + 1]); topology = Path(command[command.index('-te') + 1])
        if sha(alignment) != request['alignment_sha256'] or sha(topology) != request['topology_sha256']:
            raise ValueError('Inputs changed')
        taxa = {x.id for x in SeqIO.parse(alignment, 'fasta')}
        if set(tree_edges(folder / 'fit.treefile', taxa)) != set(tree_edges(topology, taxa)):
            raise ValueError('Topology changed')
        rates = rate_rows(folder / 'fit.rate', request['sites']); site_count += len(rates)
        block = report.split('Category  Relative_rate  Proportion', 1)[1].strip().split('\n\n', 1)[0]
        categories = {int(x[0]): (float(x[1]), float(x[2])) for line in block.splitlines() if len(x := line.split()) == 3}
        if set(categories) != {1, 2, 3, 4} or any(not math.isfinite(v) or v < 0 for pair in categories.values() for v in pair) or abs(sum(x[1] for x in categories.values()) - 1) > .00021:
            raise ValueError('Invalid fitted rate categories')
        for rate in rates:
            mean = float(rate['Rate']); modal = float(rate['Categorized_rate']); cat = int(rate['Category'])
            if abs(modal - categories[cat][0]) > .0001 or not min(x[0] for x in categories.values()) - .0001 <= mean <= max(x[0] for x in categories.values()) + .0001:
                raise ValueError('Rate/category bounds differ')
        lltext = re.search(r'Log-likelihood of the tree: ([\d.eE+-]+)', report).group(1)
        ll = float(lltext); lh = (folder / 'fit.sitelh').read_text().split()
        if lh[:3] != ['1', str(request['sites']), 'Site_Lh'] or len(lh) != request['sites'] + 3:
            raise ValueError('Likelihood grid differs')
        values = [float(x) for x in lh[3:]]
        if any(not math.isfinite(x) or x > 0 for x in values):
            raise ValueError('Invalid site likelihood')
        bound = sum(rounding_bound(x) for x in lh[3:]) + rounding_bound(lltext) + 1e-8
        if abs(math.fsum(values) - ll) > bound or ll != fit['log_likelihood']:
            raise ValueError('Likelihood accounting differs')
        row = {k: request[k] for k in ['name', 'marker', 'fit', 'start', 'optimizer']}
        row.update(log_likelihood=ll, minus_gamma=ll - request['gamma_log_likelihood'],
                   minus_original_freerate=ll - request['original_freerate_log_likelihood'])
        if any(row[k] != fit[k] for k in ['minus_gamma', 'minus_original_freerate']):
            raise ValueError('Likelihood comparison differs')
        rows.append(row); groups[request['marker'], request['fit']].append(row)
    best = []
    for key, fits in groups.items():
        if {(r['start'], r['optimizer']) for r in fits} != {(s, o) for s in ['gamma', 'freerate'] for o in ['EM', '2-BFGS']} or len(fits) != 4:
            raise ValueError('Missing initialization/optimizer combination')
        best.append(max(fits, key=lambda r: (r['log_likelihood'], r['name'])))
    a.output.mkdir(parents=True)
    write_table(a.output / 'all_fits.tsv', rows); write_table(a.output / 'best_diagnostics.tsv', best)
    result = {'status': 'passed_full_freerate_optimization_diagnostics', 'diagnostic_fits': len(rows),
              'source_fits': len(best), 'site_rate_rows': site_count,
              'source_fits_with_diagnostic_at_least_gamma': sum(r['minus_gamma'] >= 0 for r in best),
              'diagnostic_receipt_sha256': sha(a.diagnostics / 'receipt.json'),
              'script_sha256': sha(Path(__file__)),
              'interpretation': 'All diagnostic grids, hashes, model/parameter counts, fixed topologies, site-rate/category bounds and likelihood sums checked. Best diagnostics maximize observed likelihood among four explicit refits; global optimality is not established. Original outputs remain unchanged. This diagnoses initialization/optimizer sensitivity, not biological model adequacy or selection.',
              'artifacts': {f.name: sha(f) for f in a.output.iterdir()}}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
