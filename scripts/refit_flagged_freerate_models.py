#!/usr/bin/env python3
"""Diagnose lower-likelihood FreeRate fits using explicit starts and two optimizers."""
import argparse
import fcntl
import json
import math
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import read_table, sha
from estimate_paired_site_rates import rate_rows
from prepare_paired_phylogenetic_inputs import write_table
from run_paired_marker_fits import tree_edges


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['comparison', 'gamma', 'free', 'output']:
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--all-fits', action='store_true', help='Apply the same four refits to every comparison fit')
    a = p.parse_args()
    comparison = checked_receipt(a.comparison)
    for name in ['gamma', 'free']:
        if comparison['source_receipts'][name] != sha(getattr(a, name) / 'receipt.json'):
            raise ValueError('Source runs differ')
    all_rows = read_table(a.comparison / 'fit_summary.tsv')
    rows = [r for r in all_rows if a.all_fits or r['freerate_ll_lower_by_more_than_point_one'] == 'True']
    if len({(r['marker'], r['fit']) for r in rows}) != len(rows) or len(all_rows) != comparison['fits']:
        raise ValueError('Duplicate or incomplete comparison fit grid')
    if not rows:
        raise ValueError('No flagged fits')
    pins = {}; requests = []
    for row in rows:
        marker, fit = row['marker'], row['fit']
        base = json.loads((a.free / marker / (fit + '.config.json')).read_text())['command']
        model = base[base.index('-m') + 1]
        if not model.endswith('+R4'):
            raise ValueError('Expected FreeRate4 source model')
        alignment = Path(base[base.index('-s') + 1])
        for path in [alignment, Path(base[0])]:
            pins[str(path)] = sha(path)
        matrix = model.split('+')[0]
        if matrix != 'LG':
            pins[matrix] = sha(Path(matrix))
        for start, run in [('gamma', a.gamma), ('freerate', a.free)]:
            rp = run / marker / (fit + '.receipt.json')
            receipt = json.loads(rp.read_text())
            run_receipt = json.loads((run / 'receipt.json').read_text())
            if sha(rp) != run_receipt['fit_receipts'][marker + '/' + fit]:
                raise ValueError('Source fit receipt changed')
            for name, digest in receipt['artifacts'].items():
                path = run / marker / name
                if sha(path) != digest:
                    raise ValueError('Source artifact changed')
                pins[str(path)] = digest
            report = (run / marker / (fit + '.iqtree')).read_text()
            block = report.split('Category  Relative_rate  Proportion', 1)[1].strip().split('\n\n', 1)[0].split('Relative rates', 1)[0]
            categories = [line.split() for line in block.splitlines() if line.strip()]
            if [int(x[0]) for x in categories] != [1, 2, 3, 4]:
                raise ValueError('Invalid starting categories')
            weights = [float(x[2]) for x in categories]
            rates = [max(1e-6, float(x[1])) for x in categories]
            if any(not math.isfinite(x) or x <= 0 for x in weights + rates):
                raise ValueError('Invalid starting parameters')
            weights = [x / sum(weights) for x in weights]
            mean = sum(w * r for w, r in zip(weights, rates))
            rates = [x / mean for x in rates]
            params = ','.join(format(x, '.17g') for pair in zip(weights, rates) for x in pair)
            topology = (run / marker / (fit + '.treefile')).resolve()
            for optimizer in ['EM', '2-BFGS']:
                name = marker + '__' + fit + '__' + start + '__' + optimizer
                prefix = (a.output / name / 'fit').resolve()
                command = list(base)
                command[command.index('-m') + 1] = model + '{' + params + '}'
                command[command.index('-te') + 1] = str(topology)
                command[command.index('--prefix') + 1] = str(prefix)
                command.extend(['-optfromgiven', '-optalg', optimizer, '--epsilon', '0.000001'])
                requests.append({'name': name, 'marker': marker, 'fit': fit, 'start': start,
                                 'optimizer': optimizer, 'command': command,
                                 'starting_weights': weights, 'starting_rates': rates,
                                 'alignment_sha256': sha(alignment), 'topology_sha256': sha(topology),
                                 'sites': int(row['sites']), 'gamma_log_likelihood': float(row['gamma_log_likelihood']),
                                 'original_freerate_log_likelihood': float(row['freerate_log_likelihood'])})
    a.output.mkdir(parents=True, exist_ok=True)
    lock = (a.output / '.lock').open('w'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    config = {'source_receipts': {name: sha(getattr(a, name) / 'receipt.json') for name in ['comparison', 'gamma', 'free']},
              'pinned_files': pins, 'script_sha256': sha(Path(__file__)), 'requests': requests,
              'workers': 4, 'threads_per_fit': 1, 'memory_per_fit_gb': 2,
              'interpretation': 'Diagnostic refits of every FreeRate fit more than 0.1 log units below Gamma. Two starting points from reported category weights/rates and corresponding fitted trees; floor rounded zero rates at 1e-6 and normalize. Starts approximate printed estimates, not exact parameter replay. Both EM and 2-BFGS optimize from supplied values at epsilon 1e-6; no topology search. No automatic replacement or calibrated inference.'}
    cp = a.output / 'config.json'
    config['scope'] = 'all_comparison_fits' if a.all_fits else 'lower_likelihood_flagged_fits'
    if a.all_fits:
        config['interpretation'] = config['interpretation'].replace('Diagnostic refits of every FreeRate fit more than 0.1 log units below Gamma.', 'Consistent four-refit sensitivity for every FreeRate fit in the completed comparison, including originally unflagged fits.')
    if cp.exists() and json.loads(cp.read_text()) != config:
        raise ValueError('Configuration changed; use a new output')
    cp.write_text(json.dumps(config, indent=2) + '\n')
    config_hash = sha(cp)
    def run(request):
        folder = a.output / request['name']; folder.mkdir(exist_ok=True)
        rp = folder / 'receipt.json'
        if rp.exists():
            result = json.loads(rp.read_text())
            if result['config_sha256'] != config_hash or any(sha(folder / p) != h for p, h in result['artifacts'].items()):
                raise ValueError('Completed diagnostic changed')
            return result
        command = request['command']
        with (folder / 'stdout.log').open('a') as f:
            subprocess.run(command, stdout=f, stderr=subprocess.STDOUT, check=True)
        taxa = {x.id for x in SeqIO.parse(command[command.index('-s') + 1], 'fasta')}
        if set(tree_edges(folder / 'fit.treefile', taxa)) != set(tree_edges(Path(command[command.index('-te') + 1]), taxa)):
            raise ValueError('Diagnostic topology changed')
        rate_rows(folder / 'fit.rate', request['sites'])
        report = (folder / 'fit.iqtree').read_text()
        value = float(re.search(r'Log-likelihood of the tree: ([\d.eE+-]+)', report).group(1))
        if not math.isfinite(value):
            raise ValueError('Nonfinite likelihood')
        result = {'name': request['name'], 'config_sha256': config_hash,
                  'log_likelihood': value, 'minus_gamma': value - request['gamma_log_likelihood'],
                  'minus_original_freerate': value - request['original_freerate_log_likelihood'],
                  'artifacts': {p.name: sha(p) for p in folder.iterdir() if p.is_file()}}
        rp.write_text(json.dumps(result, indent=2) + '\n')
        print(request['name'], result['minus_gamma'], flush=True)
        return result
    results = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        for future in as_completed([pool.submit(run, r) for r in requests]):
            results.append(future.result())
    if any(sha(Path(p)) != h for p, h in pins.items()):
        raise ValueError('Source changed during diagnostic')
    summary = [{k: r[k] for k in ['name', 'log_likelihood', 'minus_gamma', 'minus_original_freerate']}
               for r in sorted(results, key=lambda r: r['name'])]
    write_table(a.output / 'fit_summary.tsv', summary)
    result = {'status': 'complete_flagged_freerate_diagnostic_execution', 'config_sha256': config_hash,
              'source_fits': len(rows), 'flagged_source_fits': sum(r['freerate_ll_lower_by_more_than_point_one'] == 'True' for r in rows), 'diagnostic_fits': len(results),
              'artifacts': {'fit_summary.tsv': sha(a.output / 'fit_summary.tsv')},
              'fit_receipts': {r['name']: sha(a.output / r['name'] / 'receipt.json') for r in results},
              'interpretation': config['interpretation']}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
