#!/usr/bin/env python3
"""Run isolated stricter Gamma fits from pinned alternative branch-length starts."""
import argparse
import json
import math
import os
from pathlib import Path
import re
import subprocess
from Bio import SeqIO
from audit_busco_gene_copies import sha
from run_paired_marker_fits import tree_edges
from audit_paired_site_rates import rounding_bound
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    a = p.parse_args(); c = json.loads(a.plan.read_text())
    for name, digest in c['pins'].items():
        if sha(Path(name)) != digest:
            raise ValueError('Changed source: ' + name)
    output = Path(c['output'])
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True); os.sched_setaffinity(0, c['cpu_affinity'])
    alignment = list(SeqIO.parse(c['alignment'], 'fasta'))
    taxa = {r.id for r in alignment}; width = len(alignment[0].seq)
    if len(taxa) != len(alignment) or any(len(r.seq) != width for r in alignment):
        raise ValueError('Alignment grid differs')
    topology = set(tree_edges(Path(c['reference_tree']), taxa)); results = []
    for request in c['requests']:
        prefix = output / request['name']
        with prefix.with_suffix('.stdout').open('w') as log:
            subprocess.run(request['command'], check=True, stdout=log, stderr=subprocess.STDOUT,
                env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'))
        report = prefix.with_suffix('.iqtree').read_text()
        if 'Model of substitution: LG+F+G4\n' not in report or set(tree_edges(prefix.with_suffix('.treefile'), taxa)) != topology:
            raise ValueError('Model/topology changed')
        likelihood = float(re.search(r'Log-likelihood of the tree: ([\d.eE+-]+)', report).group(1))
        alpha = float(re.search(r'Gamma shape alpha: ([\d.eE+-]+)', report).group(1))
        rate = [line.split() for line in prefix.with_suffix('.rate').read_text().splitlines() if line.strip() and not line.startswith('#')]
        if rate[0] != ['Site','Rate','Cat','C_Rate'] or [int(x[0]) for x in rate[1:]] != list(range(1,width+1)):
            raise ValueError('Rate grid differs')
        if any(not math.isfinite(float(x[1])) or float(x[1]) < 0 for x in rate[1:]):
            raise ValueError('Invalid rate')
        raw = prefix.with_suffix('.sitelh').read_text().split()
        if raw[:2] != ['1',str(width)] or len(raw) != width+3:
            raise ValueError('Site likelihood grid differs')
        error = abs(math.fsum(float(x) for x in raw[3:])-likelihood)
        if error > sum(rounding_bound(x) for x in raw[3:])+.00005+1e-8:
            raise ValueError('Likelihood sum differs')
        results.append(dict(start=request['name'], log_likelihood=likelihood, gamma_alpha=alpha,
            minus_original=likelihood-c['original_likelihood'], minus_export=likelihood-c['export_likelihood'],
            site_likelihood_sum_error=error, sites=width))
    write_table(output / 'summary.tsv', results)
    receipt = dict(status='complete_isolated_gamma_export_refits', plan_sha256=sha(a.plan),
        script_sha256=sha(Path(__file__)), fits=len(results), sites_per_fit=width,
        interpretation='Stricter fixed-topology optimization from two branch-length starts; source estimates preserved. No global-optimum claim or automatic replacement. Independent rate and branch readback and downstream propagation remain required.',
        artifacts={f.name:sha(f) for f in output.iterdir()})
    (output / 'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(results,indent=2))


if __name__ == '__main__':
    main()
