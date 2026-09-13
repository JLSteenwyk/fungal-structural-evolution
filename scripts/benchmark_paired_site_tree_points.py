#!/usr/bin/env python3
"""Join full fitted tree path point estimates to exact paired-site geometry."""
import argparse
import hashlib
import json
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import numpy as np
from Bio import Phylo, SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha, read_table
from prepare_paired_phylogenetic_inputs import write_table
from run_paired_marker_fits import tree_edges


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['inputs', 'fits', 'audit', 'geometry', 'review', 'output']:
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use immutable output')
    inputs = checked_receipt(a.inputs); audit = checked_receipt(a.audit); geometry = checked_receipt(a.geometry)
    fits = json.loads((a.fits / 'receipt.json').read_text()); config = json.loads((a.fits / 'config.json').read_text())
    if audit['fit_receipt_sha256'] != sha(a.fits / 'receipt.json') or fits['config_sha256'] != sha(a.fits / 'config.json') or config['input_receipt_sha256'] != sha(a.inputs / 'receipt.json') or geometry['source_receipts']['inputs']['sha256'] != sha(a.inputs / 'receipt.json'):
        raise ValueError('Source lineage mismatch')
    caveats = {r['marker']: r['status'] for r in json.loads(a.review.read_text())['records']}
    by_marker = defaultdict(list)
    for accepted, filename in [(True, 'paired_site_geometry.tsv'), (False, 'coverage_exclusions.tsv')]:
        for row in read_table(a.geometry / filename):
            by_marker[row['marker']].append((accepted, row))
    fitted = {r['marker']: r for r in fits['results']}
    if len(fitted) != len(fits['results']) or set(fitted) != set(by_marker) or len(fitted) != inputs['ready_markers']:
        raise ValueError('Incomplete marker universe')
    results, excluded, summaries = [], [], []
    independent_checks = 0
    labels = ['aa', '3di_af', '3di_af_empirical', '3di_llm']
    for marker in sorted(fitted):
        taxa = sorted(r.id for r in SeqIO.parse(a.inputs / marker / 'aa.faa', 'fasta'))
        index = {t: i for i, t in enumerate(taxa)}
        if len(index) != len(taxa):
            raise ValueError('Duplicate input taxon')
        expected = set(combinations(taxa, 2)); seen = set(); distances = {}; trees = {}; first_edges = None
        if sha(a.fits / marker / 'paired_branches.tsv') != fitted[marker]['paired_branches_sha256']:
            raise ValueError('Branch table changed')
        for label in labels:
            path = a.fits / marker / (label + '.treefile')
            fr = json.loads((a.fits / marker / (label + '.receipt.json')).read_text())
            if sha(path) != fr['artifacts'][path.name]:
                raise ValueError('Tree changed')
            edges = tree_edges(path, set(taxa)); trees[label] = Phylo.read(path, 'newick')
            if first_edges is not None and set(edges) != first_edges:
                raise ValueError('Topology differs')
            first_edges = set(edges)
            matrix = np.zeros((len(taxa), len(taxa)))
            for side, length in edges.items():
                inside = np.array([t in side for t in taxa])
                matrix += np.logical_xor(inside[:, None], inside[None, :]) * length
            distances[label] = matrix
        check_pairs = set(sorted(expected, key=lambda pair: hashlib.sha256((marker+'|'+ '|'.join(pair)).encode()).digest())[:10])
        accepted_count = 0
        for accepted, row in by_marker[marker]:
            pair = row['taxon_a'], row['taxon_b']
            if pair not in expected or pair in seen:
                raise ValueError('Unexpected or duplicate geometry pair')
            seen.add(pair); ia, ib = map(index.__getitem__, pair)
            additions = {label + '_tree_path_point': float(distances[label][ia, ib]) for label in labels}
            if pair in check_pairs:
                for label in labels:
                    if not np.isclose(additions[label+'_tree_path_point'], trees[label].distance(*pair), atol=1e-10, rtol=1e-10):
                        raise ValueError('Independent tree traversal differs')
                    independent_checks += 1
            output = dict(row, **additions, marker_review_status=caveats.get(marker, 'no_current_copy_review_flag'), uncertainty_status='point_only_resampling_pending')
            (results if accepted else excluded).append(output)
            accepted_count += accepted
        if seen != expected:
            raise ValueError('Incomplete per-marker pair grid')
        summaries.append({'marker': marker, 'taxa': len(taxa), 'pairs': len(expected), 'accepted_geometry_pairs': accepted_count,
                          'excluded_geometry_pairs': len(expected)-accepted_count,
                          'marker_review_status': caveats.get(marker, 'no_current_copy_review_flag')})
    if len(results) != geometry['accepted_taxon_pairs'] or len(excluded) != geometry['excluded_taxon_pairs']:
        raise ValueError('Geometry totals differ')
    a.output.mkdir(parents=True)
    for name, rows in [('path_geometry_points.tsv', results), ('geometry_exclusions_with_paths.tsv', excluded), ('marker_summary.tsv', summaries)]:
        write_table(a.output / name, rows)
    result = {'status': 'complete_paired_site_tree_path_point_benchmark',
              'source_receipts': {name: sha(getattr(a, name) / 'receipt.json') for name in ['inputs', 'fits', 'audit', 'geometry']},
              'marker_review_sha256': sha(a.review), 'script_sha256': sha(Path(__file__)),
              'markers': len(fitted), 'accepted_pairs': len(results), 'excluded_pairs': len(excluded),
              'independent_tree_path_checks': independent_checks,
              'interpretation': 'All source geometry rows retained unchanged with fixed-topology path point estimates. Geometry uses shared qualified paired sites; tree fits use the full paired alignment with taxon-specific missingness. No fitted relationship, additive-geometry assumption, confidence interval, phylogenetically adjusted coupling or acceleration claim. Copy caveats remain explicit; no flag does not establish orthology.',
              'artifacts': {f.name: sha(f) for f in a.output.iterdir()}}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n'); print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
