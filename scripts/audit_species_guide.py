#!/usr/bin/env python3
"""Freeze a completed full-taxon guide with report, tree and composition readback."""
import argparse
import csv
import io
import json
import math
from pathlib import Path
import re
import shutil
from Bio import Phylo, SeqIO
from audit_busco_gene_copies import ROOT, sha


def edges(tree, expected):
    tips = [x.name for x in tree.get_terminals()]
    if len(tips) != len(set(tips)) or set(tips) != expected:
        raise ValueError('Guide tip grid differs')
    result = {}
    for node in tree.find_clades():
        if node is tree.root:
            if node.branch_length not in [None, 0]:
                raise ValueError('Unexpected root stem')
            continue
        value = node.branch_length
        if value is None or not math.isfinite(value) or value < 0:
            raise ValueError('Invalid guide branch length')
        side = frozenset(t.name for t in node.get_terminals())
        left, right = tuple(sorted(side)), tuple(sorted(expected - side))
        key = min(left, right, key=lambda s: (len(s), s))
        if key in result:
            raise ValueError('Duplicate unrooted edge')
        result[key] = value
    if len(result) != 2 * len(expected) - 3:
        raise ValueError('Unexpected binary guide edge count')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['guide', 'matrix', 'output']:
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use a new immutable guide audit')
    source = json.loads((a.guide / 'receipt.json').read_text())
    if source['status'] != 'guide_inferred_not_final' or source['returncode'] != 0:
        raise ValueError('Completed guide execution required')
    matrix = json.loads((a.matrix / 'receipt.json').read_text())
    for name, expected in matrix['artifacts'].items():
        if sha(a.matrix / name) != expected:
            raise ValueError('Changed matrix artifact')
    if source['input_sha256'] != sha(a.matrix / 'matrix.faa'):
        raise ValueError('Guide input differs')
    manifest = ROOT / 'metadata/analysis_manifest.tsv'
    if matrix['manifest_sha256'] != sha(manifest):
        raise ValueError('Changed sampling manifest')
    with manifest.open() as handle:
        taxa = list(csv.DictReader(handle, delimiter='\t'))
    expected = {r['taxon_id'] for r in taxa}
    if len(expected) != 526 or len(taxa) != 526:
        raise ValueError('Full sampling grid differs')
    with (a.matrix / 'matrix.faa').open() as handle:
        seqs = list(SeqIO.parse(handle, 'fasta'))
    if len(seqs) != 526 or {r.id for r in seqs} != expected or any(len(r.seq) != matrix['columns'] for r in seqs):
        raise ValueError('Matrix identity or dimensions differ')
    tree_path = a.guide / 'guide.treefile'
    if sha(tree_path) != source['tree_sha256']:
        raise ValueError('Changed completed tree')
    tree = Phylo.read(tree_path, 'newick'); observed = edges(tree, expected)
    report = (a.guide / 'guide.iqtree').read_text()
    reported_trees = [line.strip() for line in report.splitlines() if line.startswith('(') and line.rstrip().endswith(';')]
    if len(reported_trees) != 1 or edges(Phylo.read(io.StringIO(reported_trees[0]), 'newick'), expected) != observed:
        raise ValueError('Report/tree edge mismatch')
    total = float(re.search(r'Total tree length \(sum of branch lengths\): ([\d.]+)', report)[1])
    if abs(math.fsum(observed.values()) - total) > .000051:
        raise ValueError('Reported total branch length differs')
    model = re.search(r'Model of substitution: (\S+)', report)[1]
    command = source['command']
    if model != 'LG+F+G4' or command[command.index('-m')+1] != model:
        raise ValueError('Guide model differs')
    likelihood = float(re.search(r'Log-likelihood of the tree: ([-\d.]+)', report)[1])
    if not math.isfinite(likelihood):
        raise ValueError('Nonfinite likelihood')
    log = (a.guide / 'guide.log').read_text()
    composition = {}
    for number, taxon, gaps, status, pvalue in re.findall(r'^\s*(\d+)\s+(\S+)\s+([\d.]+)%\s+(passed|failed)\s+([\d.]+)%\s*$', log, re.M):
        if taxon in composition:
            raise ValueError('Repeated composition row')
        composition[taxon] = {'gap_ambiguity_percent_reported': gaps,
                              'nominal_composition_status': status,
                              'composition_pvalue_percent_reported': pvalue}
    if set(composition) != expected:
        raise ValueError('Composition report taxon grid differs')
    failed = sum(r['nominal_composition_status'] == 'failed' for r in composition.values())
    if int(re.search(r'(\d+) sequences failed composition chi2 test', log)[1]) != failed:
        raise ValueError('Composition aggregate differs')
    a.output.mkdir(parents=True)
    with (a.output / 'taxon_qc.tsv').open('w') as handle:
        rows = [dict(taxon_id=r['taxon_id'], species_name=r['species_name'], study_role=r['study_role'],
                     **composition[r['taxon_id']]) for r in taxa]
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    for name in ['guide.treefile','guide.iqtree','guide.log','run_config.json']:
        shutil.copyfile(a.guide / name, a.output / name)
    result = {'status': 'passed_full_species_guide_readback', 'taxa': len(expected),
        'ingroup': sum(r['study_role']=='ingroup' for r in taxa),
        'outgroup': sum(r['study_role']=='outgroup' for r in taxa), 'columns': matrix['columns'],
        'edges': len(observed), 'total_branch_length': math.fsum(observed.values()),
        'model': model, 'log_likelihood': likelihood, 'nominal_composition_failures': failed,
        'guide_receipt_sha256': sha(a.guide / 'receipt.json'), 'matrix_receipt_sha256': sha(a.matrix / 'receipt.json'),
        'script_sha256': sha(Path(__file__)),
        'interpretation': 'Unrooted homogeneous-model guide without support estimates. Full tip identities, finite edges and report/tree agreement checked using Bio.Phylo; no independent tree inference. Nominal composition screening is not a calibrated topology test. Richer models, support, discordance and sampling/marker sensitivities remain required.',
        'artifacts': {p.name: sha(p) for p in a.output.iterdir()}}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
