#!/usr/bin/env python3
"""Read back PMSF profiles, trees, report and empirical bootstrap split support."""
import argparse
from collections import Counter
import csv
from decimal import Decimal
import io
import json
import math
from pathlib import Path
import re
from Bio import Phylo, SeqIO
from audit_busco_gene_copies import sha


def canonical(mask, all_mask):
    other = all_mask ^ mask
    return min((mask, other), key=lambda x: (x.bit_count(), x))


def tree_edges(parsed, taxa):
    names = [n.name for n in parsed.get_terminals()]
    if len(names) != len(taxa) or set(names) != set(taxa):
        raise ValueError('Tree tip universe differs')
    index = {name: 1 << i for i, name in enumerate(taxa)}
    all_mask = (1 << len(taxa)) - 1
    masks, edges = {}, {}
    for node in parsed.find_clades(order='postorder'):
        if node.is_terminal():
            mask = index[node.name]
        else:
            mask = 0
            for child in node.clades:
                mask |= masks[child]
        masks[node] = mask
        if node is parsed.root:
            if node.branch_length not in (None, 0):
                raise ValueError('Unexpected root stem')
            continue
        if node.branch_length is None or not math.isfinite(node.branch_length) or node.branch_length < 0:
            raise ValueError('Missing, negative or nonfinite branch length')
        key = canonical(mask, all_mask)
        if key in edges:
            raise ValueError('Repeated unrooted split')
        edges[key] = (node.branch_length, node.name, node.confidence)
    return edges


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new immutable audit directory')
    source = json.loads((args.run / 'receipt.json').read_text())
    config = json.loads((args.run / 'config.json').read_text())
    if source['status'] != 'complete_pmsf_execution_pending_full_audit' or source['returncode'] != 0:
        raise ValueError('Completed inference required')

    def verify():
        if sha(args.run / 'config.json') != source['config_sha256']:
            raise ValueError('Changed execution config')
        for name, digest in source['artifacts'].items():
            if sha(args.run / name) != digest:
                raise ValueError('Changed output: ' + name)
        for name, digest in config['pinned_files'].items():
            if sha(Path(name)) != digest:
                raise ValueError('Changed input: ' + name)

    verify()
    command = config['command']
    for flag, expected in [('-m', 'LG+C20+F+G4'), ('--alrt', '1000'), ('-B', '1000')]:
        if command[command.index(flag) + 1] != expected:
            raise ValueError('Unexpected inference settings')
    if '--bnni' not in command or '--boot-trees' not in command or '--tree-freq' not in command:
        raise ValueError('Required PMSF/bootstrap settings missing')
    alignment = Path(command[command.index('-s') + 1])
    records = list(SeqIO.parse(alignment, 'fasta'))
    taxa = sorted(r.id for r in records)
    if (len(taxa) != 526 or len(set(taxa)) != 526
            or any(len(r.seq) != source['columns'] for r in records)):
        raise ValueError('Alignment universe differs')
    all_mask = (1 << len(taxa)) - 1
    profiles = 0
    max_sum_error = 0.0
    with (args.run / 'pmsf.sitefreq').open() as handle:
        for site, line in enumerate(handle, 1):
            fields = line.split()
            if len(fields) != 21 or int(fields[0]) != site:
                raise ValueError('Site-profile dimensions or order differ')
            values = list(map(float, fields[1:]))
            if any(not math.isfinite(x) or x < 0 or x > 1 for x in values):
                raise ValueError('Invalid amino-acid profile')
            error = abs(math.fsum(values) - 1)
            if error > 1e-5:
                raise ValueError('Profile frequencies do not sum to one at output precision')
            max_sum_error = max(max_sum_error, error)
            profiles += 1
    if profiles != source['columns']:
        raise ValueError('Missing site profiles')
    report = (args.run / 'pmsf.iqtree').read_text()
    log = (args.run / 'pmsf.log').read_text()
    if (f'Input data: 526 sequences with {profiles} amino-acid sites' not in report
            or 'Model of substitution: LG+SSF+F+G4' not in report
            or 'SH-aLRT support (%) / ultrafast bootstrap support (%)' not in report
            or 'Computing posterior mean site frequencies' not in log):
        raise ValueError('Model or report dimensions differ')
    trees = {name: tree_edges(Phylo.read(args.run / filename, 'newick'), taxa)
             for name, filename in [('ml', 'pmsf.treefile'), ('consensus', 'pmsf.contree'),
                                    ('guide', 'input_guide.treefile')]}
    if len(trees['ml']) != 2 * len(taxa) - 3:
        raise ValueError('Nonbinary ML tree')
    reported = [line.strip() for line in report.splitlines()
                if line.startswith('(') and line.rstrip().endswith(';')]
    if len(reported) != 2:
        raise ValueError('Unexpected number of reported trees')
    for name, text in zip(['ml', 'consensus'], reported):
        if tree_edges(Phylo.read(io.StringIO(text), 'newick'), taxa) != trees[name]:
            raise ValueError('Report and saved tree differ: ' + name)
    counts = Counter()
    replicates = 0
    with (args.run / 'pmsf.ufboot').open() as handle:
        for parsed in Phylo.parse(handle, 'newick'):
            edges = tree_edges(parsed, taxa)
            counts.update(edges.keys())
            replicates += 1
    if replicates != 1000:
        raise ValueError('Bootstrap count differs')
    nexus = (args.run / 'pmsf.splits.nex').read_text()
    labels = {int(i): name for i, name in re.findall(r"^\[(\d+)\] '([^']+)'$", nexus, re.M)}
    if set(labels) != set(range(1, 527)) or set(labels.values()) != set(taxa):
        raise ValueError('NEXUS taxon mapping differs')
    index = {name: 1 << i for i, name in enumerate(taxa)}
    matrix = re.search(r'\bMATRIX\s*\n(.*?)\n\s*;', nexus, re.S)[1]
    nexus_splits = {}
    nexus_precision = {}
    for line in matrix.splitlines():
        parts = line.strip().rstrip(',').split()
        if not parts:
            continue
        weight = float(parts[0])
        members = list(map(int, parts[1:]))
        if len(members) != len(set(members)) or not set(members) <= set(labels):
            raise ValueError('Malformed NEXUS split')
        mask = 0
        for i in members:
            mask |= index[labels[i]]
        key = canonical(mask, all_mask)
        if key in nexus_splits or not math.isfinite(weight):
            raise ValueError('Duplicate or nonfinite NEXUS split')
        nexus_splits[key] = weight
        nexus_precision[key] = .5 * 10 ** Decimal(parts[0]).as_tuple().exponent
    if set(nexus_splits) != set(counts):
        raise ValueError('NEXUS and empirical bootstrap split universes differ')
    differences = [(weight, 100 * counts[key] / replicates) for key, weight in nexus_splits.items()
                   if abs(weight - 100 * counts[key] / replicates) > nexus_precision[key] + 1e-8]
    if differences:
        raise ValueError('NEXUS weights differ from empirical bootstrap counts: ' + repr(differences[:10]))
    rows = []
    for name in ['ml', 'consensus']:
        for mask, (length, label, confidence) in trees[name].items():
            if mask.bit_count() == 1:
                continue
            if name == 'ml':
                alrt, boot = map(float, (label or '').split('/'))
                if not math.isfinite(alrt) or not 0 <= alrt <= 100:
                    raise ValueError('Invalid SH-aLRT label')
            else:
                alrt, boot = None, confidence
            empirical = 100 * counts[mask] / replicates
            if boot is None or not math.isfinite(boot) or not 0 <= boot <= 100 or abs(boot - empirical) > .500001:
                raise ValueError('Tree support differs from rounded empirical bootstrap frequency')
            rows.append(dict(tree=name, split_taxa_json=json.dumps([taxon for i, taxon in enumerate(taxa) if mask & (1 << i)]),
                             branch_length=length, sh_alrt_percent=alrt,
                             reported_ufboot_percent=boot, empirical_ufboot_percent=empirical))
    rf = len(set(trees['ml']) ^ set(trees['consensus']))
    reported_rf = int(re.search(r'Robinson-Foulds distance between ML tree and consensus tree: (\d+)', report)[1])
    if rf != reported_rf:
        raise ValueError('Reported topology difference disagrees')
    ml_lnl = float(re.search(r'Log-likelihood of the tree: ([-\d.]+)', report)[1])
    consensus_lnl = float(re.search(r'Log-likelihood of consensus tree: ([-\d.]+)', report)[1])
    if not all(math.isfinite(x) for x in [ml_lnl, consensus_lnl]):
        raise ValueError('Nonfinite reported likelihood')
    total = float(re.search(r'Total tree length \(sum of branch lengths\): ([\d.]+)', report)[1])
    if abs(math.fsum(x[0] for x in trees['ml'].values()) - total) > .000051:
        raise ValueError('ML total branch length differs')
    verify()
    args.output.mkdir(parents=True)
    with (args.output / 'branch_support.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    result = dict(status='passed_pmsf_profile_tree_and_bootstrap_readback',
                  taxa=len(taxa), sites=profiles, bootstrap_trees=replicates,
                  empirical_splits=len(counts), internal_support_rows=len(rows),
                  nexus_support_rounding_tolerance_max=max(nexus_precision.values()),
                  max_profile_sum_rounding_error=max_sum_error, ml_consensus_rf=rf,
                  reported_ml_log_likelihood=ml_lnl, reported_consensus_log_likelihood=consensus_lnl,
                  consensus_likelihood_improvement=consensus_lnl - ml_lnl,
                  warnings=[line for line in log.splitlines() if 'WARNING' in line],
                  source_receipt_sha256=sha(args.run / 'receipt.json'),
                  script_sha256=sha(Path(__file__)),
                  artifacts={'branch_support.tsv': sha(args.output / 'branch_support.tsv')},
                  scope='All saved profiles, report trees, 1000 bootstrap tip/branch grids and empirical split supports checked. SH-aLRT labels validated but not independently recomputed; likelihoods read back, not recomputed. Guide dependence, model adequacy, rooting and discordance remain separate. Not the final species phylogeny.')
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
