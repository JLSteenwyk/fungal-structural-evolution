#!/usr/bin/env python3
"""Audit every reported duplication's identities and clade; not event inference."""
import argparse
from collections import Counter
import csv
import json
import math
from pathlib import Path
from Bio import Phylo
from assess_small_family_output_exposure import groups, sha

HEADER = ['Orthogroup', 'Species Tree Node', 'Gene Tree Node', 'Support', 'Type', 'Genes 1', 'Genes 2']


def check_event(row, lookup, descendants, terminal):
    if len(row) != 7:
        raise ValueError('Malformed event')
    family, node, gene_node, support, kind, left, right = row
    if node not in descendants or not gene_node:
        raise ValueError('Unknown species node or empty gene node')
    value = float(support)
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError('Invalid support')
    allowed = {'Terminal'} if node in terminal else {'Non-Terminal', 'Non-Terminal: STRIDE'}
    if kind not in allowed:
        raise ValueError('Type disagrees with species node')
    # Native multifurcations can emit all descendants in Genes 1 and empty Genes 2.
    sides = [cell.split(', ') if cell else [] for cell in (left, right)]
    seen = set()
    for side in sides:
        for label in side:
            if label in seen:
                raise ValueError('Repeated gene within event')
            seen.add(label)
            if label not in lookup:
                raise ValueError('Unknown protein')
            taxon, expected_family = lookup[label]
            if expected_family != family:
                raise ValueError('Wrong family')
            if taxon not in descendants[node]:
                raise ValueError('Gene outside species clade')
    if not sides[0] or len(seen) < 2:
        raise ValueError('Insufficient event genes')
    return len(seen), not bool(sides[1])


def audit(entry, output):
    source, result = Path(entry['source']), Path(entry['result'])
    species = {}
    for line in (source/'SpeciesIDs.txt').read_text().splitlines():
        native, name = line.split(': ', 1)
        if native in species:
            raise ValueError('Repeated species')
        species[native] = name.rsplit('.', 1)[0]
    if len(set(species.values())) != len(species):
        raise ValueError('Ambiguous species')
    membership = {}
    families = set()
    for family, genes in groups(source/'clusters_OrthoFinder.txt_id_pairs.txt'):
        families.add(family)
        for gene in genes:
            if gene in membership:
                raise ValueError('Repeated source gene')
            membership[gene] = family
    source_genes = len(membership)
    lookup = {}
    with (source/'SequenceIDs.txt').open() as f:
        for line in f:
            gene, label = line.rstrip('\n').split(': ', 1)
            if label in lookup or len(label.split()) != 1:
                raise ValueError('Ambiguous protein label')
            lookup[label] = (species[gene.split('_')[0]], membership.pop(gene))
    if membership:
        raise ValueError('Unmapped source genes')
    tree = Phylo.read(result/'Species_Tree/SpeciesTree_rooted_node_labels.txt', 'newick')
    descendants = {}
    terminal = set()
    for node in tree.find_clades():
        if not node.name or node.name in descendants:
            raise ValueError('Ambiguous tree label')
        descendants[node.name] = {tip.name for tip in node.get_terminals()}
        if node.is_terminal():
            terminal.add(node.name)
    if terminal != set(species.values()):
        raise ValueError('Species tree coverage differs')
    counts = Counter(); by_node = Counter(); event_keys = set(); event_families = set()
    path = result/'Gene_Duplication_Events/Duplications.tsv'
    with path.open() as f:
        reader = csv.reader(f, delimiter='\t')
        if next(reader) != HEADER:
            raise ValueError('Unexpected event header')
        for row in reader:
            n, empty = check_event(row, lookup, descendants, terminal)
            key = (row[0], row[2])
            if key in event_keys:
                raise ValueError('Duplicate family/node event')
            event_keys.add(key); event_families.add(row[0])
            counts['events'] += 1; counts['gene_assignments_across_events'] += n
            counts['empty_second_side_events'] += int(empty)
            counts[row[4]] += 1; by_node[(row[1], row[4])] += 1
    table = output/(entry['guide']+'_node_counts.tsv')
    with table.open('w') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['species_tree_node', 'event_type', 'reported_events'])
        for (node, kind), n in sorted(by_node.items()):
            w.writerow([node, kind, n])
    return dict(guide=entry['guide'], source_genes=source_genes, source_families=len(families),
                event_families=len(event_families), counts=dict(counts),
                artifact=str(table), artifact_sha256=sha(table))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    a = p.parse_args(); plan = json.loads(a.plan.read_text()); plan_hash = sha(a.plan)
    def verify():
        if sha(a.plan) != plan_hash:
            raise ValueError('Changed plan')
        for path, digest in plan['pins'].items():
            if sha(path) != digest:
                raise ValueError('Changed source: '+path)
    verify()
    out = Path(plan['output']); out.mkdir(parents=True, exist_ok=False)
    csv.field_size_limit(32*1024*1024)
    summaries = []
    for entry in plan['guides']:
        (out/'state.json').write_text(json.dumps({'status':'auditing', 'guide':entry['guide']}))
        summaries.append(audit(entry, out))
    verify()
    receipt = dict(status='passed_full_duplication_event_identity_audit', plan_sha256=plan_hash,
                   guides=summaries, scope='All reported event genes checked for source identity, family, species-clade containment and uniqueness within event; finite bounded support, node type and unique family/gene-node event keys checked. Empty second sides retained as native multifurcation representation. Does not validate gene-tree node membership, event completeness, support calculation, rooting, timing, biological duplication/loss or independence of events. Counts are nested reported events, not independent observations.')
    (out/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt), flush=True)


if __name__ == '__main__':
    main()
