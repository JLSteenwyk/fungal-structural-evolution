"""Describe all aphelid BUSCO-copy placements in hash-validated retained trees.

Unrooted split exclusivity is descriptive, not evidence of duplication timing.
Missing/revised families remain explicit rather than borrowing another tree.
"""
import csv
import itertools
import json
from collections import Counter, defaultdict
from pathlib import Path
from Bio import Phylo
from readback_whole_proteome_catalog import sha
from ancestral_chain_attempt import write_json


def read(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write(path, rows):
    with open(path, 'w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    assert read(path) == [{k: str(v) for k, v in row.items()} for row in rows]


def main():
    base = Path('results/ecology/aphelid-marker-copy-inventory-20260928-v1')
    receipt = json.loads((base / 'receipt.json').read_text())
    for name, digest in receipt['artifacts'].items():
        assert sha(base / name) == digest
    proteins = read(base / 'marker_protein_copies.tsv')
    assignments = read(base / 'family_assignments.tsv')
    native = {r['protein_id']: r['native_gene_id'] for r in proteins}
    assert len(native) == 241
    catalog = Path('results/orthology/retained-tree-catalog-v1/retained_tree_catalog.tsv')
    catalog_receipt = json.loads((catalog.parent / 'receipt.json').read_text())
    assert catalog_receipt['status'] == 'completed_retained_tree_catalog'
    assert sha(catalog) == catalog_receipt['catalog_sha256']
    # Bind the catalog itself; each tree receives its own checksum verification.
    index = {}
    for row in read(catalog):
        for guide in ('mafft', 'profile'):
            key = (guide, row[guide + '_family'])
            assert key not in index
            index[key] = row
    grouped = defaultdict(list)
    for row in assignments:
        assert row['status'] == 'assigned'
        grouped[row['marker'], row['guide'], row['family']].append(row['protein_id'])
    trees = {}
    outputs, pairs = [], []
    for (marker, guide, family), ids in sorted(grouped.items()):
        row = index.get((guide, family))
        result = dict(marker=marker, guide=guide, family=family,
                      copy_count=len(ids), protein_ids_json=json.dumps(sorted(ids)),
                      status='family_not_in_retained_catalog', tree_sha256='',
                      tree_tips='', exclusive_unrooted_split='')
        if row is not None:
            result['status'] = row['status']
            if row['status'] == 'validated_exact_membership':
                path = row['tree_path']
                if path not in trees:
                    assert sha(path) == row['tree_sha256']
                    tree = Phylo.read(path, 'newick')
                    leaves = tree.get_terminals()
                    names = {leaf.name for leaf in leaves}
                    assert len(names) == len(leaves) == int(row['proteins'])
                    descendants = {}
                    for clade in tree.find_clades(order='postorder'):
                        descendants[clade] = ({clade.name} if clade.is_terminal()
                                              else set().union(*(descendants[c] for c in clade.clades)))
                    trees[path] = tree, names, descendants
                tree, names, descendants = trees[path]
                selected = {native[p] for p in ids}
                assert len(selected) == len(ids) and selected <= names
                exclusive = len(ids) > 1 and any(
                    tips == selected or names - tips == selected
                    for clade, tips in descendants.items() if clade is not tree.root)
                result.update(status='all_copies_found_in_hash_validated_tree',
                              tree_sha256=row['tree_sha256'], tree_tips=len(names),
                              exclusive_unrooted_split=str(exclusive) if len(ids) > 1 else 'not_applicable_single_copy')
                for a, b in itertools.combinations(sorted(ids), 2):
                    pairs.append(dict(marker=marker, guide=guide, family=family,
                                      protein_a=a, protein_b=b,
                                      patristic_distance=tree.distance(native[a], native[b]),
                                      tree_sha256=row['tree_sha256']))
        outputs.append(result)
    # Label-invariant agreement of grouping among all 241 target proteins.
    mapping = {g: {r['protein_id']: r['family'] for r in assignments if r['guide'] == g}
               for g in ('mafft', 'profile')}
    disagreements = sum((mapping['mafft'][a] == mapping['mafft'][b]) !=
                        (mapping['profile'][a] == mapping['profile'][b])
                        for a, b in itertools.combinations(sorted(native), 2))
    out = Path('results/ecology/aphelid-copy-tree-placements-20260928-v2')
    out.mkdir(parents=True, exist_ok=False)
    write(out / 'marker_family_placements.tsv', outputs)
    if pairs:
        write(out / 'within_family_copy_distances.tsv', pairs)
    summary = dict(status='complete_retained_tree_copy_placement_screen',
                   source_inventory_receipt_sha256=sha(base / 'receipt.json'),
                   source_catalog_sha256=sha(catalog),
                   source_catalog_receipt_sha256=sha(catalog.parent / 'receipt.json'),
                   script_sha256=sha(__file__), target_proteins=len(native),
                   target_pair_partition_comparisons=241 * 240 // 2,
                   target_pair_partition_disagreements=disagreements,
                   unique_trees_checked=len(trees), marker_guide_family_rows=len(outputs),
                   dispositions=dict(Counter(r['status'] for r in outputs)),
                   within_family_pair_rows=len(pairs),
                   exclusive_split_counts=dict(Counter(r['exclusive_unrooted_split'] for r in outputs)),
                   artifacts={p.name: sha(p) for p in out.iterdir()},
                   limitations='Retained trees only. Target-subset partition agreement does not establish full-family agreement. Split exclusivity is unrooted and descriptive; branch distances use source-tree units. No duplication timing, orthology, hybrid origin or structural divergence inference.')
    write_json(out / 'receipt.json', summary)
    write_json(Path('metadata/aphelid_copy_tree_placements_20260928.json'), summary)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
