#!/usr/bin/env python3
"""Reproduce a root-HOG leaf omission in the installed native writer, without files."""
import hashlib
import json
from pathlib import Path
from orthofinder.gene_tree_inference import trees2ologs_of as native


def run(duplication):
    species = native.tree_lib.Tree('(0,1)N0;', format=1)
    tree = native.tree_lib.Tree('(0_a,((0_b,1_b)n2,(0_c,1_c)n3)n1)n0;', format=1)
    writer = native.HogWriter(species, ['N0'], {n.name:n.name for n in tree},
                              {'0':'s0','1':'s1'}, [0,1], False, write_output=False)
    for node in tree.traverse():
        if not node.is_leaf():
            node.add_feature('sp_node', 'N0')
            node.add_feature('dup', duplication and node.name == 'n1')
    writer.mark_dups_below(tree)
    rows = []
    for node in tree.traverse('preorder'):
        rows.extend(row for level,row in writer.write_clade_v2(node,'OGfixture') if level == 'N0')
    genes = {g for row in rows for cell in row[2:] for g in cell.split(', ') if g}
    return rows, genes


def main():
    rows, genes = run(True)
    assert genes == {'0_b','1_b','0_c','1_c'}
    assert {r[1] for r in rows} == {'n2','n3'}
    control, all_genes = run(False)
    assert all_genes == genes | {'0_a'}
    assert len(control) == 1 and control[0][1] == 'n0'
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    print(json.dumps({'status':'passed_native_hog_leaf_omission_reproduction',
        'native_module':str(Path(native.__file__).resolve()),
        'native_module_sha256':digest(native.__file__),
        'script_sha256':digest(__file__),
        'with_nested_root_duplication':{'parent_clades':[r[1] for r in rows], 'omitted':['0_a']},
        'without_duplication':{'parent_clades':[r[1] for r in control], 'omitted':[]},
        'scope':'Synthetic assigned reconciliation features exercise the installed HOG writer. Demonstrates a possible omission mechanism and counterfactual, not independent reconciliation or proof of the cause of every production omission. No native outputs written.'}, indent=2))


if __name__ == '__main__':
    main()
