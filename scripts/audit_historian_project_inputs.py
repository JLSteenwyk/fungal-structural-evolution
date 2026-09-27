#!/usr/bin/env python3
"""Inventory every selected alignment against pinned Historian input requirements."""
import json
from pathlib import Path
from Bio import Phylo, SeqIO
from prepare_case_ancestral_neighborhoods import sha


def main():
    source_plan = Path('metadata/ancestral_indel_coding_plan_20260927.json')
    source = Path('data/software_audits/historian-20260927/source')
    plan = json.loads(source_plan.read_text())
    tree_dir = Path('results/ancestral/case-local-trees-20260927-v1')
    tree_receipt = json.loads((tree_dir / 'receipt.json').read_text())
    pins = {str(p): sha(p) for p in [source_plan, tree_dir / 'receipt.json',
            Path(__file__), source / 'src/tree.cpp', source / 'src/tree.h',
            source / 'src/recon.cpp', source / 'src/sumprod.h']}
    rows = []
    for job in plan['jobs']:
        alignment = Path(job['alignment'])
        assert sha(alignment) == job['alignment_sha256']
        dataset = 'whole' if job['boundary'] == 'whole' else 'domain'
        path = tree_dir / ('profile-' + job['family'] + '-' + dataset + '.nwk')
        assert sha(path) == tree_receipt['artifacts'][path.name]
        pins[str(path)] = sha(path)
        pins[str(alignment)] = sha(alignment)
        records = list(SeqIO.parse(alignment, 'fasta'))
        seq = {r.id: str(r.seq) for r in records}
        tree = Phylo.read(path, 'newick')
        tips = tree.get_terminals()
        assert len(records) == len(seq) == len(tips) == job['proteins']
        assert set(seq) == {n.name for n in tips}
        assert {len(s) for s in seq.values()} == {job['columns']}
        nodes = list(tree.find_clades())
        assert len({n.name for n in nodes}) == len(nodes)
        nonbinary = [dict(node=n.name, children=len(n.clades), root=n is tree.root)
                     for n in tree.get_nonterminals() if len(n.clades) != 2]
        floor = [dict(node=n.name, source_length=n.branch_length)
                 for n in nodes if n is not tree.root and n.branch_length < 1e-9]
        rows.append(dict(input_id=job['input_id'], family=job['family'],
                         proteins=len(seq), columns=job['columns'], tree=str(path),
                         root=tree.root.name, nonbinary_nodes=nonbinary,
                         parser_branch_floor_changes=floor,
                         unknown_X_residues=sum(s.count('X') for s in seq.values()),
                         maximum_ungapped_length=max(len(s.replace('-', '')) for s in seq.values()),
                         tip_root_distance_range=[min(tree.depths()[n] for n in tips),
                                                  max(tree.depths()[n] for n in tips)],
                         status='requires_polytomy_treatment' if nonbinary else 'binary_input_compatible_capacity_untested'))
    assert len(rows) == 78
    out = dict(status='complete_full_78_input_inventory_not_inference', inputs=rows,
               counts=dict(alignments=len(rows), families=len({r['family'] for r in rows}),
                           maximum_proteins=max(r['proteins'] for r in rows),
                           nonbinary_inputs=sum(bool(r['nonbinary_nodes']) for r in rows),
                           branch_floor_inputs=sum(bool(r['parser_branch_floor_changes']) for r in rows)),
               pins=pins,
               scope='All selected alignments retained. No tree resolution, branch flooring, sequence modification or inference performed. Binary topology is only an input precondition; rooted history, model fitting, full capacity and uncertainty remain unqualified.')
    target = Path('metadata/historian_project_input_audit_20260927.json')
    target.write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out['counts']))


if __name__ == '__main__':
    main()
