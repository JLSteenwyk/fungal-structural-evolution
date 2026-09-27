#!/usr/bin/env python3
"""Preserve all tips and enumerate every resolution of the observed four-tip star."""
import copy, itertools, json
from pathlib import Path
from Bio import Phylo, SeqIO
from Bio.Phylo.BaseTree import Clade
from prepare_case_ancestral_neighborhoods import sha


def canonical(value):
    if isinstance(value, str): return value
    return '(' + ','.join(sorted(canonical(v) for v in value)) + ')'


def resolutions(labels):
    if len(labels) == 1: return [labels[0]]
    result = []
    for size in range(1, len(labels)):
        for subset in itertools.combinations(labels[1:], size - 1):
            left = (labels[0],) + subset
            right = tuple(x for x in labels if x not in left)
            for a in resolutions(left):
                for b in resolutions(right): result.append((a, b))
    return result


def verify_resolution_set(labels, variants):
    # Independent enumeration of the two possible four-leaf rooted shapes.
    expected = set()
    for a, b, c, d in itertools.permutations(labels):
        expected.add(canonical(((a, b), (c, d))))
        expected.add(canonical((a, (b, (c, d)))))
    actual = [canonical(v) for v in variants]
    assert len(actual) == len(set(actual)) == 15
    assert set(actual) == expected


def main():
    audit_path = Path('metadata/historian_project_input_audit_20260927.json')
    audit = json.loads(audit_path.read_text())
    for p, h in audit['pins'].items(): assert sha(p) == h, p
    pp = Path('metadata/ancestral_indel_coding_plan_20260927.json')
    old = json.loads(pp.read_text())
    out = Path('results/ancestral/historian-sensitivity-inputs-20260927-v1')
    out.mkdir(parents=True, exist_ok=False)
    trees = {}; records = []; jobs = []
    for row in audit['inputs']:
        source_path = Path(row['tree'])
        if str(source_path) in trees: continue
        source = Phylo.read(source_path, 'newick')
        original = {n.name: n for n in source.find_clades()}
        nonbinary = [n for n in source.get_nonterminals() if len(n.clades) != 2]
        variants = [None]; star = None
        if nonbinary:
            assert len(nonbinary) == 1
            star = nonbinary[0]
            assert star.name == 'n1790' and len(star.clades) == 4
            assert all(c.is_terminal() and c.branch_length == 0 for c in star.clades)
            labels = tuple(sorted(c.name for c in star.clades))
            variants = sorted(resolutions(labels), key=canonical)
            verify_resolution_set(labels, variants)
        derived = []
        for floor in [1e-9, 1e-7]:
            for index, variant in enumerate(variants):
                tree = copy.deepcopy(source)
                nodes = {n.name: n for n in tree.find_clades()}
                introduced = []
                if variant is not None:
                    def build(part):
                        if isinstance(part, str): return nodes[part]
                        name = 'historian_artificial_' + str(len(introduced) + 1)
                        assert name not in nodes
                        introduced.append(name)
                        return Clade(branch_length=floor, name=name, clades=[build(c) for c in part])
                    nodes[star.name].clades = [build(c) for c in variant]
                    assert len(introduced) == 2
                for node in tree.find_clades():
                    if node is not tree.root: node.branch_length = max(node.branch_length, floor)
                name = source_path.stem + '-floor' + str(floor) + '-resolution' + str(index)
                path = out / (name + '.nwk')
                Phylo.write(tree, path, 'newick', format_branch_length='%.17g')
                check = Phylo.read(path, 'newick')
                cnodes = {n.name: n for n in check.find_clades()}
                assert check.root.name == source.root.name
                assert {n.name for n in check.get_terminals()} == {n.name for n in source.get_terminals()}
                assert all(len(n.clades) == 2 for n in check.get_nonterminals())
                assert set(cnodes) == set(original) | set(introduced)
                def visible_children(n):
                    return set().union(*(visible_children(c) if c.name in introduced else {c.name} for c in n.clades))
                for key, node in original.items():
                    assert visible_children(cnodes[key]) == {c.name for c in node.clades}
                    expected_length = node.branch_length if node is source.root else max(node.branch_length, floor)
                    assert cnodes[key].branch_length == expected_length
                # Check every root-to-tip distance against its explicitly changed edges.
                errors = []
                for tip in source.get_terminals():
                    new_path = check.get_path(cnodes[tip.name])
                    increment = sum(floor if n.name in introduced else max(0., floor-original[n.name].branch_length) for n in new_path)
                    error = abs(check.distance(cnodes[tip.name]) - source.distance(tip) - increment)
                    assert error < 1e-12
                    errors.append(error)
                record = dict(tree_id=name, source_tree=str(source_path), source_sha256=sha(source_path),
                              path=str(path), sha256=sha(path), floor=floor, resolution=index,
                              resolution_shape=canonical(variant) if variant else None,
                              artificial_nodes=introduced, root=tree.root.name,
                              maximum_path_readback_error=max(errors), proteins=len(check.get_terminals()))
                records.append(record); derived.append(record)
        trees[str(source_path)] = derived
    audit_rows = {r['input_id']: r for r in audit['inputs']}
    identical_checks = []
    for job in old['jobs']:
        row = audit_rows[job['input_id']]
        if row['nonbinary_nodes']:
            seq = {r.id: str(r.seq) for r in SeqIO.parse(job['alignment'], 'fasta')}
            source = Phylo.read(row['tree'], 'newick')
            star = next(n for n in source.find_clades() if n.name == 'n1790')
            assert len({seq[n.name] for n in star.clades}) == 1
            identical_checks.append(job['input_id'])
        for tree in trees[row['tree']]:
            jobs.append(dict(job_id=job['input_id']+'-floor'+str(tree['floor'])+'-resolution'+str(tree['resolution']),
                             input_id=job['input_id'], family=job['family'], alignment=job['alignment'],
                             alignment_sha256=job['alignment_sha256'], proteins=job['proteins'],
                             tree=tree['path'], tree_sha256=tree['sha256'], root=tree['root'],
                             artificial_nodes=tree['artificial_nodes']))
    assert len(records) == 108 and len(jobs) == 324 and len(identical_checks) == 6
    result = dict(status='complete_full_input_sensitivity_grid', trees=records, jobs=jobs,
                  identical_star_alignment_checks=identical_checks,
                  pins={str(p):sha(p) for p in [audit_path, pp, Path(__file__)]},
                  scope='All78 alignments and all tips retained; two explicit zero-edge floors and all15 resolutions of the identical four-tip star. These are computational sensitivities, not supported new branches or posterior-weighted histories. Original roots retained as assumptions. No reconstruction yet.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(len(records),'trees',len(jobs),'jobs prepared')

if __name__ == '__main__': main()
