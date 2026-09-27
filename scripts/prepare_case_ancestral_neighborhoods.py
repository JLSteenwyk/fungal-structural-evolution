#!/usr/bin/env python3
"""Extract guide-dependent ancestral neighborhoods for every selected case."""
import csv
import hashlib
import json
from io import StringIO
from pathlib import Path
from Bio import Phylo


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def main():
    case_root = Path('results/structural_comparisons/whole-domain-case-dossiers-20260927-v1')
    case_path = case_root/'case_dossiers.tsv'
    receipt = json.loads((case_root/'receipt.json').read_text())
    assert sha(case_path) == receipt['artifacts'][case_path.name]
    cases = read(case_path)
    assert len(cases) == 13 and len({r['family'] for r in cases}) == 13
    plan_path = Path('metadata/duplication_candidate_tree_review_plan_20260926.json')
    plan = json.loads(plan_path.read_text())
    wanted = {r['family']:r for r in cases}
    output = Path('results/ancestral/case-neighborhoods-20260927-v1')
    output.mkdir(parents=True, exist_ok=False)
    pins = {str(p):sha(p) for p in [case_path, case_root/'receipt.json', plan_path]}
    summaries, descendants, neighborhoods = [], [], {}
    for source in plan['guides']:
        path = Path(source['trees'])
        proof_path = Path(source['tree_readback'])
        proof = json.loads(proof_path.read_text())
        assert proof['status'] == 'passed_complete_resolved_tree_membership_readback'
        assert sha(path) == proof['resolved_tree_file_sha256']
        pins[str(path)] = sha(path)
        pins[str(proof_path)] = sha(proof_path)
        seen = set()
        with path.open() as handle:
            for line in handle:
                family, newick = line.rstrip().split(': ', 1)
                if family not in wanted:
                    continue
                assert family not in seen
                seen.add(family)
                case = wanted[family]
                tree = Phylo.read(StringIO(newick), 'newick')
                nodes = list(tree.find_clades())
                names = [n.name for n in nodes]
                assert all(names) and len(names) == len(set(names))
                tips = {n.name:n for n in tree.get_terminals()}
                a, b = tips[case['gene_a']], tips[case['gene_b']]
                ancestor = tree.common_ancestor(a, b)
                parents = {child:node for node in nodes for child in node.clades}
                treepath = output/(source['guide']+'-'+family+'.nwk')
                treepath.write_text(newick+'\n')
                assert {t.name for t in Phylo.read(treepath,'newick').get_terminals()} == set(tips)
                for level in range(4):
                    genes = sorted(t.name for t in ancestor.get_terminals())
                    assert case['gene_a'] in genes and case['gene_b'] in genes
                    guide = source['guide']
                    neighborhoods[guide,family,level] = set(genes)
                    summaries.append(dict(guide=guide,family=family,species_name=case['species_name'],gene_a=case['gene_a'],gene_b=case['gene_b'],level=level,node=ancestor.name,proteins=len(genes),taxa=len({g.split('_',1)[0] for g in genes}),whole_family_proteins=len(tips),is_root=int(ancestor is tree.root),focal_pair_exact=int(level==0 and set(genes)=={case['gene_a'],case['gene_b']}),direct_children=len(ancestor.clades),tree_file=treepath.name))
                    descendants.extend(dict(guide=guide,family=family,level=level,node=ancestor.name,gene=g) for g in genes)
                    if ancestor not in parents:
                        break
                    ancestor = parents[ancestor]
        assert seen == set(wanted)
    comparisons = []
    for family in wanted:
        for level in range(4):
            a = neighborhoods.get(('profile',family,level))
            b = neighborhoods.get(('mafft',family,level))
            comparisons.append(dict(family=family,level=level,profile_proteins=len(a) if a is not None else '',mafft_proteins=len(b) if b is not None else '',status='missing_level' if a is None or b is None else 'exact_descendant_set' if a==b else 'different_descendant_set',shared_proteins=len(a&b) if a is not None and b is not None else '',union_proteins=len(a|b) if a is not None and b is not None else ''))
    for filename, rows in [('neighborhoods.tsv',summaries),('descendants.tsv',descendants),('guide_comparison.tsv',comparisons)]:
        with (output/filename).open('w') as f:
            w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
        assert read(output/filename)==[{k:str(v) for k,v in r.items()} for r in rows]
    for path,h in pins.items():
        assert sha(path)==h
    result=dict(status='complete_selected_case_ancestral_neighborhood_inventory',cases=13,guide_trees=26,neighborhoods=len(summaries),descendant_records=len(descendants),comparisons=comparisons,source_hashes=pins,script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in output.iterdir()},scope='Exact extracted resolved trees and focal MRCA plus up to three ancestral levels for every selected case under both guides. Same level is descriptive and does not assert node homology. Reconciliation-derived topology is not independent duplication support. No ancestral sequence, alignment-quality, root certainty, posterior probability, species delimitation or prediction result.')
    (output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','cases','guide_trees','neighborhoods','descendant_records']}),flush=True)
    for row in comparisons:
        if row['level']<2:print(row,flush=True)


if __name__=='__main__':
    main()
