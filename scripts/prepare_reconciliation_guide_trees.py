#!/usr/bin/env python3
"""Root audited guide sensitivities and preserve exact native species-ID mapping."""
import argparse
import copy
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
from Bio import Phylo

ROOT = Path(__file__).resolve().parents[1]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def edges(tree):
    universe = frozenset(t.name for t in tree.get_terminals())
    result = {}
    for node in tree.find_clades():
        if node is tree.root:
            continue
        side = frozenset(t.name for t in node.get_terminals())
        key = min([tuple(sorted(side)), tuple(sorted(universe - side))], key=lambda x: (len(x), x))
        if node.branch_length is None or not math.isfinite(node.branch_length) or node.branch_length < 0:
            raise ValueError('Invalid edge length')
        result[key] = result.get(key, 0.) + node.branch_length
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    a = ap.parse_args()
    plan = json.loads(a.plan.read_text())
    for name, expected in plan['pins'].items():
        if sha(Path(name)) != expected:
            raise ValueError('Changed pinned input ' + name)
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    with Path(plan['manifest']).open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    roles = {r['taxon_id']: r['study_role'] for r in rows}
    if len(roles) != len(rows) or set(roles.values()) != {'ingroup', 'outgroup'}:
        raise ValueError('Invalid taxon-role universe')
    groups = {role: frozenset(t for t, r in roles.items() if r == role) for role in ['ingroup', 'outgroup']}
    if len(groups['ingroup']) != 501 or len(groups['outgroup']) != 25:
        raise ValueError('Expected 501/25 taxon-entry design')
    native_to_taxon = {}
    for line in Path(plan['species_ids']).read_text().splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        native, filename = line.split(': ', 1)
        taxon = filename.rsplit('.', 1)[0]
        if native in native_to_taxon:
            raise ValueError('Duplicate native ID')
        native_to_taxon[native] = taxon
    if len(set(native_to_taxon.values())) != 526 or set(native_to_taxon.values()) != set(roles):
        raise ValueError('Species ID/taxon mismatch')
    taxon_to_native = {v: k for k, v in native_to_taxon.items()}
    (out / 'species_id_mapping.json').write_text(json.dumps(native_to_taxon, indent=2) + '\n')
    results = []
    for name, guide in plan['guides'].items():
        path = Path(guide)
        source = Phylo.read(path, 'newick')
        tips = [t.name for t in source.get_terminals()]
        if len(tips) != 526 or set(tips) != set(roles):
            raise ValueError('Guide taxon universe mismatch')
        original = edges(source)
        candidates = [n for n in source.find_clades() if n is not source.root
                      and frozenset(t.name for t in n.get_terminals()) in groups.values()]
        if len(candidates) != 1:
            raise ValueError('Require one explicit fungal/outgroup separating edge')
        edge_length = candidates[0].branch_length
        source.root_with_outgroup(candidates[0], outgroup_branch_length=edge_length / 2)
        if {frozenset(t.name for t in n.get_terminals()) for n in source.root.clades} != set(groups.values()):
            raise ValueError('Root partition differs from required role split')
        target = out / name
        target.mkdir()
        named = target / 'species_tree_taxa.nwk'
        Phylo.write(source, named, 'newick', format_branch_length='%.15g')
        rooted = Phylo.read(named, 'newick')
        current = edges(rooted)
        if current.keys() != original.keys() or any(abs(current[k] - original[k]) > 1e-12 for k in original):
            raise ValueError('Rooting changed original unrooted edge lengths or topology')
        native_tree = copy.deepcopy(rooted)
        for tip in native_tree.get_terminals():
            tip.name = taxon_to_native[tip.name]
        native_path = target / 'species_tree_ids.nwk'
        Phylo.write(native_tree, native_path, 'newick', format_branch_length='%.15g')
        # Run the actual native converter as an independent mapping check. Its
        # default serializer rounds branch lengths, so preserve our high-precision
        # IDs file as the intended input and report native round-trip error.
        native_check = target / 'native_converter_check.nwk'
        code = ('import json,sys; from orthofinder.gene_tree_inference.infer_trees import ConvertUserSpeciesTree; '
                'ConvertUserSpeciesTree(sys.argv[1],json.load(open(sys.argv[2])),sys.argv[3])')
        subprocess.run([plan['native_python'], '-c', code, str(named.resolve()),
                        str((out / 'species_id_mapping.json').resolve()), str(native_check.resolve())], check=True)
        check = Phylo.read(native_check, 'newick')
        check_tips = [t.name for t in check.get_terminals()]
        if len(check_tips) != 526 or set(check_tips) != set(native_to_taxon):
            raise ValueError('Native converter taxon universe differs')
        for tip in check.get_terminals():
            tip.name = native_to_taxon[tip.name]
        native_edges = edges(check)
        if native_edges.keys() != original.keys():
            raise ValueError('Native conversion changes topology')
        errors = [abs(native_edges[k] - original[k]) for k in original]
        if any(abs(native_edges[k] - original[k]) > 1e-5 * max(original[k], 1e-6) for k in original):
            raise ValueError('Native conversion length discrepancy exceeds serialization tolerance')
        results.append(dict(guide=name, tips=526, unrooted_edges=len(original),
                            root_split_taxa=[501, 25], root_edge_original_length=edge_length,
                            root_child_lengths=[n.branch_length for n in rooted.root.clades],
                            maximum_native_converter_edge_rounding_error=max(errors),
                            artifacts={p.name: sha(p) for p in target.iterdir()}))
    for name, expected in plan['pins'].items():
        if sha(Path(name)) != expected:
            raise ValueError('Pinned input changed during preparation')
    receipt = dict(status='prepared_rooted_reconciliation_guide_sensitivities', guides=results,
                   plan_sha256=sha(a.plan), script_sha256=sha(Path(__file__)),
                   mapping_sha256=sha(out / 'species_id_mapping.json'),
                   interpretation='Both homogeneous guide alternatives rooted on their shared 501-fungal/25-outgroup separating edge. Equal edge halves are a serialization convention, not inferred root position or branch duration. Unrooted topology/lengths and 526 native IDs checked. No final supported topology, time calibration or completed reconciliation is claimed.')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
