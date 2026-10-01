#!/usr/bin/env python3
"""Full taxon character coverage and nearest declared-role boundary diagnostics."""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from Bio import SeqIO
import dendropy
from four_run_pmsf_sources import load_sources
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha

AA = 'ACDEFGHIKLMNPQRSTVWY'
CHAR_FIELDS = ['matrix', 'taxon_id', 'species_name', 'study_role', 'lineage', 'sites',
               'canonical_aa', 'gaps', 'noncanonical_or_ambiguous', 'canonical_fraction',
               'gap_fraction', 'noncanonical_fraction', 'other_characters_json']
ROLE_FIELDS = ['view', 'minimum_symmetric_difference_to_declared_outgroups', 'nearest_boundary_taxa_json',
               'missing_outgroups_json', 'added_ingroup_json', 'missing_species_json', 'added_species_json',
               'branch_length', 'sh_alrt', 'empirical_ufb']


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', required=True, type=Path); args = p.parse_args()
    plan = json.loads(args.plan.read_text()); source_plan = json.loads(Path(plan['source_plan']).read_text()); c = json.loads(Path(plan['completion']).read_text())
    assert c['status'] == 'complete_verified_full_four_run_pmsf_ML_consensus_sensitivity' and len(c['services']) == 2
    assert c['source_hashes'][plan['source_plan']] == sha(plan['source_plan'])
    sources, manifest, bindings = load_sources(source_plan, plan['source_plan'])
    for path, digest in c['source_hashes'].items(): bind(bindings, path, digest)
    for path, digest in plan['pins'].items(): bind(bindings, path, digest)
    bind(bindings, args.plan); bind(bindings, plan['completion'])
    taxa = {row['taxon_id']: row for row in manifest}; outgroups = {k for k, v in taxa.items() if v['study_role'] == 'outgroup'}
    matrices = {}
    for label, spec in source_plan['runs'].items():
        config = json.loads((Path(spec['run']) / 'config.json').read_text()); command = config['command']; path = command[command.index('-s') + 1]
        matrix = label.split('_')[0]; sites = sources[label]['audit']['sites']
        if matrix in matrices: assert matrices[matrix] == (path, sites)
        else: matrices[matrix] = path, sites
    assert set(matrices) == {'profile', 'mafft'}
    characters = []
    for matrix, (path, sites) in matrices.items():
        seen = set()
        for record in SeqIO.parse(path, 'fasta'):
            assert record.id in taxa and record.id not in seen; seen.add(record.id)
            seq = str(record.seq); assert len(seq) == sites and seq == seq.upper(); counts = Counter(seq)
            observed = sum(counts[a] for a in AA); gap = counts['-']; other = {a: count for a, count in sorted(counts.items()) if a not in AA + '-'}
            remaining = sites - observed - gap; assert sum(other.values()) == remaining
            source = taxa[record.id]
            characters.append(dict(matrix=matrix, taxon_id=record.id, species_name=source['species_name'], study_role=source['study_role'], lineage=source['lineage'], sites=sites,
                                   canonical_aa=observed, gaps=gap, noncanonical_or_ambiguous=remaining, canonical_fraction=observed / sites,
                                   gap_fraction=gap / sites, noncanonical_fraction=remaining / sites, other_characters_json=json.dumps(other, sort_keys=True)))
        assert seen == set(taxa)
    root = Path(source_plan['output'])
    with (root / 'split_presence.tsv').open() as f: presence = list(csv.DictReader(f, delimiter='\t'))
    views = list(dict.fromkeys(r['view'] for r in presence)); role_rows = []
    for view in views:
        choices = [r for r in presence if r['view'] == view and r['present'] == 'True']; assert len(choices) == 523
        minimum = min(len(set(json.loads(r['split_taxa_json'])) ^ outgroups) for r in choices)
        for row in choices:
            side = set(json.loads(row['split_taxa_json']))
            if len(side ^ outgroups) != minimum: continue
            missing, added = sorted(outgroups - side), sorted(side - outgroups)
            role_rows.append(dict(view=view, minimum_symmetric_difference_to_declared_outgroups=minimum, nearest_boundary_taxa_json=row['split_taxa_json'],
                                  missing_outgroups_json=json.dumps(missing), added_ingroup_json=json.dumps(added),
                                  missing_species_json=json.dumps([taxa[k]['species_name'] for k in missing]), added_species_json=json.dumps([taxa[k]['species_name'] for k in added]),
                                  branch_length=row['branch_length'], sh_alrt=row['sh_alrt'], empirical_ufb=row['empirical_ufb']))
    output = Path(plan['output']); output.mkdir(exist_ok=False, parents=True)
    for name, rows, fields in [('taxon_character_coverage.tsv', characters, CHAR_FIELDS), ('nearest_role_boundary.tsv', role_rows, ROLE_FIELDS)]:
        with (output / name).open('x') as f:
            writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n'); writer.writeheader(); writer.writerows(rows)
    # Independently parse complete FASTA streams without SeqIO/Counter, then
    # reconstruct every produced character count with per-character string counts.
    actual = {(r['matrix'], r['taxon_id']): r for r in csv.DictReader((output / 'taxon_character_coverage.tsv').open(), delimiter='\t')}
    assert len(actual) == len(characters) == 1052
    for matrix, (path, sites) in matrices.items():
        sequences = {}; name = None; chunks = []
        for line in Path(path).read_text().splitlines():
            if line.startswith('>'):
                if name is not None: assert name not in sequences; sequences[name] = ''.join(chunks)
                name = line[1:].split()[0]; chunks = []
            else: chunks.append(line.strip())
        assert name not in sequences; sequences[name] = ''.join(chunks); assert set(sequences) == set(taxa)
        for tid, seq in sequences.items():
            row = actual[matrix, tid]; total = sum(seq.count(a) for a in AA); gaps = seq.count('-'); other = {a: seq.count(a) for a in sorted(set(seq) - set(AA + '-'))}
            assert len(seq) == sites and row['sites'] == str(sites)
            assert row['canonical_aa'] == str(total) and row['gaps'] == str(gaps) and row['noncanonical_or_ambiguous'] == str(sites - total - gaps)
            assert float(row['canonical_fraction']) == total / sites and float(row['gap_fraction']) == gaps / sites and float(row['noncanonical_fraction']) == (sites - total - gaps) / sites
            assert json.loads(row['other_characters_json']) == other
            assert all(row[field] == taxa[tid][field] for field in ['species_name', 'study_role', 'lineage'])
    # Independent raw DendroPy topology parser exhausts all8views/523edges,
    # retaining every equally-near boundary rather than selecting a favored one.
    observed = {(r['view'], tuple(json.loads(r['nearest_boundary_taxa_json']))): r for r in csv.DictReader((output / 'nearest_role_boundary.tsv').open(), delimiter='\t')}
    assert len(observed) == len(role_rows); expected = set(); ns = dendropy.TaxonNamespace()
    for view in views:
        label, kind = view.split(':'); spec = source_plan['runs'][label]; file = 'pmsf.treefile' if kind == 'ml' else 'pmsf.contree'
        tree = dendropy.Tree.get(path=str(Path(spec['run']) / file), schema='newick', rooting='force-unrooted', preserve_underscores=True, taxon_namespace=ns)
        assert {n.taxon.label for n in tree.leaf_node_iter()} == set(taxa); edges = {}
        for node in tree.preorder_node_iter():
            if node is tree.seed_node or node.is_leaf(): continue
            side = {n.taxon.label for n in node.leaf_iter()}; key = min(tuple(sorted(side)), tuple(sorted(set(taxa) - side)), key=lambda v: (len(v), v))
            assert key not in edges; edges[key] = node
        assert len(edges) == 523; minimum = min(len(set(side) ^ outgroups) for side in edges)
        audit = {}
        for row in sources[label]['rows']:
            if row['tree'] != kind: continue
            side = set(json.loads(row['split_taxa_json'])); key = min(tuple(sorted(side)), tuple(sorted(set(taxa) - side)), key=lambda v: (len(v), v)); audit[key] = row
        for side, node in edges.items():
            if len(set(side) ^ outgroups) != minimum: continue
            state = view, side; expected.add(state); row = observed[state]; missing = sorted(outgroups.difference(side)); added = sorted(set(side).difference(outgroups))
            assert row['minimum_symmetric_difference_to_declared_outgroups'] == str(minimum)
            assert json.loads(row['missing_outgroups_json']) == missing and json.loads(row['added_ingroup_json']) == added
            assert json.loads(row['missing_species_json']) == [taxa[k]['species_name'] for k in missing] and json.loads(row['added_species_json']) == [taxa[k]['species_name'] for k in added]
            assert float(row['branch_length']) == node.edge.length
            assert row['sh_alrt'] == audit[side]['sh_alrt_percent'] and row['empirical_ufb'] == audit[side]['empirical_ufboot_percent']
    assert set(observed) == expected; verify(bindings)
    result = dict(status='passed_full_pmsf_role_boundary_and_character_diagnostics', plan_sha256=sha(args.plan), taxon_matrix_rows=1052, tree_views=8, boundary_rows=len(role_rows),
                  matrix_sites={k: sites for k, (_, sites) in matrices.items()},
                  character_totals={k: {field: sum(row[field] for row in characters if row['matrix'] == k) for field in ['canonical_aa', 'gaps', 'noncanonical_or_ambiguous']} for k in matrices},
                  greater_than_half_noncanonical_or_gap_taxa={k: sum(row['canonical_fraction'] < .5 for row in characters if row['matrix'] == k) for k in matrices},
                  source_hashes=bindings, artifacts={name: sha(output / name) for name in ['taxon_character_coverage.tsv', 'nearest_role_boundary.tsv']},
                  scientific_eligibility=False, scope=plan['scope'])
    with (output / 'receipt.json').open('x') as f: f.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True)


if __name__ == '__main__': main()
