#!/usr/bin/env python3
"""Check native restart fixtures against their known orthology and duplication truth."""
import argparse
import csv
import gzip
import hashlib
from io import StringIO
import json
from pathlib import Path
from Bio import Phylo


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def table(p):
    compressed = p.with_suffix(p.suffix + '.gz')
    if p.exists() and compressed.exists():
        raise ValueError('Ambiguous plain and compressed table: ' + str(p))
    with (p.open() if p.exists() else gzip.open(compressed, 'rt')) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--fixture', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    r = json.loads((a.fixture / 'receipt.json').read_text())
    cases = {x['case']: x for x in r['cases']}
    assert set(cases) == {'missing_ids_tree', 'prepared_ids_tree'}
    for name, c in cases.items():
        assert c['returncode'] == 0 and not c['timed_out']
        assert not c['source_files_changed'] and not c['source_files_added']
        assert c['source_before'] == c['source_after']
        source = a.fixture / name / 'Source'
        now = {str(p.relative_to(source)): sha(p) for p in source.rglob('*') if p.is_file()}
        assert now == c['source_before']
        assert sha(a.fixture / name / 'stdout.log') == c['log_sha256']
    missing = a.fixture / 'missing_ids_tree/Results_native_fixture'
    assert not (missing / 'Phylogenetic_Hierarchical_Orthogroups/N0.tsv').exists()
    assert not (missing / 'Comparative_Genomics_Statistics/Statistics_Overall.tsv').exists()
    folder = a.fixture / 'prepared_ids_tree/Results_native_fixture'
    expected = {f'OG{g:07d}': {f'protein_{i}_{g}' for i in range(4)} | ({'protein_0_3'} if g == 2 else set()) for g in range(3)}
    hogs = table(folder / 'Phylogenetic_Hierarchical_Orthogroups/N0.tsv')
    assert len(hogs) == 3 and {x['OG'] for x in hogs} == set(expected)
    for row in hogs:
        genes = [g.strip() for i in range(4) for g in row[f'Taxon{i}'].split(',')]
        assert len(genes) == len(set(genes)) and set(genes) == expected[row['OG']]
        for i in range(4):
            assert all(g.strip().startswith(f'protein_{i}_') for g in row[f'Taxon{i}'].split(','))
    duplicates = table(folder / 'Gene_Duplication_Events/Duplications.tsv')
    assert len(duplicates) == 1
    d = duplicates[0]
    assert d['Orthogroup'] == 'OG0000002' and d['Species Tree Node'] == 'Taxon0'
    assert d['Type'] == 'Terminal' and float(d['Support']) == 1
    assert {d['Genes 1'], d['Genes 2']} == {'Taxon0_protein_0_2', 'Taxon0_protein_0_3'}
    actual_pairs = []
    for i in range(4):
        for row in table(folder / f'Orthologues/Taxon{i}.tsv'):
            j = int(row['Species'].removeprefix('Taxon'))
            assert i != j
            actual_pairs.extend((row['Orthogroup'], left.strip(), right.strip())
                                for left in row[f'Taxon{i}'].split(',') for right in row['Orthologs'].split(','))
    expected_pairs = {(og, left, right) for og, genes in expected.items()
                      for left in genes for right in genes if left.split('_')[1] != right.split('_')[1]}
    assert len(actual_pairs) == len(set(actual_pairs)) == 42 and set(actual_pairs) == expected_pairs
    seen = set()
    for line in (folder / 'Resolved_Gene_Trees/Resolved_Gene_Trees.txt').read_text().splitlines():
        og, text = line.split(': ', 1)
        assert og not in seen; seen.add(og)
        tips = [x.name for x in Phylo.read(StringIO(text), 'newick').get_terminals()]
        target = {f'Taxon{g.split("_")[1]}_{g}' for g in expected[og]}
        assert len(tips) == len(set(tips)) and set(tips) == target
    assert seen == set(expected)
    species = Phylo.read(folder / 'Species_Tree/SpeciesTree_rooted.txt', 'newick')
    assert {x.name for x in species.get_terminals()} == {f'Taxon{i}' for i in range(4)}
    assert {frozenset(x.name for x in clade.get_terminals()) for clade in species.root.clades} == {frozenset(['Taxon0', 'Taxon1']), frozenset(['Taxon2', 'Taxon3'])}
    with (folder / 'Comparative_Genomics_Statistics/Statistics_Overall.tsv').open() as handle:
        stats = {x[0]: x[1] for x in csv.reader(handle, delimiter='\t') if len(x) >= 2}
    assert all(int(stats[k]) == v for k, v in {'Number of species': 4, 'Number of genes': 13, 'Number of orthogroups': 3, 'Number of unassigned genes': 0}.items())
    root = Path(__file__).resolve().parents[1]
    package = root / '.cache/envs/orthofinder/lib/python3.12/site-packages/orthofinder'
    result = dict(status='passed_native_restart_fixture_output_readback',
                  genes_checked=13, families_checked=3, directed_ortholog_pairs_checked=42,
                  terminal_duplications_checked=1, source_files_unchanged=True,
                  missing_ids_case_exit_zero_despite_incomplete_outputs=True,
                  fixture_receipt_sha256=sha(a.fixture / 'receipt.json'), script_sha256=sha(Path(__file__)),
                  installed_sources={str(p.relative_to(root)): sha(p) for p in [package / 'run/main.py', package / 'run/process_args.py', package / 'utils/files.py', package / 'comparative_genomics/orthologues.py']},
                  checked_output_hashes={str(p.relative_to(folder)): sha(p) for p in folder.rglob('*') if p.is_file()},
                  interpretation='Actual installed --from-trees --no-fix-files completes this prepared-ID fixture and recovers known outputs without source changes. Missing IDs fixture returns zero but has no final HOG/statistics outputs. This is software contract evidence, not scientific validation or proof of source protection under other flags/datasets.')
    a.output.write_text(json.dumps(result, indent=2) + '\n')
    print('Verified 13 genes, 3 families, 42 directed ortholog pairs, 1 duplication; false-zero failure detected')


if __name__ == '__main__':
    main()
