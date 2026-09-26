#!/usr/bin/env python3
"""Stream grouped ortholog tables and validate every protein's source family.

Counts are directed incidences, not deduplicated pairs or validated orthology.
"""
import argparse
import csv
import hashlib
import json
import shutil
from collections import Counter
from pathlib import Path

from assess_small_family_output_exposure import groups, sha


def labels(source):
    species = {}
    for line in (source / 'SpeciesIDs.txt').read_text().splitlines():
        native, label = line.split(': ', 1)
        if native in species:
            raise ValueError('Duplicate species ID')
        species[native] = label.rsplit('.', 1)[0]
    if len(set(species.values())) != len(species):
        raise ValueError('Ambiguous species labels')
    family_by_gene = {}
    for family, genes in groups(source / 'clusters_OrthoFinder.txt_id_pairs.txt'):
        for gene in genes:
            if gene in family_by_gene:
                raise ValueError('Repeated gene in partition')
            family_by_gene[gene] = family
    lookup = {taxon: {} for taxon in species.values()}
    with (source / 'SequenceIDs.txt').open() as handle:
        for line in handle:
            native, protein = line.rstrip('\n').split(': ', 1)
            if len(protein.split()) != 1 or any(c in protein for c in ',:()'):
                raise ValueError('Unsupported protein naming mode')
            taxon = species[native.split('_')[0]]
            if protein in lookup[taxon] or native not in family_by_gene:
                raise ValueError('Ambiguous or unassigned protein')
            lookup[taxon][protein] = family_by_gene.pop(native)
    if family_by_gene:
        raise ValueError('Missing source protein labels')
    return lookup


def inspect_row(raw, taxon, lookup):
    row = next(csv.reader([raw.decode()], delimiter='\t'))
    if len(row) != 4:
        raise ValueError('Malformed grouped ortholog row')
    family, target, left, right = row
    if target == taxon or target not in lookup:
        raise ValueError('Invalid target species')
    sides = [left.split(', '), right.split(', ')]
    for species, proteins in zip([taxon, target], sides):
        if len(proteins) != len(set(proteins)):
            raise ValueError('Repeated protein inside grouped row')
        if any(lookup[species].get(protein) != family for protein in proteins):
            raise ValueError('Unknown protein or wrong family/species')
    return family, len(sides[0]) * len(sides[1])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    for name, digest in plan['pins'].items():
        if sha(name) != digest:
            raise ValueError('Changed plan input: ' + name)
    output = Path(plan['output'])
    if output.exists():
        raise FileExistsError(output)
    if shutil.disk_usage(output.parent).free < plan['resources']['minimum_free_disk_gib'] * 2**30:
        raise ValueError('Insufficient disk')
    snapshot = json.loads(Path(plan['supplement_receipt']).read_text())
    if snapshot['status'] != 'complete_separate_small_family_ortholog_supplement':
        raise ValueError('Completed source snapshot required')
    for path, digest in snapshot['input_hashes'].items():
        if sha(path) != digest:
            raise ValueError('Changed source identity input')
    roots = {Path(p).parent for p in snapshot['input_hashes']}
    if len(roots) != 1:
        raise ValueError('Ambiguous source directory')
    lookup = labels(roots.pop())
    paths = {Path(p).stem: Path(p) for p in snapshot['native_ortholog_hashes']}
    if paths.keys() != lookup.keys() or len(paths) != len(snapshot['native_ortholog_hashes']):
        raise ValueError('Incomplete native species tables')
    output.mkdir()
    totals = Counter()
    incidences = Counter()
    files = []
    for taxon, path in sorted(paths.items()):
        digest = hashlib.sha256()
        before = path.stat()
        nrows = npairs = 0
        with path.open('rb') as handle:
            header = handle.readline()
            digest.update(header)
            if header.decode().rstrip('\r\n').split('\t') != ['Orthogroup', 'Species', taxon, 'Orthologs']:
                raise ValueError('Unexpected header')
            for raw in handle:
                digest.update(raw)
                family, count = inspect_row(raw, taxon, lookup)
                totals[family] += 1
                incidences[family] += count
                nrows += 1
                npairs += count
        after = path.stat()
        if (before.st_size, before.st_mtime_ns, before.st_ino) != (after.st_size, after.st_mtime_ns, after.st_ino):
            raise ValueError('Native table changed while reading')
        if digest.hexdigest() != snapshot['native_ortholog_hashes'][str(path)]:
            raise ValueError('Native table differs from completed snapshot')
        files.append(dict(taxon=taxon, rows=nrows, directed_incidences=npairs, sha256=digest.hexdigest()))
        (output / 'state.json').write_text(json.dumps(dict(completed_tables=len(files), total_tables=len(paths), rows=sum(r['rows'] for r in files))) + '\n')
        print(json.dumps(files[-1]), flush=True)
    if sum(totals.values()) != snapshot['native_rows_scanned']:
        raise ValueError('Native row count differs')
    for path, digest in snapshot['input_hashes'].items():
        if sha(path) != digest:
            raise ValueError('Source identity input changed during scan')
    for filename, rows, fields in [
        ('family_counts.tsv', [dict(family=k, rows=totals[k], directed_incidences=incidences[k]) for k in sorted(totals)], ['family', 'rows', 'directed_incidences']),
        ('table_counts.tsv', files, ['taxon', 'rows', 'directed_incidences', 'sha256'])]:
        with (output / filename).open('x') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)
    result = dict(status='complete_grouped_ortholog_protein_family_identity_audit',
                  plan_sha256=sha(args.plan), script_sha256=sha(__file__),
                  taxa=len(lookup), proteins=sum(map(len, lookup.values())),
                  rows=sum(totals.values()), directed_incidences=sum(incidences.values()),
                  families_with_rows=len(totals),
                  artifacts={name: sha(output / name) for name in ['family_counts.tsv', 'table_counts.tsv']},
                  scope='Every native grouped row has valid species/protein/family identities and no within-side repeats; native file hashes match the completed supplement snapshot. Counts are directed row incidences, not unique pairs. Cross-row overlap, reciprocal equality, missing larger-family pairs, reconciliation semantics, biological orthology and corrected integrated statistics remain unvalidated.')
    (output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
