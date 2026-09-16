#!/usr/bin/env python3
"""Account for every native protein across the retained family partition."""
import argparse
import csv
import hashlib
import json
import time
from collections import Counter
from pathlib import Path
from audit_orthology_family_universe import groups


def sha(p):
    with p.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    for name, expected in plan['pins'].items():
        if sha(Path(name)) != expected:
            raise ValueError(f'Pin changed: {name}')
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    wd = Path(plan['working_directory'])
    originals = Path(plan['original_working_directory'])
    with Path(plan['staging_manifest']).open() as handle:
        staged = {Path(r['relative_path']).name: r['sha256'] for r in csv.DictReader(handle, delimiter='\t')}
    species = {int(line.split(':',1)[0]): line.split(':',1)[1].strip().removesuffix('.faa')
               for line in (wd / 'SpeciesIDs.txt').read_text().splitlines() if line.strip()}
    flags = {sid: bytearray() for sid in species}
    start = time.time()
    with (wd / 'SequenceIDs.txt').open() as handle:
        for line in handle:
            native, sep, label = line.rstrip('\n').partition(': ')
            sid, index = map(int, native.split('_'))
            if not sep or not label or sid not in flags or index != len(flags[sid]):
                raise ValueError('Unexpected or duplicate sequence mapping')
            flags[sid].append(0)
    families = Counter()
    for number, genes in groups(wd / 'clusters_OrthoFinder.txt_id_pairs.txt'):
        kind = min(len(genes), 3)
        if not kind:
            raise ValueError('Empty family')
        families[kind] += 1
        for native in genes:
            sid, index = map(int, native.split('_'))
            if sid not in flags or index < 0 or index >= len(flags[sid]) or flags[sid][index]:
                raise ValueError('Unknown or multiply assigned family protein')
            flags[sid][index] = kind
    print('Mapped all sequence IDs and retained family memberships', flush=True)
    rows, provenance = [], {}
    for sid in sorted(species):
        source = wd / f'Species{sid}.fa'
        n = 0
        with source.open() as handle:
            for line in handle:
                if line.startswith('>'):
                    if line[1:].strip() != f'{sid}_{n}':
                        raise ValueError('FASTA/mapping identity or order mismatch')
                    n += 1
        if n != len(flags[sid]) or sha(source) != staged[source.name]:
            raise ValueError('FASTA identity coverage or checksum mismatch')
        legacy = originals / f'Unassigned.Species{sid}.fa'
        seen = bytearray(n)
        legacy_counts = Counter()
        before = sha(legacy)
        with legacy.open() as handle:
            for line in handle:
                if line.startswith('>'):
                    s, i = map(int, line[1:].strip().split('_'))
                    if s != sid or i < 0 or i >= n or seen[i]:
                        raise ValueError('Unknown/duplicate intermediate unassigned ID')
                    seen[i] = 1
                    legacy_counts[flags[sid][i]] += 1
        if sha(legacy) != before:
            raise ValueError('Intermediate input changed')
        provenance[str(legacy)] = before
        c = Counter(flags[sid])
        rows.append(dict(species_id=sid, taxon_id=species[sid], proteins=n,
                         outside_retained_partition=c[0], singleton_family_proteins=c[1],
                         two_member_family_proteins=c[2], tree_eligible_family_proteins=c[3],
                         outside_fraction=c[0]/n,
                         intermediate_unassigned_proteins=sum(legacy_counts.values()),
                         intermediate_unassigned_now_in_partition=sum(legacy_counts[k] for k in [1,2,3]),
                         intermediate_unassigned_still_outside=legacy_counts[0],
                         outside_not_in_intermediate=c[0]-legacy_counts[0]))
    for name in ['SpeciesIDs.txt','SequenceIDs.txt','clusters_OrthoFinder.txt_id_pairs.txt']:
        if sha(wd / name) != staged[name]:
            raise ValueError('Staged mapping/partition changed')
    with (out / 'taxon_coverage.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(rows)
    omitted = 0
    with (out / 'outside_retained_partition.tsv').open('w') as output, (wd / 'SequenceIDs.txt').open() as handle:
        writer = csv.writer(output, delimiter='\t')
        writer.writerow(['native_id','taxon_id','protein_id'])
        for line in handle:
            native, label = line.rstrip('\n').split(': ',1)
            sid, i = map(int,native.split('_'))
            if flags[sid][i] == 0:
                writer.writerow([native,species[sid],label]); omitted += 1
    totals = {key:sum(r[key] for r in rows) for key in rows[0] if key not in ['species_id','taxon_id','outside_fraction']}
    if omitted != totals['outside_retained_partition'] or totals['proteins'] != sum(totals[k] for k in ['outside_retained_partition','singleton_family_proteins','two_member_family_proteins','tree_eligible_family_proteins']):
        raise ValueError('Protein accounting does not close')
    (out / 'intermediate_source_checksums.json').write_text(json.dumps(provenance,indent=2)+'\n')
    result = dict(status='complete_full_protein_family_coverage_accounting',taxa=len(rows),totals=totals,
                  family_counts_by_size_class=dict(families),outside_fraction=omitted/totals['proteins'],
                  elapsed_seconds=time.time()-start,plan_sha256=sha(args.plan),script_sha256=sha(Path(__file__)),
                  artifacts={p.name:sha(p) for p in out.iterdir() if p.is_file()},
                  interpretation='Exact native-ID accounting against the saved family partition and full species FASTAs. Historical Unassigned.Species files are compared as intermediate inputs, not presumed final assignments. Outside-partition proteins are retained for atlas/coverage work; this status is not evidence of biological novelty, missing homologs or gene loss. No family memberships changed.')
    (out / 'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__ == '__main__':
    main()
