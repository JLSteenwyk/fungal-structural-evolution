#!/usr/bin/env python3
"""Partition all recovered marker links by availability in two frozen catalogs."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def rows(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write(path, data):
    with path.open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(data[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(data)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    a = p.parse_args()
    plan = json.loads(a.plan.read_text())
    pins = {str(a.plan): sha(a.plan), **plan['pins']}
    def verify():
        for path, digest in pins.items():
            if sha(path) != digest:
                raise ValueError('Changed input: ' + path)
    verify()
    out = Path(plan['output'])
    if out.exists():
        raise FileExistsError(out)
    links = rows(plan['global_links'])
    manifest = rows(plan['manifest'])
    keys = lambda r: (r['marker'], r['taxon_id'], r['protein_id'], r['sequence_sha256'])
    universe = {keys(r) for r in links}
    if len(universe) != len(links):
        raise ValueError('Duplicate global link')
    known = {r['taxon_id'] for r in manifest}
    if len(known) != len(manifest) or any(r['taxon_id'] not in known for r in links):
        raise ValueError('Invalid taxon universe')
    sources = {}
    for name in ['esmfold', 'alphafold']:
        rr = rows(plan[name + '_links'])
        sources[name] = {keys(r) for r in rr}
        if len(sources[name]) != len(rr) or not sources[name] <= universe:
            raise ValueError('Invalid source links: ' + name)
    er = json.loads(Path(plan['esmfold_readback']).read_text())
    ar = json.loads(Path(plan['alphafold_receipt']).read_text())
    if (er['status'] != 'passed_complete_disjoint_prediction_inventory_readback'
            or ar['status'] != 'complete_source_specific_model_catalog'
            or len(sources['esmfold']) != er['matched_marker_protein_links']
            or len(sources['alphafold']) != ar['marker_proteins_linked']
            or ar['marker_mapping_sha256'] != sha(plan['global_links'])
            or ar['artifacts']['marker_model_links.tsv'] != sha(plan['alphafold_links'])
            or er['artifacts']['expected_global_links.tsv'] != sha(plan['esmfold_links'])):
        raise ValueError('Source receipts do not bind the supplied catalogs')
    counts = Counter()
    per_taxon = {t: Counter() for t in known}
    states = ['both', 'esmfold_only', 'alphafold_only', 'neither_catalog']
    labeled = []
    for row in links:
        key = keys(row)
        e, f = key in sources['esmfold'], key in sources['alphafold']
        state = 'both' if e and f else 'esmfold_only' if e else 'alphafold_only' if f else 'neither_catalog'
        counts[state] += 1
        per_taxon[row['taxon_id']][state] += 1
        labeled.append(dict(row, availability=state))
    nmarkers = len({r['marker'] for r in links})
    taxa = []
    for r in manifest:
        c = per_taxon[r['taxon_id']]
        recovered = sum(c.values())
        if recovered > nmarkers:
            raise ValueError('Multiple recovered proteins per marker/taxon')
        taxa.append(dict(taxon_id=r['taxon_id'], species_name=r['species_name'],
                         study_role=r['study_role'], marker_slots=nmarkers,
                         recovered_sequence_links=recovered,
                         missing_marker_sequences=nmarkers-recovered,
                         **{s: c[s] for s in states}))
    verify()
    out.mkdir(parents=True)
    write(out/'marker_availability.tsv', labeled)
    write(out/'taxon_availability.tsv', taxa)
    result = dict(status='complete_frozen_catalog_marker_availability',
                  taxa=len(manifest), markers=nmarkers, marker_links=len(links),
                  distinct_sequences=len({r['sequence_sha256'] for r in links}),
                  link_counts={s: counts[s] for s in states},
                  links_with_model=len(sources['esmfold'] | sources['alphafold']),
                  taxa_with_model=sum(sum(per_taxon[t][s] for s in states[:3]) > 0 for t in known),
                  plan_sha256=sha(a.plan), script_sha256=sha(__file__),
                  artifacts={p.name: sha(p) for p in out.iterdir()},
                  scope='Exact marker/protein links in the two specified frozen catalogs. Neither catalog is not proof of absence from current caches or public databases. Counts precede confidence qualification and are not full-proteome coverage. Identical sequences and shared models are not independent observations.')
    (out/'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
