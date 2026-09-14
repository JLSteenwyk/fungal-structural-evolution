#!/usr/bin/env python3
"""Inventory assembly alternatives for the two EXCLUDE-affected marker taxa.

Read-only review: never replace selected assemblies or interpret availability as QC.
"""
import argparse
import csv
import datetime
import hashlib
import json
from pathlib import Path
import urllib.parse
import urllib.request


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    manifest = Path('metadata/analysis_manifest.tsv')
    with manifest.open() as handle:
        targets = {r['species_taxid']: r for r in csv.DictReader(handle, delimiter='\t')
                   if r['taxon_id'] in {'F610337', 'F27376'}}
    assert len(targets) == 2
    sources = {}
    matches = {taxid: [] for taxid in targets}
    for name in ['genbank_fungi', 'refseq_fungi']:
        path = Path('data/raw') / (name + '_assembly_summary.txt')
        sources[str(path)] = sha(path)
        with path.open() as handle:
            header = next(line.lstrip('# ').rstrip().split('\t') for line in handle
                          if line.lstrip('# ').startswith('assembly_accession\t'))
            for row in csv.DictReader(handle, fieldnames=header, delimiter='\t'):
                if row['species_taxid'] in targets:
                    matches[row['species_taxid']].append({'catalog': name, **row})
    (args.output / 'frozen_catalog_matches.json').write_text(json.dumps(matches, indent=2) + '\n')
    summaries, requests = [], []
    for taxid, selected in sorted(targets.items()):
        reports, token, seen = [], '', set()
        while True:
            query = {'page_size': 100}
            if token:
                query['page_token'] = token
            url = ('https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/' + taxid
                   + '/dataset_report?' + urllib.parse.urlencode(query))
            with urllib.request.urlopen(url, timeout=60) as response:
                raw = response.read()
            path = args.output / f'{taxid}_page_{len(seen)}.json'
            path.write_bytes(raw)
            data = json.loads(raw)
            requests.append({'url': url, 'file': str(path), 'sha256': sha(path),
                             'retrieved_utc': datetime.datetime.now(datetime.timezone.utc).isoformat()})
            if 'error' in data or 'total_count' not in data:
                raise ValueError('Unrecognized API response')
            reports.extend(data.get('reports', []))
            token = data.get('next_page_token', '')
            if not token:
                break
            if token in seen:
                raise ValueError('Repeated pagination token')
            seen.add(token)
        if len(reports) != int(data['total_count']):
            raise ValueError('Incomplete API inventory')
        accessions = [r['accession'] for r in reports]
        if len(set(accessions)) != len(accessions):
            raise ValueError('Duplicate assembly records')
        summaries.append({'taxon_id': selected['taxon_id'], 'species_name': selected['species_name'],
                          'species_taxid': taxid, 'selected_assembly': selected['assembly_accession'],
                          'frozen_catalog_accessions': [r['assembly_accession'] for r in matches[taxid]],
                          'live_api_accessions': accessions,
                          'selected_found_in_live_api': selected['assembly_accession'] in accessions,
                          'other_accession_candidates': [a for a in accessions if a != selected['assembly_accession']]})
    receipt = {'status': 'complete_scoped_assembly_availability_inventory',
               'script_sha256': sha(Path(__file__)), 'manifest_sha256': sha(manifest),
               'catalog_sha256': sources, 'requests': requests, 'taxa': summaries,
               'artifacts': {p.name: sha(p) for p in args.output.iterdir()},
               'limits': 'Frozen fungal GenBank/RefSeq catalogs and current NCBI Datasets taxon endpoint only. '
                         'No claim about all external archives, unpublished data, historical assemblies or taxonomically mislabelled accessions. '
                         'Other accessions can be paired records or versions, not independent samples. '
                         'No sequence provenance resolved and no assembly replacement made.'}
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(summaries, indent=2))


if __name__ == '__main__':
    main()
