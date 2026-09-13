#!/usr/bin/env python3
"""Retrieve NCBI catalog snapshots and inventory assemblies without selecting taxa."""
import argparse
import csv
import hashlib
import io
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]


def parse_catalog(text):
    lines = text.splitlines()
    header = next((i for i, line in enumerate(lines) if line.lstrip('# ').startswith('assembly_accession\t')), None)
    if header is None:
        raise ValueError('Missing NCBI assembly catalog header')
    lines[header] = lines[header].lstrip('# ')
    rows = list(csv.DictReader(io.StringIO('\n'.join(lines[header:])), delimiter='\t'))
    if not rows or any(not r.get('assembly_accession') or None in r for r in rows):
        raise ValueError('Empty or malformed catalog')
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--refresh', action='store_true', help='Explicitly replace cached snapshots')
    args = parser.parse_args()
    config = json.loads((ROOT / 'config/project.json').read_text())
    receipt_path = ROOT / 'metadata/source_receipts.json'
    receipts = json.loads(receipt_path.read_text()) if receipt_path.exists() else {}
    records = []
    for name, url in config['sources'].items():
        path = ROOT / 'data/raw' / (name + '_assembly_summary.txt')
        path.parent.mkdir(parents=True, exist_ok=True)
        if args.refresh or not path.exists():
            with urlopen(url, timeout=120) as response:
                content = response.read()
            parse_catalog(content.decode())
            temporary = path.with_suffix('.tmp')
            temporary.write_bytes(content)
            temporary.replace(path)
            receipts[name] = {'retrieved_at_utc': datetime.now(timezone.utc).isoformat()}
        content = path.read_bytes()
        digest = hashlib.sha256(content).hexdigest()
        old = receipts.get(name, {})
        if old.get('sha256') and old['sha256'] != digest and not args.refresh:
            raise ValueError(f'Cached snapshot changed unexpectedly: {name}')
        receipts[name] = dict(old, url=url, path=str(path.relative_to(ROOT)), sha256=digest, bytes=len(content))
        if 'retrieved_at_utc' not in receipts[name]:
            receipts[name]['retrieved_at_utc'] = None
            receipts[name]['provenance_note'] = 'Initial curl retrieval during project initialization; exact retrieval timestamp not logged. File mtime recorded separately.'
            receipts[name]['file_mtime_utc'] = datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat()
        rows = parse_catalog(content.decode())
        for row in rows:
            row['catalog_source'] = name
        records.extend(rows)
    fields = ['catalog_source'] + [f for f in records[0] if f != 'catalog_source']
    with (ROOT / 'metadata/assembly_candidates.tsv').open('w') as out:
        writer = csv.DictWriter(out, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(records)
    species = {r['species_taxid'] for r in records if r['species_taxid'] not in ('', 'na')}
    annotated = [r for r in records if r.get('annotation_name', '') not in ('', 'na')]
    summary = {'assembly_records':len(records), 'distinct_species_taxids':len(species),
               'records_by_source':dict(Counter(r['catalog_source'] for r in records)),
               'records_with_annotation_name':len(annotated),
               'distinct_species_with_annotation_name':len({r['species_taxid'] for r in annotated}),
               'assembly_levels':dict(Counter(r['assembly_level'] for r in records)),
               'selected_taxa':0, 'warning':'Catalog taxonomy and annotation metadata are unvalidated; paired GCA/GCF records and strains are not independent species.'}
    receipt_path.write_text(json.dumps(receipts, indent=2)+'\n')
    (ROOT / 'metadata/catalog_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    main()
