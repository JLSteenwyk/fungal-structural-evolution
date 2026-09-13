#!/usr/bin/env python3
"""Join assembly records to an NCBI ranked-lineage snapshot; audit coverage."""
import csv
import hashlib
import json
import tarfile
from collections import defaultdict
from pathlib import Path
from datetime import datetime, timezone
ROOT = Path(__file__).resolve().parents[1]


def fields(line):
    return [p.strip() for p in line.rstrip('\n').split('|')][:-1]


def main():
    archive = ROOT / 'data/raw/new_taxdump.tar.gz'
    records = list(csv.DictReader((ROOT/'metadata/assembly_candidates.tsv').open(), delimiter='\t'))
    wanted = {r['species_taxid'] for r in records}
    taxonomy = {}
    with tarfile.open(archive) as tar:
        with tar.extractfile('rankedlineage.dmp') as stream:
            for line in stream:
                parts = fields(line.decode())
                if parts[0] in wanted:
                    taxonomy[parts[0]] = dict(zip(['species_taxid','taxon_name','species','genus','family','order','class','phylum','kingdom','superkingdom'], parts))
    extra = ['taxon_name','genus','family','order','class','phylum','kingdom','superkingdom']
    missing = []
    coverage = defaultdict(lambda: {'all':set(),'annotated':set()})
    output = ROOT/'metadata/assembly_candidates_taxonomy.tsv'
    with output.open('w') as handle:
        writer=csv.DictWriter(handle,list(records[0])+extra,delimiter='\t',lineterminator='\n')
        writer.writeheader()
        for row in records:
            tax = taxonomy.get(row['species_taxid'])
            if tax is None: missing.append(row['species_taxid'])
            row.update({k: (tax or {}).get(k,'') for k in extra})
            writer.writerow(row)
            group=row['phylum'] or 'unresolved_phylum'
            coverage[group]['all'].add(row['species_taxid'])
            if row['annotation_name'] not in ('','na'): coverage[group]['annotated'].add(row['species_taxid'])
    summary={k:{'species':len(v['all']),'species_with_annotation_name':len(v['annotated'])} for k,v in sorted(coverage.items())}
    report={'coverage_by_ncbi_phylum':summary,'missing_species_taxids':sorted(set(missing)), 'caution':'NCBI ranks are source taxonomy, not the project ingroup/outgroup definition. Annotation name does not verify proteome availability.'}
    (ROOT/'metadata/taxonomy_coverage.json').write_text(json.dumps(report,indent=2)+'\n')
    receipts=ROOT/'metadata/source_receipts.json'
    d=json.loads(receipts.read_text())
    d['ncbi_taxonomy']={'url':'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz','path':str(archive.relative_to(ROOT)), 'sha256':hashlib.file_digest(archive.open('rb'),'sha256').hexdigest() if hasattr(hashlib,'file_digest') else hashlib.sha256(archive.read_bytes()).hexdigest(),'bytes':archive.stat().st_size,'file_mtime_utc':datetime.fromtimestamp(archive.stat().st_mtime,timezone.utc).isoformat(),'retrieved_at_utc':None,'provenance_note':'Initial curl retrieval; file mtime recorded. Preserve this snapshot; future downloads may change.'}
    receipts.write_text(json.dumps(d,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
