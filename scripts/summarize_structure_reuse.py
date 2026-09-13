#!/usr/bin/env python3
"""Snapshot structure retrieval without interpreting queued candidates as coverage."""
import json,csv
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
latest={}
for line in (ROOT/'data/raw/afdb_models.jsonl').read_text().splitlines():
 try:r=json.loads(line)
 except json.JSONDecodeError:continue
 latest[r['uniprot_accession']]=r
models=[dict(uniprot_accession=r['uniprot_accession'],**m) for r in latest.values() if r['status']=='verified' for m in r['models']]
matched={(m['uniprot_accession'],m['sequence_sha256']) for m in models}
proteins=set()
for row in csv.DictReader((ROOT/'data/structures/afdb/input_links.tsv').open(),delimiter='\t'):
 if (row['uniprot_accession'],row['sequence_sha256']) in matched:proteins.add((row['taxon_id'],row['protein_id']))
summary={'accession_status_counts':dict(Counter(r['status'] for r in latest.values())),'verified_models':len(models),'distinct_input_proteins_linked':len(proteins),'distinct_input_taxa_linked':len({x[0] for x in proteins}),'pae_downloaded_models':sum(m['pae_downloaded'] for m in models),'note':'Partial retrieval snapshot; full-length exact sequence and CA coverage required. Confidence filtering, domain parsing and PAE assessment remain pending.'}
(ROOT/'metadata/structure_reuse_snapshot.json').write_text(json.dumps(summary,indent=2)+'\n')
# Full coordinate and API receipts remain outside Git; store manageable per-model audit snapshot.
(ROOT/'metadata/verified_structure_receipts.json').write_text(json.dumps(models,indent=2)+'\n')
print(json.dumps(summary,indent=2))
