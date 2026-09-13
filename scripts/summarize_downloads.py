#!/usr/bin/env python3
"""Snapshot completed QC records without treating running downloads as failures."""
import csv,json
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader((ROOT/'metadata/fungal_sampling_draft.tsv').open(),delimiter='\t'))
by_url={}
for line in (ROOT/'data/raw/proteome_downloads.jsonl').read_text().splitlines():
 try:r=json.loads(line)
 except json.JSONDecodeError:continue # writer may be appending the final record
 by_url[r['url']]=r
receipts=[by_url[r['proteome_url']] for r in rows if r['proteome_url'] in by_url]
valid=[r for r in receipts if r['status']=='validated']
summary={'target':len(rows),'recorded':len(receipts),'pending_without_terminal_record':len(rows)-len(receipts),'status_counts':dict(Counter(r['status'] for r in receipts)),'validated_proteins':sum(r['proteins'] for r in valid),'validated_residues':sum(r['residues'] for r in valid),'compressed_bytes':sum(r['compressed_bytes'] for r in valid),'note':'Snapshot only; per-process liveness must be checked through execution handle. Biological completeness and isoform filtering not assessed.'}
(ROOT/'metadata/proteome_qc_snapshot.json').write_text(json.dumps(summary,indent=2)+'\n')
(ROOT/'metadata/proteome_download_receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')
print(json.dumps(summary,indent=2))
