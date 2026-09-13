#!/usr/bin/env python3
"""Combine reviewed sampling drafts and expose unresolved annotation availability."""
import csv,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
rows=[]
for name in ('fungal_sampling_expanded_draft.tsv','outgroup_sampling_draft.tsv'):
 rows.extend(csv.DictReader((ROOT/'metadata'/name).open(),delimiter='\t'))
seen={}
for line in (ROOT/'data/raw/proteome_downloads.jsonl').read_text().splitlines():
 r=json.loads(line);seen[r['url']]=r['status']
for r in json.loads((ROOT/'metadata/external_genome_receipts.json').read_text()):
 if 'proteins' in r:seen[r['url']]='validated' if r.get('fasta_status')=='validated' else 'requires_normalization'
for r in rows:
 s=seen.get(r['proteome_url'],'pending')
 r['status']={'validated':'provisional_file_QC_passed','requires_normalization':'provisional_normalization_required','error':'annotation_unavailable'}.get(s,'provisional_pending_acquisition')
assert len(rows)==527
assert sum(r['study_role']=='outgroup' for r in rows)==25
assert len({r['taxon_id'] for r in rows})==527
assert len({r['species_name'] for r in rows})==527
with (ROOT/'metadata/sampling_manifest.tsv').open('w') as f:
 w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
print('Manifest: 502 fungal candidates + 25 outgroups; status records QC and unresolved annotation')
