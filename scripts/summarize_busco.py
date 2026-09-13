#!/usr/bin/env python3
"""Summarize only successful BUSCO jobs with internally consistent marker counts."""
import csv,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
rows=[]
for receipt in sorted((ROOT/'results/busco').glob('*.receipt.json')):
 r=json.loads(receipt.read_text())
 if r.get('returncode')!=0:continue
 summaries=list((ROOT/'results/busco'/r['taxon_id']).glob('short_summary.specific.*.json'))
 if len(summaries)!=1:raise ValueError('Expected one summary for '+r['taxon_id'])
 d=json.loads(summaries[0].read_text());s=d['results']
 assert s['Complete BUSCOs']==s['Single copy BUSCOs']+s['Multi copy BUSCOs']
 assert s['n_markers']==s['Complete BUSCOs']+s['Fragmented BUSCOs']+s['Missing BUSCOs']
 rows.append({'taxon_id':r['taxon_id'],'input_sha256':r['input_sha256'],'complete_percent':s['Complete percentage'],'single_copy_percent':s['Single copy percentage'],'duplicated_percent':s['Multi copy percentage'],'fragmented_percent':s['Fragmented percentage'],'missing_percent':s['Missing percentage'],'markers':s['n_markers'],'summary_sha256':hashlib.sha256(summaries[0].read_bytes()).hexdigest(),'summary_path':str(summaries[0].relative_to(ROOT))})
if rows:
 with (ROOT/'metadata/busco_eukaryota_qc.tsv').open('w') as out:
  w=csv.DictWriter(out,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
print(len(rows),'successful jobs summarized; raw-proteome QC only')
