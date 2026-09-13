#!/usr/bin/env python3
"""Record the available full-scale analysis set, retaining excluded candidate evidence."""
import csv,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
candidates=list(csv.DictReader((ROOT/'metadata/sampling_manifest.tsv').open(),delimiter='\t'))
inputs={r['taxon_id']:r for r in json.loads((ROOT/'metadata/qc_input_receipts.json').read_text())}
selected=[];excluded=[]
for r in candidates:
 if r['taxon_id'] in inputs:selected.append(r)
 else:excluded.append({'taxon_id':r['taxon_id'],'species_name':r['species_name'],'assembly_accession':r['assembly_accession'],'reason':'No usable annotated proteome acquired; retained in candidate manifest for future recovery','evidence':'metadata/proteome_download_receipts.json; docs/progress.md'})
assert len(selected)==526
assert sum(r['study_role']=='ingroup' for r in selected)==501
assert sum(r['study_role']=='outgroup' for r in selected)==25
with (ROOT/'metadata/analysis_manifest.tsv').open('w') as out:
 w=csv.DictWriter(out,list(selected[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(selected)
(ROOT/'metadata/sampling_exclusions.json').write_text(json.dumps(excluded,indent=2)+'\n')
print('Working analysis set: 501 fungi + 25 outgroups; 1 source-limited candidate explicitly excluded')
