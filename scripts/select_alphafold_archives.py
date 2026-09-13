#!/usr/bin/env python3
"""Deduplicate AFDB bulk shard versions without claiming sequence-level coverage."""
import json,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
d=json.loads((ROOT/'metadata/alphafold_bulk_objects.json').read_text());best={}
for obj in d['objects']:
 m=re.fullmatch(r'(proteomes/proteome-tax_id-\d+-\d+)_v(\d+)\.tar',obj['name'])
 if not m:raise ValueError('Unexpected bulk archive name: '+obj['name'])
 key,version=m[1],int(m[2])
 if key not in best or version>best[key][0]:best[key]=(version,obj)
selected=[o for _,o in sorted(best.values(),key=lambda pair:pair[1]['name'])]
result={'bucket':d['bucket'],'selection_rule':'Highest available suffix version for each exact taxid/shard; bucket name is not a version assertion','selected_objects':selected,'archive_bytes':sum(int(o['size']) for o in selected),'warning':'Taxon archives include unverified annotation/strain matches. Exact sequence matching required before reuse; no archive at queried IDs does not establish absence from current AFDB.'}
(ROOT/'metadata/alphafold_selected_archives.json').write_text(json.dumps(result,indent=2)+'\n')
print(len(selected),'archives;',round(result['archive_bytes']/1e9,2),'GB')
