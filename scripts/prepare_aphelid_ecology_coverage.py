"""Bind full47-taxon evidence to both existing qualified alignment collections."""
import csv
import json
from pathlib import Path
from ancestral_chain_attempt import sha,write_json

source=Path('metadata/aphelid_ecology_evidence_receipt_20260928.json')
r=json.loads(source.read_text());table='metadata/species_ecology_evidence_aphelids_20260928.tsv'
assert r['status']=='complete_two_assembly_linked_aphelid_ecology_additions'
for name,h in r['artifacts'].items():assert sha(name)==h
with open(table) as f:rows=list(csv.DictReader(f,delimiter='\t'))
with open('metadata/species_ecology_evidence_reviewed_20260927.tsv') as f:old=list(csv.DictReader(f,delimiter='\t'))
assert len(rows)==47 and rows[:45]==old and len({r['taxon_id'] for r in rows})==47
binding='metadata/aphelid_ecology_coverage_evidence_binding_20260928.json'
write_json(Path(binding),dict(status='verified47_taxon_table_binding_for_coverage',table_sha256=sha(table),source_receipt_sha256=sha(source),taxa=47,preserved_rows=45,scope='Evidence-table linkage for coverage, not ecological classification validation.'))
for method in ['afdb','esmfold']:
    p=json.loads(Path(f'metadata/qualified_{method}_reviewed45_ecology_plan_20260927.json').read_text())
    p['pins'].pop(p['evidence']);p['pins'].pop(p['evidence_receipt'])
    p['evidence']=table;p['evidence_receipt']=binding
    p['output']=f'results/ecology/qualified-{method}-aphelids47-overlap-20260928-v1'
    for name in [table,binding,str(source)]:p['pins'][name]=sha(name)
    for name,h in p['pins'].items():assert sha(name)==h,name
    p['resources']['basis']='All47 taxa and1081 unordered pairs, full source-specific qualified AA/3Di masks; no prediction or cross-source pooling.'
    write_json(Path(f'metadata/qualified_{method}_aphelids47_ecology_plan_20260928.json'),p)
print('Prepared47 taxa and1081 pairs for each prediction source')
