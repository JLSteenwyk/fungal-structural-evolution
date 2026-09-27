#!/usr/bin/env python3
"""Check preserved records, all curation decisions and exact source identities."""
import csv,json,hashlib
from pathlib import Path

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
c=json.loads(Path('config/ecology_species_source_review_20260927.json').read_text())
for p,h in c['pins'].items():assert sha(p)==h
r=json.loads(Path(c['receipt']).read_text())
for p,h in r['artifacts'].items():assert sha(p)==h
old=read(c['previous']);new=read(c['output']);source={r['species_name']:r for r in read(c['source'])};decisions=read(c['dispositions'])
assert new[:len(old)]==old and len(new)==45 and len({r['taxon_id'] for r in new})==45
assert len(decisions)==15 and {r['species_name'] for r in decisions}==set(c['decisions'])
accepted={r['species_name'] for r in decisions if r['disposition']=='retain_published_species_classification'}
assert {r['species_name'] for r in new[len(old):]}==accepted and len(accepted)==13
for row in new[len(old):]:
    s=source[row['species_name']];d=c['decisions'][row['species_name']]
    for a,b in [('taxon_id','taxon_id'),('assembly_accession','selected_assembly_accession'),('source_species_name','source_species_name'),('source_doi','source_doi'),('source_url','source_url'),('genus_candidate_primary_lifestyle','genus_candidate_primary_lifestyle'),('conflicts_with_genus_ECM_candidate','genus_ECM_disagreement')]:assert row[a]==s[b]
    assert row['state']==d['state'] and s['source_ecology']==d['source_label']
    assert row['trait']=='ecological_guild_as_classified_by_source' and row['selected_isolate_experimentally_verified']=='False'
    assert row['provisional_transition_group']=='unassigned_'+row['taxon_id']
    assert 'No binary ECM absence inferred' in row['scope_note'] and 'pending_' in row['confirmatory_test_status']
    assert 'Excel row '+s['source_excel_row'] in row['source_locator']
for row in decisions:
    d=c['decisions'][row['species_name']];s=source[row['species_name']]
    assert row['source_label']==s['source_ecology']==d['source_label'] and row['disposition']==d['disposition'] and row['retained_state']==d.get('state','') and row['note']==d['note']
result={'status':'passed_reviewed_species_ecology_preservation_and_source_joins','producer_receipt_sha256':sha(c['receipt']),'preserved_statements':len(old),'added_statements':len(accepted),'reviewed_dispositions':len(decisions),'script_sha256':sha(__file__),'scope':'All original rows preserved and new identities/source fields/states/dispositions reconstructed from independently audited imported records and explicit review configuration. Checks transcription and source linkage, not biological truth or isolate phenotype.'}
Path('metadata/ecology_species_source_review_readback_20260927.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
