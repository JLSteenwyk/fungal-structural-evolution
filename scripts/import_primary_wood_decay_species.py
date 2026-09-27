"""Preserve complete reviewed Figure 1 classifications and exact-name sample links."""
import json
from pathlib import Path
import pandas as pd
from screen_duplication_domain_alignment_coverage import sha

cp=Path('config/wood_decay_primary_species_20260927.json');c=json.loads(cp.read_text())
assert sha(c['source_pdf'])==c['source_sha256'] and sha(c['manifest'])==c['manifest_sha256']
source=pd.DataFrame(c['source_classifications'])
assert len(source)==source.source_species_name.nunique()==22
assert source.state.value_counts().to_dict()==dict(white_rot=12,brown_rot=7,uncertain_decay_mode=3)
manifest=pd.read_csv(c['manifest'],sep='\t',keep_default_na=False)
assert manifest.taxon_id.is_unique
print('Manifest name column candidates:',[x for x in manifest if 'name' in x])
rows=[]
for r in source.itertuples(index=False):
    hits=manifest[(manifest.species_name==r.source_species_name)&(manifest.study_role=='ingroup')]
    status='exact_unique_selected_name' if len(hits)==1 else ('no_exact_selected_name' if len(hits)==0 else 'ambiguous_selected_name')
    row=dict(source_species_name=r.source_species_name,state=r.state,match_status=status,selected_matches=len(hits),taxon_id='',assembly_accession='',source_doi=c['source_doi'],source_locator=c['source_locator'],trait='wood_decay_mode_as_classified_by_source',selected_isolate_verified=False,ecm_status='not_inferred')
    if len(hits)==1:row.update(taxon_id=hits.iloc[0].taxon_id,assembly_accession=hits.iloc[0].assembly_accession)
    rows.append(row)
df=pd.DataFrame(rows);out=Path('metadata/wood_decay_primary_species_20260927.tsv')
assert not out.exists();df.to_csv(out,sep='\t',index=False)
back=pd.read_csv(out,sep='\t',keep_default_na=False);pd.testing.assert_frame_equal(back,df,check_dtype=False)
# Reconstruct exact-name multiplicities independently by scalar source matching.
records=manifest.to_dict('records')
for row in back.to_dict('records'):
    expected=[m for m in records if m['species_name']==row['source_species_name'] and m['study_role']=='ingroup']
    assert len(expected)==row['selected_matches']
    assert row['taxon_id']==(expected[0]['taxon_id'] if len(expected)==1 else '')
    assert row['assembly_accession']==(expected[0]['assembly_accession'] if len(expected)==1 else '')
receipt=dict(status='complete_primary_wood_decay_source_transcription_and_exact_name_join',source_rows=len(df),match_counts=df.match_status.value_counts().to_dict(),matched_state_counts=df[df.match_status.eq('exact_unique_selected_name')].state.value_counts().to_dict(),config_sha256=sha(cp),source_pdf_sha256=c['source_sha256'],manifest_sha256=c['manifest_sha256'],table_sha256=sha(out),script_sha256=sha(__file__),scope=c['interpretation'])
Path('metadata/wood_decay_primary_species_receipt_20260927.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
