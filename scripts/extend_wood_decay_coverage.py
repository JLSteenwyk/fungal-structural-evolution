"""Add the independently reviewed named Dacryopinax classification, retaining v1."""
import json
from pathlib import Path
import pandas as pd
from screen_duplication_domain_alignment_coverage import sha

source=Path('metadata/wood_decay_coverage_evidence_20260927.tsv');rp=Path('metadata/wood_decay_coverage_evidence_receipt_20260927.json')
r=json.loads(rp.read_text());assert sha(source)==r['table_sha256']
reviewpath=Path('metadata/wood_decay_name_review_20260927.json');review=json.loads(reviewpath.read_text())
assert review['status']=='complete_targeted_wood_decay_name_review'
assert sha('metadata/analysis_manifest.tsv')==review['manifest_sha256']
accepted=[x for x in review['decisions'] if x['decision']=='independent_named_species_primary_genome_description_available']
assert len(accepted)==1;addition=accepted[0]
assert addition['taxon_id']=='F1858805' and addition['state']=='brown_rot'
d=pd.read_csv(source,sep='\t',keep_default_na=False);assert len(d)==7 and addition['taxon_id'] not in set(d.taxon_id)
row={k:'' for k in d.columns}
row.update(source_species_name=addition['selected_species'],state=addition['state'],match_status='independent_named_species_genome_description',selected_matches=1,taxon_id=addition['taxon_id'],assembly_accession=addition['selected_assembly'],source_doi='',source_locator=review['sources']['dacryopinax_jgi']['url']+'; genome description',trait='wood_decay_mode_as_classified_by_source',selected_isolate_verified=False,ecm_status='not_inferred',species_name=addition['selected_species'],provisional_transition_group='source_category_brown_rot')
expanded=pd.concat([d,pd.DataFrame([row])],ignore_index=True);assert expanded.taxon_id.is_unique
pd.testing.assert_frame_equal(expanded.iloc[:7].reset_index(drop=True),d)
p=Path('metadata/wood_decay_coverage_evidence_v2_20260927.tsv');assert not p.exists();expanded.to_csv(p,sep='\t',index=False)
pd.testing.assert_frame_equal(pd.read_csv(p,sep='\t',keep_default_na=False),expanded)
er=Path('metadata/wood_decay_coverage_evidence_v2_receipt_20260927.json')
er.write_text(json.dumps(dict(status='complete_expanded_primary_wood_decay_coverage_input',taxa=8,table_sha256=sha(p),original_table_sha256=sha(source),source_review_sha256=sha(reviewpath),scope='Seven original rows preserved exactly; one independent named-species primary genome classification added. No inferred ecological origins or binary ECM coding.'),indent=2)+'\n')
for source_name in ['afdb','esmfold']:
    old=Path(f'metadata/qualified_{source_name}_wood_decay_plan_20260927.json');plan=json.loads(old.read_text())
    plan.update(evidence=str(p),evidence_receipt=str(er),output=f'results/ecology/qualified-{source_name}-wood-decay-overlap-20260927-v2')
    for f in [str(p),str(er),str(old),str(reviewpath),'scripts/extend_wood_decay_coverage.py']:plan['pins'][f]=sha(f)
    Path(f'metadata/qualified_{source_name}_wood_decay_v2_plan_20260927.json').write_text(json.dumps(plan,indent=2)+'\n')
print('Prepared eight-taxon source plans, retaining all prior evidence')
