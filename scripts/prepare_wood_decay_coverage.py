"""Prepare source-specific coverage plans for every exact-name decay classification."""
import json
from pathlib import Path
import pandas as pd
from screen_duplication_domain_alignment_coverage import sha

source=Path('metadata/wood_decay_primary_species_20260927.tsv');rp=Path('metadata/wood_decay_primary_species_receipt_20260927.json');r=json.loads(rp.read_text())
assert r['status']=='complete_primary_wood_decay_source_transcription_and_exact_name_join' and sha(source)==r['table_sha256']
d=pd.read_csv(source,sep='\t',keep_default_na=False);d=d[d.match_status.eq('exact_unique_selected_name')].copy();assert len(d)==7
# Group labels describe shared source categories, never independent transitions.
d['species_name']=d.source_species_name;d['provisional_transition_group']='source_category_'+d.state
p=Path('metadata/wood_decay_coverage_evidence_20260927.tsv');assert not p.exists();d.to_csv(p,sep='\t',index=False)
er=Path('metadata/wood_decay_coverage_evidence_receipt_20260927.json');er.write_text(json.dumps(dict(status='complete_exact_name_wood_decay_coverage_input',taxa=7,table_sha256=sha(p),source_receipt_sha256=sha(rp),scope='Source categories retained, including uncertain decay modes. Groups describe coverage only, not origins or biological replicates.'),indent=2)+'\n')
for source in ['afdb','esmfold']:
    template=json.loads(Path(f'metadata/qualified_{source}_reviewed45_ecology_plan_20260927.json').read_text())
    template.update(evidence=str(p),evidence_receipt=str(er),output=f'results/ecology/qualified-{source}-wood-decay-overlap-20260927-v1')
    pins={f:sha(f) for f in [str(p),str(er),str(rp),'metadata/wood_decay_primary_species_20260927.tsv','config/wood_decay_primary_species_20260927.json',template['input_readback'],str(Path(template['inputs'])/'receipt.json'),'scripts/readback_qualified_ecology_overlap.py','scripts/readback_whole_proteome_family_coverage.py','scripts/prepare_wood_decay_coverage.py']}
    script='scripts/assess_qualified_afdb_ecology_overlap.py' if source=='afdb' else 'scripts/assess_qualified_ecology_overlap.py';pins[script]=sha(script)
    template['pins']=pins
    Path(f'metadata/qualified_{source}_wood_decay_plan_20260927.json').write_text(json.dumps(template,indent=2)+'\n')
print('Prepared seven-taxon plans: one CPU, 4 GiB, 1–10 minutes per source, no GPU')
