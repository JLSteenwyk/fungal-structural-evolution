"""Exercise gains, losses, replacements, duplicates and sequence mismatches."""
import json
from pathlib import Path
import pandas as pd
from compare_whole_proteome_catalogs import compare_links
from readback_whole_proteome_catalog import sha

def row(protein,model):
    return dict(taxon_id='fixture',protein_id=protein,sequence_sha256=protein,model_id=model,version='1',model_path=model+'.cif')
a=pd.DataFrame([row('same','a'),row('replace','b'),row('lost','c')])
b=pd.DataFrame([row('same','a'),row('replace','d'),row('new','e')])
r=compare_links(a,b)
assert r.set_index('protein_id').disposition.to_dict()==dict(same='unchanged_model',replace='changed_selected_model',lost='lost_catalog_link',new='new_catalog_link')
for bad in [pd.concat([b,b.iloc[[0]]],ignore_index=True),b.assign(sequence_sha256='changed')]:
    try:compare_links(a,bad)
    except AssertionError:pass
    else:raise AssertionError('Invalid source accepted')
p=Path('metadata/whole_proteome_catalog_comparison_checks_20260928.json')
p.write_text(json.dumps(dict(status='passed_gain_loss_replacement_and_invalid_input_fixtures',
    dispositions=4,duplicate_rejected=True,changed_sequence_rejected=True,
    pins={str(f):sha(f) for f in [Path(__file__),Path('scripts/compare_whole_proteome_catalogs.py')]}),indent=2)+'\n')
print('Passed catalog comparison fixtures')
