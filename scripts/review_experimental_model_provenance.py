#!/usr/bin/env python3
"""Flag entry-level starting models and chronology without asserting training independence."""
import argparse,datetime,json
from collections import Counter
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha,read_table
from prepare_paired_phylogenetic_inputs import write_table


def date(value):
    if not value or value=='unknown':return None
    return datetime.date.fromisoformat(value[:10])


def classify(records):
    if not records:return 'starting_model_annotation_missing'
    if any('alphafold' in str(x.get('source_name','')).lower() for x in records):return 'explicit_AlphaFold_starting_model_entry_level'
    if any(x.get('type','').lower()=='in silico model' for x in records):return 'other_computational_starting_model_entry_level'
    if all(x.get('type','').lower()=='experimental model' for x in records):return 'only_experimental_starting_models_reported_not_independence_proof'
    return 'other_or_incomplete_starting_model_annotation'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['references','predictions','output']:p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use new immutable provenance review')
    checked_receipt(a.references);checked_receipt(a.predictions)
    refs=read_table(a.references/'entries.tsv');raw={r['entry_id']:r for r in json.loads((a.references/'method_specific_metadata.json').read_text())}
    models={r['model_id']:r for r in json.loads((a.predictions/'model_provenance.json').read_text())};entries=[];lookup={}
    for r in refs:
        if r['methodology']!='experimental':continue
        record=raw[r['entry_id']]['fields'].get('pdbx_initial_refinement_model',[]);category=classify(record)
        row={'entry_id':r['entry_id'],'experimental_method':r['methods'],'starting_model_review':category,'starting_models_json':json.dumps(record,sort_keys=True),'target_chain_attribution':'unknown_requires_methods_or_coordinate_annotation_review','training_independence':'unresolved','initial_release_date':r['initial_release_date']}
        entries.append(row);lookup[r['entry_id']]=row
    chronology=[]
    for r in read_table(a.references/'entities.tsv'):
        if r['entry_id'] not in lookup:continue
        model=models[r['model_id']]
        if model['sequence_sha256']!=r['model_sequence_sha256']:raise ValueError('Prediction identity mismatch')
        released=date(lookup[r['entry_id']]['initial_release_date']);created=date(model.get('model_created_date'))
        status='missing_date' if released is None or created is None else ('released_before_prediction_creation' if released<created else 'same_day_as_prediction_creation' if released==created else 'released_after_prediction_creation')
        chronology.append({'entry_id':r['entry_id'],'entity_id':r['entity_id'],'model_id':r['model_id'],'sequence_sha256':r['model_sequence_sha256'],'entry_initial_release_date':lookup[r['entry_id']]['initial_release_date'],'prediction_model_created_date':model.get('model_created_date','unknown'),'chronology':status,'starting_model_review':lookup[r['entry_id']]['starting_model_review'],'training_independence':'unresolved_release_and_creation_dates_are_not_training_or_template_cutoffs'})
    a.output.mkdir(parents=True);write_table(a.output/'entry_starting_model_review.tsv',entries);write_table(a.output/'entity_prediction_chronology.tsv',chronology)
    result={'status':'complete_entry_starting_model_and_prediction_chronology_review','reference_receipt_sha256':sha(a.references/'receipt.json'),'prediction_receipt_sha256':sha(a.predictions/'receipt.json'),'script_sha256':sha(Path(__file__)),'experimental_entries':len(entries),'entry_review_counts':dict(Counter(r['starting_model_review'] for r in entries)),'entity_model_rows':len(chronology),'entity_model_chronology_counts':dict(Counter(r['chronology'] for r in chronology)),'interpretation':'Starting models describe whole entries, not necessarily the exact matched target chain. Missing annotations do not establish absence of prediction-assisted refinement. Reported experimental starting models may themselves have dependencies. Prediction creation dates are not training/template cutoffs; later entry releases do not prove unseen sequences, structures or homologs. All references retain unresolved independence.','artifacts':{f.name:sha(f) for f in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
