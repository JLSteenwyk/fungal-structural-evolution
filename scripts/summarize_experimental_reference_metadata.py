#!/usr/bin/env python3
"""Freeze method-specific experimental metadata without declaring benchmark eligibility."""
import argparse,json
from collections import Counter,defaultdict
from pathlib import Path
from audit_busco_gene_copies import sha,read_table
from assess_pae_sensitivity import checked_receipt
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--screen',type=Path,required=True)
    p.add_argument('--metadata',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable output')
    checked_receipt(a.screen)
    screen=json.loads((a.screen/'receipt.json').read_text())
    if screen['status']!='complete_entity_sequence_screen':raise ValueError('Complete sequence screen required')
    config_hash=sha(a.metadata/'config.json')
    if screen['metadata_config_sha256']!=config_hash:raise ValueError('Metadata lineage differs')
    retrieval=json.loads((a.metadata/'receipt.json').read_text())
    pins={(r['kind'],r['identifier']):r for r in retrieval['responses']}
    def response(kind,identifier):
        path=a.metadata/kind/(identifier+'.json');rp=path.with_suffix('.receipt.json')
        r=json.loads(rp.read_text());pin=pins[kind,identifier]
        if r['config_sha256']!=config_hash or r['response_sha256']!=sha(path):raise ValueError('Changed response')
        if pin['receipt_sha256']!=sha(rp):raise ValueError('Changed receipt')
        value=json.loads(path.read_text())
        if value['rcsb_id']!=identifier:raise ValueError('Response identity differs')
        return value
    rows=[r for r in read_table(a.screen/'sequence_correspondence.tsv') if r['sequence_class']=='exact_full_sequence']
    grouped=defaultdict(list)
    for r in rows:grouped[r['entity_id'].rsplit('_',1)[0]].append(r)
    entries=[];entities=[];raw=[]
    fields=['exptl','refine','em_3d_reconstruction','pdbx_nmr_ensemble','pdbx_vrpt_summary_geometry','pdbx_vrpt_summary_em','pdbx_vrpt_summary_diffraction','pdbx_vrpt_summary_nmr','pdbx_initial_refinement_model','em_3d_fitting','em_3d_fitting_list','rcsb_entry_info','rcsb_accession_info']
    for entry,matched in sorted(grouped.items()):
        d=response('entry',entry);info=d.get('rcsb_entry_info',{});acc=d.get('rcsb_accession_info',{})
        methods=sorted({x['method'] for x in d.get('exptl',[])})
        entries.append({'entry_id':entry,'methodology':info.get('structure_determination_methodology','unknown'),'methods':';'.join(methods) or 'unknown','resolution_combined_json':json.dumps(info.get('resolution_combined',None)),'initial_release_date':acc.get('initial_release_date','unknown'),'revision_date':acc.get('revision_date','unknown'),'experimental_data_released':acc.get('has_released_experimental_data','unknown'),'polymer_composition':info.get('polymer_composition','unknown'),'polymer_entities':info.get('polymer_entity_count','unknown'),'deposited_models':info.get('deposited_model_count','unknown'),'exact_entities':len({r['entity_id'] for r in matched}),'exact_predicted_proteins':len({r['model_id'] for r in matched}),'entry_response_sha256':sha(a.metadata/'entry'/(entry+'.json')),'benchmark_eligibility':'pending_method_specific_quality_coordinate_context_and_training_review'})
        raw.append({'entry_id':entry,'fields':{key:d[key] for key in fields if key in d},'absent_fields':[key for key in fields if key not in d]})
        for row in matched:
            entity=response('polymer_entity',row['entity_id'])
            if sha(a.metadata/'polymer_entity'/(row['entity_id']+'.receipt.json'))!=screen['entity_receipt_sha256'][row['entity_id']]:raise ValueError('Screen entity lineage differs')
            entities.append({**row,'entry_id':entry,'source_organisms_json':json.dumps(entity.get('rcsb_entity_source_organism',None),sort_keys=True),'host_organisms_json':json.dumps(entity.get('rcsb_entity_host_organism',None),sort_keys=True),'entity_description':entity.get('rcsb_polymer_entity',{}).get('pdbx_description','unknown')})
    a.output.mkdir(parents=True)
    write_table(a.output/'entries.tsv',entries);write_table(a.output/'entities.tsv',entities)
    (a.output/'method_specific_metadata.json').write_text(json.dumps(raw,indent=2)+'\n')
    result={'status':'complete_exact_sequence_reference_metadata_inventory','entries':len(entries),'entity_model_rows':len(entities),'distinct_predicted_proteins':len({r['model_id'] for r in rows}),'entry_methodology_counts':dict(Counter(r['methodology'] for r in entries)),'entry_method_counts':dict(Counter(r['methods'] for r in entries)),'entry_composition_counts':dict(Counter(r['polymer_composition'] for r in entries)),'entries_missing_resolution':sum(r['resolution_combined_json']=='null' for r in entries),'screen_receipt_sha256':sha(a.screen/'receipt.json'),'metadata_receipt_sha256':sha(a.metadata/'receipt.json'),'script_sha256':sha(Path(__file__)),'interpretation':'All exact full canonical sequence candidates retained without choosing on prediction agreement. Resolution is recorded with experimental method, never substituted for local residue reliability. Entry-level geometry and map metrics are not target-chain metrics. Missing fields remain unknown. Dates do not establish training independence; refinement starting models and homolog/template overlap require review. Repeated entries/chains are not independent proteins.','artifacts':{f.name:sha(f) for f in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
