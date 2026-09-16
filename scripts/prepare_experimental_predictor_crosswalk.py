#!/usr/bin/env python3
"""Bind experimental controls to original reference targets by exact sequence identity."""
import argparse
import hashlib
import json
from pathlib import Path
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Data.IUPACData import protein_letters_3to1
from audit_joint_path_uncertainty import checked, rows, sha
from prepare_paired_phylogenetic_inputs import write_table


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['controls','control-audit','predictions','screen','references','mapping','output']:
        ap.add_argument('--'+name,type=Path,required=True)
    a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    receipts={name:checked(getattr(a,name)) for name in ['controls','control_audit','predictions','screen','references']}
    if receipts['controls']['source_audit_receipt_sha256']!=sha(a.control_audit/'receipt.json') or receipts['screen']['mapping_receipt_sha256']!=sha(a.predictions/'receipt.json'):
        raise ValueError('Prediction lineage differs')
    mr=json.loads((a.mapping/'receipt.json').read_text());mc=json.loads((a.mapping/'config.json').read_text())
    if mr['status']!='complete_experimental_CA_correspondence_screen' or mr['config_sha256']!=sha(a.mapping/'config.json') or mc['screen_receipt_sha256']!=sha(a.screen/'receipt.json'):
        raise ValueError('Coordinate mapping lineage differs')
    models=json.loads((a.controls/'model_provenance.json').read_text())
    controls={m['sequence_sha256']:m for m in models}
    if len(controls)!=len(models) or len(models)!=receipts['controls']['models']:
        raise ValueError('Ambiguous control sequences')
    afrows=json.loads((a.predictions/'model_provenance.json').read_text());af={m['model_id']:m for m in afrows}
    if len(af)!=len(afrows):raise ValueError('Duplicate reference prediction identifier')
    refs={(r['entity_id'],r['model_id']):r for r in rows(a.references/'entities.tsv')}
    screened={(r['entity_id'],r['model_id']):r for r in rows(a.screen/'sequence_correspondence.tsv')}
    links=list(rows(a.control_audit/'experimental_reference_links.tsv'))
    if len(links)!=receipts['control_audit']['experimental_reference_links']:
        raise ValueError('Control link count differs')
    base=Path(__file__).resolve().parents[1];three={k.upper():v for k,v in protein_letters_3to1.items()};verified={}
    def sequence(model):
        name=model['model_id']
        if name in verified:return verified[name]
        p=base/model['path']
        if sha(p)!=model['sha256']:raise ValueError('Prediction file changed')
        d=MMCIF2Dict(str(p));ca={}
        for atom,pos,aa in zip(d['_atom_site.label_atom_id'],d['_atom_site.label_seq_id'],d['_atom_site.label_comp_id']):
            if atom=='CA':
                i=int(pos)
                if i in ca:raise ValueError('Ambiguous predicted residue')
                ca[i]=three.get(aa,'?')
        if sorted(ca)!=list(range(1,model['length']+1)):
            raise ValueError('Incomplete predicted sequence')
        seq=''.join(ca[i] for i in sorted(ca))
        if hashlib.sha256(seq.encode()).hexdigest()!=model['sequence_sha256']:
            raise ValueError('Coordinate sequence differs')
        verified[name]=seq;return seq
    output=[];seen=set();used=set()
    for link in links:
        key=link['entity_id'],link['predicted_model_id'];digest=link['sequence_sha256']
        if key in seen:raise ValueError('Duplicate reference link')
        seen.add(key);s=screened[key];ref=refs[key];control=controls[digest];prediction=af[key[1]]
        if s['sequence_class']!='exact_full_sequence' or s['model_sequence_sha256']!=digest or s['entity_sequence_sha256']!=digest or ref['model_sequence_sha256']!=digest or prediction['sequence_sha256']!=digest:
            raise ValueError('Nonexact reference correspondence')
        if link['sequence_id']!='S'+digest or ref['entry_id']!=link['entry_id'] or sequence(control)!=sequence(prediction):
            raise ValueError('Sequence or entry identity mismatch')
        used.add(digest)
        output.append({'entry_id':link['entry_id'],'entity_id':link['entity_id'],'sequence_sha256':digest,'sequence_length':control['length'],'alphafold_model_id':prediction['model_id'],'esmfold_model_id':control['model_id'],'alphafold_path':prediction['path'],'alphafold_sha256':prediction['sha256'],'esmfold_path':control['path'],'esmfold_sha256':control['sha256'],'coordinate_mapping_available':link['entry_id'] in mr['entry_receipt_sha256'],'sequence_class':'exact_full_sequence'})
    if used!=set(controls):raise ValueError('Unlinked control sequence')
    a.output.mkdir(parents=True);write_table(a.output/'crosswalk.tsv',output)
    result={'status':'complete_exact_sequence_experimental_predictor_crosswalk','control_sequences':len(controls),'alphafold_models':len({r['alphafold_model_id'] for r in output}),'reference_links':len(output),'entries':len({r['entry_id'] for r in output}),'links_with_coordinate_mapping':sum(r['coordinate_mapping_available'] for r in output),'coordinate_sequences_checked':len(verified),'source_receipts':{name:sha(getattr(a,name)/'receipt.json') for name in ['controls','control_audit','predictions','screen','references','mapping']},'script_sha256':sha(Path(__file__)),'artifacts':{'crosswalk.tsv':sha(a.output/'crosswalk.tsv')},'interpretation':'Exact canonical sequence correspondence with both prediction coordinate sequences re-read. All linked entries retained; no agreement-based selection. Coordinate masks, matched geometry, experimental context and training/template independence require separate checks. Entries/entities are not independent protein replicates.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
