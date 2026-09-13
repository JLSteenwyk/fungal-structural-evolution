#!/usr/bin/env python3
"""Export an audited ESMFold snapshot as sequence-explicit mmCIF and local inventory."""
import argparse
import io
import json
from pathlib import Path
import numpy as np
from Bio.PDB import MMCIFIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT,sha,read_table


def cif_text(pdb, sequence, name):
    atoms=[line for line in pdb.splitlines() if line.startswith('ATOM  ')]
    if any(line.startswith('HETATM') for line in pdb.splitlines()):raise ValueError('Unexpected heteroatoms')
    ca=[line for line in atoms if line[12:16].strip()=='CA']
    if ''.join(seq1(line[17:20]) for line in ca)!=sequence or [int(line[22:26]) for line in ca]!=list(range(1,len(sequence)+1)) or len({line[21] for line in atoms})!=1:
        raise ValueError('PDB sequence, chain or residue numbering differs')
    if any(line[16]!=' ' or line[26]!=' ' for line in atoms):raise ValueError('Alternate locations or insertion codes require explicit handling')
    coordinates=np.array([[float(line[a:b]) for a,b in [(30,38),(38,46),(46,54),(60,66)]] for line in atoms])
    if not np.isfinite(coordinates).all() or (coordinates[:,3]<0).any() or (coordinates[:,3]>100).any():raise ValueError('Invalid PDB coordinates/confidence')
    d={'data_':name,'_entity.id':['1'],'_entity.type':['polymer'],'_entity_poly.entity_id':['1'],
       '_entity_poly.type':['polypeptide(L)'],'_entity_poly.pdbx_seq_one_letter_code_can':[sequence],
       '_entity_poly.pdbx_seq_one_letter_code':[sequence],'_entity_poly.pdbx_strand_id':[atoms[0][21]]}
    fields={'group_PDB':lambda x:'ATOM','id':lambda x:x[6:11].strip(),'type_symbol':lambda x:x[76:78].strip(),
            'label_atom_id':lambda x:x[12:16].strip(),'label_alt_id':lambda x:'.','label_comp_id':lambda x:x[17:20],
            'label_asym_id':lambda x:x[21],'label_entity_id':lambda x:'1','label_seq_id':lambda x:x[22:26].strip(),
            'pdbx_PDB_ins_code':lambda x:'?','Cartn_x':lambda x:x[30:38].strip(),'Cartn_y':lambda x:x[38:46].strip(),
            'Cartn_z':lambda x:x[46:54].strip(),'occupancy':lambda x:x[54:60].strip(),'B_iso_or_equiv':lambda x:x[60:66].strip(),
            'auth_seq_id':lambda x:x[22:26].strip(),'auth_asym_id':lambda x:x[21],'pdbx_PDB_model_num':lambda x:'1'}
    for field,extract in fields.items():d['_atom_site.'+field]=[extract(x) for x in atoms]
    writer=MMCIFIO();writer.set_dict(d);buffer=io.StringIO();writer.save(buffer);text=buffer.getvalue()
    restored=MMCIF2Dict(io.StringIO(text))
    for field in d:
        if restored[field]!=d[field]:raise ValueError('mmCIF roundtrip changed '+field)
    return text


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['audit','predictions','output']:p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable local structure snapshot')
    audit=checked_receipt(a.audit)
    if audit['status'] not in ['complete_artifact_readback','complete_artifact_readback_of_partial_prediction_snapshot']:raise ValueError('Independent completed artifact audit required')
    config_path=a.predictions/'config.json';config_digest=sha(config_path)
    if audit['config_sha256']!=config_digest:raise ValueError('Different prediction configuration')
    rows=read_table(a.audit/'predictions.tsv')
    if len(rows)!=audit['predictions'] or len({r['sequence_id'] for r in rows})!=len(rows):raise ValueError('Audited prediction universe differs')
    a.output=a.output.resolve();a.output.mkdir(parents=True);folder=a.output/'models';folder.mkdir()
    inventory=[];provenance=[]
    for row in rows:
        sid=row['sequence_id'];receipt_path=a.predictions/(sid+'.json')
        if sha(receipt_path)!=row['prediction_receipt_sha256']:raise ValueError('Changed audited prediction receipt')
        r=json.loads(receipt_path.read_text())
        if r['config_sha256']!=config_digest or r['status']!='verified_prediction':raise ValueError('Prediction configuration/status differs')
        for name,h in r['artifacts'].items():
            if sha(a.predictions/name)!=h:raise ValueError('Changed prediction artifact')
        with np.load(a.predictions/(sid+'.npz'),allow_pickle=False) as data:sequence=str(data['sequence'])
        model_id='ESM-'+sid+'-'+config_digest[:12];target=folder/(model_id+'-v1.cif')
        target.write_text(cif_text((a.predictions/(sid+'.pdb')).read_text(),sequence,model_id))
        model={'model_id':model_id,'version':1,'sequence_sha256':r['sequence_sha256'],'length':r['length'],
               'mean_ca_plddt':r['mean_ca_plddt'],'fraction_ca_plddt_below50':r['fraction_ca_plddt_below50'],
               'provider':'local','tool':'ESMFold v1','path':str(target.relative_to(ROOT)),'sha256':sha(target),
               'prediction_config_sha256':config_digest,'prediction_receipt_path':str(receipt_path.resolve().relative_to(ROOT)),
               'prediction_receipt_sha256':sha(receipt_path),'local_pae_npz_path':str((a.predictions/(sid+'.npz')).resolve().relative_to(ROOT)),
               'local_pae_npz_sha256':r['artifacts'][sid+'.npz'],'representation':'mmCIF converted from audited PDB with all atom fields and coordinates preserved; version 1 is local representation version, not AFDB release.'}
        inventory.append({'record_id':sid,'status':'verified','models':[model]});provenance.append(model)
        if len(inventory)%500==0:print('Converted',len(inventory),'of',len(rows),flush=True)
    path=a.output/'inventory.jsonl';path.write_text(''.join(json.dumps(r)+'\n' for r in inventory))
    (a.output/'model_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    result={'status':'complete_audited_local_model_conversion','models':len(rows),'source_audit_receipt_sha256':sha(a.audit/'receipt.json'),
        'prediction_config_sha256':config_digest,'script_sha256':sha(Path(__file__)),
        'interpretation':'Local ESMFold snapshot, kept separate from GDM models. Converted all audited model atom records and canonical sequence without coordinate fitting or relaxation. Mapping, local PAE binding and native feature qualification remain pending.',
        'artifacts':{'inventory.jsonl':sha(path),'model_provenance.json':sha(a.output/'model_provenance.json')}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
