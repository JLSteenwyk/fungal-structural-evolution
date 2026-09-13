#!/usr/bin/env python3
"""Map experimental CA observations without collapsing missing/alternate/modified residues."""
import argparse,csv,fcntl,gzip,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Data.IUPACData import protein_letters_3to1
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha,read_table
THREE={k.upper():v for k,v in protein_letters_3to1.items()}


def classify_atoms(atoms,aa):
    if not atoms:return 'no_CA_observation'
    if len(atoms)!=1:return 'multiple_CA_records'
    atom=atoms[0]
    if THREE.get(atom['label_comp_id'])!=aa:return 'nonstandard_or_mismatching_monomer'
    try:values=[float(atom[k]) for k in ['Cartn_x','Cartn_y','Cartn_z','occupancy']]
    except ValueError:return 'invalid_coordinates_or_occupancy'
    if not all(math.isfinite(x) for x in values) or not 0<values[3]<=1:return 'invalid_coordinates_or_occupancy'
    if atom['label_alt_id'] not in ('.','?') or values[3]!=1:return 'alternate_or_partial_occupancy'
    return 'unambiguous_full_occupancy_CA'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['screen','coordinates','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();checked_receipt(a.screen);cr=json.loads((a.coordinates/'receipt.json').read_text());cc=json.loads((a.coordinates/'config.json').read_text())
    if cr['config_sha256']!=sha(a.coordinates/'config.json') or cc['screen_receipt_sha256']!=sha(a.screen/'receipt.json'):raise ValueError('Coordinate/screen lineage differs')
    selected=defaultdict(list)
    for r in read_table(a.screen/'sequence_correspondence.tsv'):
        if r['sequence_class']=='exact_full_sequence':selected[r['entity_id'].rsplit('_',1)[0]].append(r)
    if set(selected)!={r['entry_id'] for r in cr['results']}:raise ValueError('Coordinate entry set differs')
    a.output.mkdir(parents=True,exist_ok=True);lock=(a.output/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    config={'screen_receipt_sha256':sha(a.screen/'receipt.json'),'coordinate_receipt_sha256':sha(a.coordinates/'receipt.json'),'script_sha256':sha(Path(__file__)),'policy':'All deposited model IDs and selected struct_asym chains; no alternate-location selection; full sequence positions retained. CA status is a correspondence screen, not experimental quality.'}
    cp=a.output/'config.json'
    if cp.exists() and json.loads(cp.read_text())!=config:raise ValueError('Changed mapping configuration')
    cp.write_text(json.dumps(config,indent=2)+'\n');summaries=[]
    for item in cr['results']:
        entry=item['entry_id'];rp=a.output/(entry+'.receipt.json');out=a.output/(entry+'.residues.tsv.gz')
        if rp.exists():
            r=json.loads(rp.read_text())
            if r['config_sha256']!=sha(cp) or r['table_sha256']!=sha(out):raise ValueError('Changed cached mapping')
            summaries.append(r);continue
        path=a.coordinates/(entry+'.cif.gz')
        if sha(path)!=item['gzip_sha256']:raise ValueError('Changed coordinate archive')
        with gzip.open(path,'rt') as f:d=MMCIF2Dict(f)
        seqs={i:''.join(s.split()) for i,s in zip(d['_entity_poly.entity_id'],d['_entity_poly.pdbx_seq_one_letter_code_can'])}
        targets={}
        for row in selected[entry]:
            eid=row['entity_id'].rsplit('_',1)[1];seq=seqs[eid]
            if hashlib.sha256(seq.encode()).hexdigest()!=row['model_sequence_sha256']:raise ValueError('Coordinate polymer sequence differs from exact target')
            targets[eid]=seq
        chains={c:e for c,e in zip(d['_struct_asym.id'],d['_struct_asym.entity_id']) if e in targets}
        names=['label_atom_id','label_entity_id','label_asym_id','label_seq_id','label_comp_id','label_alt_id','pdbx_PDB_model_num','Cartn_x','Cartn_y','Cartn_z','occupancy','B_iso_or_equiv','auth_asym_id','auth_seq_id','pdbx_PDB_ins_code']
        columns=[d['_atom_site.'+n] for n in names]
        if len({len(c) for c in columns})!=1:raise ValueError('Nonrectangular atom table')
        observations=defaultdict(list);model_ids=set(d['_atom_site.pdbx_PDB_model_num'])
        for values in zip(*columns):
            atom=dict(zip(names,values));eid=atom['label_entity_id']
            if eid not in targets or atom['label_atom_id']!='CA':continue
            pos=int(atom['label_seq_id'])
            if chains.get(atom['label_asym_id'])!=eid or not 1<=pos<=len(targets[eid]):raise ValueError('CA entity/position mismatch')
            observations[atom['pdbx_PDB_model_num'],atom['label_asym_id'],pos].append(atom)
        counts=Counter();n=0
        with gzip.open(out,'wt') as f:
            writer=None
            for model in sorted(model_ids,key=int):
                for chain,eid in sorted(chains.items()):
                    for pos,aa in enumerate(targets[eid],1):
                        atoms=observations.get((model,chain,pos),[]);status=classify_atoms(atoms,aa);counts[status]+=1;n+=1
                        row={'entry_id':entry,'entity_id':eid,'model_number':model,'label_asym_id':chain,'label_seq_id':pos,'target_aa':aa,'CA_status':status,'CA_record_count':len(atoms),'atom_records_json':json.dumps(atoms,separators=(',',':'))}
                        if writer is None:writer=csv.DictWriter(f,list(row),delimiter='\t',lineterminator='\n');writer.writeheader()
                        writer.writerow(row)
        r={'entry_id':entry,'config_sha256':sha(cp),'table_sha256':sha(out),'deposited_models':len(model_ids),'selected_chains':len(chains),'position_rows':n,'status_counts':dict(counts)}
        rp.write_text(json.dumps(r,indent=2)+'\n');summaries.append(r);print(entry,n,dict(counts),flush=True)
    totals=Counter()
    for r in summaries:totals.update(r['status_counts'])
    result={'status':'complete_experimental_CA_correspondence_screen','config_sha256':sha(cp),'entries':len(summaries),'position_rows':sum(r['position_rows'] for r in summaries),'status_counts':dict(totals),'entry_receipt_sha256':{r['entry_id']:sha(a.output/(r['entry_id']+'.receipt.json')) for r in summaries},'interpretation':'Deposited exact canonical sequences checked against coordinate entity sequences. Every model/selected chain/full-sequence position retained, including missing or ambiguous CA observations. Raw B factors retained only as experimental annotations, never pLDDT. Experimental method/quality, other atoms, biological assemblies and training independence remain separate requirements.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
