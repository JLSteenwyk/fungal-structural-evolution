#!/usr/bin/env python3
"""Compare AF, ESMFold and experiment on identical residues and local distance pairs."""
import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Data.IUPACData import protein_letters_3to1
from scipy.spatial.distance import pdist
from scipy.spatial.transform import Rotation
from audit_joint_path_uncertainty import checked, rows, sha
from compare_marker_structures import geometry
from prepare_paired_phylogenetic_inputs import write_table


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['crosswalk','mapping','references','output']:
        ap.add_argument('--'+name,type=Path,required=True)
    a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    cr=checked(a.crosswalk);checked(a.references)
    for name in ['mapping','references']:
        if cr['source_receipts'][name]!=sha(getattr(a,name)/'receipt.json'):
            raise ValueError('Source lineage differs')
    mr=json.loads((a.mapping/'receipt.json').read_text())
    refs={r['entry_id']:r for r in rows(a.references/'entries.tsv')}
    links=list(rows(a.crosswalk/'crosswalk.tsv'));targets=defaultdict(list)
    for link in links:targets[link['entity_id']].append(link)
    root=Path(__file__).resolve().parents[1];three={k.upper():v for k,v in protein_letters_3to1.items()};cache={}
    def model(link,label):
        name=link[label+'_model_id']
        if name in cache:return cache[name]
        path=root/link[label+'_path']
        if sha(path)!=link[label+'_sha256']:raise ValueError('Changed predicted model')
        d=MMCIF2Dict(str(path));out={}
        for atom,pos,aa,x,y,z,b in zip(d['_atom_site.label_atom_id'],d['_atom_site.label_seq_id'],d['_atom_site.label_comp_id'],d['_atom_site.Cartn_x'],d['_atom_site.Cartn_y'],d['_atom_site.Cartn_z'],d['_atom_site.B_iso_or_equiv']):
            if atom=='CA':
                i=int(pos)
                if i in out:raise ValueError('Duplicate predicted CA')
                out[i]=(three.get(aa,'?'),np.array([float(x),float(y),float(z)]),float(b))
        if sorted(out)!=list(range(1,int(link['sequence_length'])+1)) or hashlib.sha256(''.join(out[i][0] for i in sorted(out)).encode()).hexdigest()!=link['sequence_sha256']:
            raise ValueError('Predicted sequence differs')
        if any(not np.isfinite(v[1]).all() or not math.isfinite(v[2]) or not 0<=v[2]<=100 for v in out.values()):raise ValueError('Invalid coordinates/confidence')
        cache[name]=out;return out
    accepted=[];excluded=[];seen=set();independent=0;chain_keys=set()
    for entry in sorted({r['entry_id'] for r in links}):
        if entry not in mr['entry_receipt_sha256']:continue
        rp=a.mapping/(entry+'.receipt.json');receipt=json.loads(rp.read_text());path=a.mapping/(entry+'.residues.tsv.gz')
        if sha(rp)!=mr['entry_receipt_sha256'][entry] or receipt['config_sha256']!=mr['config_sha256'] or sha(path)!=receipt['table_sha256']:
            raise ValueError('Changed experimental mapping')
        if refs[entry]['methodology']!='experimental':raise ValueError('Nonexperimental reference')
        groups=defaultdict(list)
        with gzip.open(path,'rt') as stream:
            for row in csv.DictReader(stream,delimiter='\t'):
                if entry+'_'+row['entity_id'] in targets:groups[row['entity_id'],row['model_number'],row['label_asym_id']].append(row)
        for (entity,number,chain),group in sorted(groups.items()):
            group.sort(key=lambda r:int(r['label_seq_id']));positions=[int(r['label_seq_id']) for r in group]
            if positions!=list(range(1,len(group)+1)):raise ValueError('Incomplete experimental residue grid')
            digest=hashlib.sha256(''.join(r['target_aa'] for r in group).encode()).hexdigest();observed={}
            for row in group:
                if row['CA_status']=='unambiguous_full_occupancy_CA':
                    atoms=json.loads(row['atom_records_json'])
                    if len(atoms)!=1:raise ValueError('Ambiguous accepted experimental CA')
                    observed[int(row['label_seq_id'])]=np.array([float(atoms[0][k]) for k in ['Cartn_x','Cartn_y','Cartn_z']])
            for link in targets[entry+'_'+entity]:
                if digest!=link['sequence_sha256'] or len(group)!=int(link['sequence_length']):raise ValueError('Experimental sequence differs')
                af=model(link,'alphafold');esm=model(link,'esmfold');identity=(entry,entity,number,chain,link['alphafold_model_id'],link['esmfold_model_id'])
                if identity in chain_keys:raise ValueError('Duplicate chain/model comparison')
                chain_keys.add(identity);seen.add((link['entity_id'],link['alphafold_model_id'],link['esmfold_model_id']))
                for cutoff in [0,70,90]:
                    selected=sorted(i for i in observed if min(af[i][2],esm[i][2])>=cutoff);n=len(selected)
                    row={'entry_id':entry,'entity_id':entity,'deposited_model':number,'label_asym_id':chain,'alphafold_model_id':link['alphafold_model_id'],'esmfold_model_id':link['esmfold_model_id'],'sequence_sha256':digest,'sequence_length':len(group),'joint_predicted_plddt_cutoff':cutoff,'experimental_unambiguous_CA':len(observed),'matched_residues':n,'fraction_full_sequence':n/len(group),'residue_positions_json':json.dumps(selected,separators=(',',':')),'experimental_method':refs[entry]['methods'],'resolution_combined_json':refs[entry]['resolution_combined_json'],'initial_release_date':refs[entry]['initial_release_date']}
                    if n<50 or n*2<len(group):excluded.append(row|{'reason':'fewer_than_50_CA_or_below_half_full_sequence'});continue
                    xyz=[np.array([af[i][1] for i in selected]),np.array([esm[i][1] for i in selected]),np.array([observed[i] for i in selected])];pos=np.array(selected)
                    distances=[pdist(x) for x in xyz];ii,jj=np.triu_indices(n,1);local=(np.minimum.reduce(distances)<=15)&(abs(pos[ii]-pos[jj])>=3)
                    row['common_local_distance_pairs']=int(local.sum())
                    for label,xidx,yidx in [('af_experiment',0,2),('esm_experiment',1,2),('af_esm',0,1)]:
                        x,y=xyz[xidx],xyz[yidx];rms=geometry(x,y,pos,pos)['ca_superposition_rmsd_angstrom'];xc=x-x.mean(0);yc=y-y.mean(0);rotation,_=Rotation.align_vectors(yc,xc)
                        check=float(np.sqrt(np.mean(np.sum((rotation.apply(xc)-yc)**2,axis=1))))
                        if not math.isclose(rms,check,rel_tol=1e-8,abs_tol=1e-8):raise ValueError('Independent superposition differs')
                        independent+=1;delta=distances[xidx][local]-distances[yidx][local]
                        row[label+'_rmsd_angstrom']=rms;row[label+'_local_mean_absolute_difference_angstrom']=float(np.mean(abs(delta))) if len(delta) else '';row[label+'_local_rms_difference_angstrom']=float(np.sqrt(np.mean(delta**2))) if len(delta) else ''
                    accepted.append(row)
        print(entry,'matched predictor controls compared',flush=True)
    unmapped=[r for r in links if (r['entity_id'],r['alphafold_model_id'],r['esmfold_model_id']) not in seen]
    if len(accepted)+len(excluded)!=3*len(chain_keys):raise ValueError('Incomplete threshold grid')
    a.output.mkdir(parents=True)
    for name,values in [('comparisons.tsv',accepted),('exclusions.tsv',excluded),('unmapped_links.tsv',unmapped)]:
        if values:write_table(a.output/name,values)
    result={'status':'complete_matched_experimental_predictor_control_geometry','reference_links':len(links),'links_with_chain_grids':len(seen),'links_without_chain_grids':len(unmapped),'chain_model_grids':len(chain_keys),'accepted_rows':len(accepted),'excluded_rows':len(excluded),'independent_superposition_checks':independent,'source_receipts':{name:sha(getattr(a,name)/'receipt.json') for name in ['crosswalk','mapping','references']},'script_sha256':sha(Path(__file__)),'geometry_helper_sha256':sha(Path(__file__).with_name('compare_marker_structures.py')),'interpretation':'Same observed residues for all three comparisons, using joint predicted focal pLDDT thresholds 0/70/90. Common distance-pair mask: separation >=3 and distance <=15A in any of the three structures. All experimental chains/models retained; no best-agreement selection. Thresholds alter residues/cohorts. Protein-weighted summary, complete independent mask/geometry readback and experimental/training/context qualification remain pending; not unbiased accuracy or evolutionary change.','artifacts':{p.name:sha(p) for p in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
