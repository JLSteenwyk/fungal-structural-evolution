#!/usr/bin/env python3
"""Reconstruct every common residue/pair mask and all three control geometries."""
import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1
from Bio.SVDSuperimposer import SVDSuperimposer
from scipy.spatial.distance import cdist
from audit_joint_path_uncertainty import checked, rows, sha


def close(a,b):
    if not math.isclose(float(a),float(b),rel_tol=1e-7,abs_tol=1e-6):
        raise ValueError(f'Numerical disagreement: {a}, {b}')


def identity(row):
    return tuple(row[k] for k in ['entry_id','entity_id','deposited_model','label_asym_id','alphafold_model_id','esmfold_model_id','joint_predicted_plddt_cutoff'])


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['comparisons','crosswalk','mapping','references','output']:
        ap.add_argument('--'+name,type=Path,required=True)
    a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    r=checked(a.comparisons);checked(a.crosswalk);checked(a.references)
    if r['status']!='complete_matched_experimental_predictor_control_geometry':raise ValueError('Incomplete producer')
    for name in ['crosswalk','mapping','references']:
        if sha(getattr(a,name)/'receipt.json')!=r['source_receipts'][name]:raise ValueError('Source changed')
    stored={};counts={True:0,False:0}
    for accepted,name in [(True,'comparisons.tsv'),(False,'exclusions.tsv')]:
        p=a.comparisons/name
        for row in rows(p) if p.exists() else []:
            key=identity(row)
            if key in stored:raise ValueError('Duplicate output identity')
            stored[key]=(accepted,row);counts[accepted]+=1
    if counts[True]!=r['accepted_rows'] or counts[False]!=r['excluded_rows']:raise ValueError('Output counts differ')
    mr=json.loads((a.mapping/'receipt.json').read_text());refs={x['entry_id']:x for x in rows(a.references/'entries.tsv')}
    links=list(rows(a.crosswalk/'crosswalk.tsv'));targets=defaultdict(list)
    for link in links:targets[link['entity_id']].append(link)
    cache={};root=Path(__file__).resolve().parents[1];parser=MMCIFParser(QUIET=True,auth_chains=False,auth_residues=False)
    def model(link,label):
        name=link[label+'_model_id']
        if name in cache:return cache[name]
        path=root/link[label+'_path']
        if sha(path)!=link[label+'_sha256']:raise ValueError('Changed coordinate file')
        structure=parser.get_structure(name,str(path));models=list(structure);chains=list(structure.get_chains());residues=list(structure.get_residues())
        if len(models)!=1 or len(chains)!=1 or [x.id[1] for x in residues]!=list(range(1,int(link['sequence_length'])+1)):
            raise ValueError('Predicted sequence numbering differs')
        sequence=''.join(seq1(x.resname) for x in residues)
        if hashlib.sha256(sequence.encode()).hexdigest()!=link['sequence_sha256']:raise ValueError('Predicted canonical sequence differs')
        # MMCIFParser stores float32 coordinates. Read decimal fields in float64
        # so rounding cannot move an exactly 15-A pair across the cutoff.
        d=MMCIF2Dict(str(path));data={}
        for atom,pos,x,y,z,b in zip(d['_atom_site.label_atom_id'],d['_atom_site.label_seq_id'],d['_atom_site.Cartn_x'],d['_atom_site.Cartn_y'],d['_atom_site.Cartn_z'],d['_atom_site.B_iso_or_equiv']):
            if atom=='CA':
                i=int(pos)
                if i in data:raise ValueError('Duplicate CA identity')
                data[i]=(np.asarray([float(x),float(y),float(z)]),float(b))
        if set(data)!={x.id[1] for x in residues} or any(not np.allclose(data[x.id[1]][0],x['CA'].coord,rtol=1e-7,atol=1e-5) for x in residues):raise ValueError('Coordinate parser disagreement')
        if any(not np.isfinite(v[0]).all() or not math.isfinite(v[1]) or not 0<=v[1]<=100 for v in data.values()):raise ValueError('Invalid prediction data')
        cache[name]=data;return data
    seen=set();seen_links=set();geometries=0;pair_masks=0
    for entry in sorted({x['entry_id'] for x in links}):
        if entry not in mr['entry_receipt_sha256']:continue
        rp=a.mapping/(entry+'.receipt.json');er=json.loads(rp.read_text());path=a.mapping/(entry+'.residues.tsv.gz')
        if sha(rp)!=mr['entry_receipt_sha256'][entry] or er['config_sha256']!=mr['config_sha256'] or sha(path)!=er['table_sha256']:raise ValueError('Mapping changed')
        groups=defaultdict(list)
        with gzip.open(path,'rt') as stream:
            for x in csv.DictReader(stream,delimiter='\t'):
                if entry+'_'+x['entity_id'] in targets:groups[x['entity_id'],x['model_number'],x['label_asym_id']].append(x)
        for (entity,number,chain),records in groups.items():
            records.sort(key=lambda x:int(x['label_seq_id']))
            if [int(x['label_seq_id']) for x in records]!=list(range(1,len(records)+1)):raise ValueError('Incomplete experimental grid')
            observed={}
            for x in records:
                if x['CA_status']=='unambiguous_full_occupancy_CA':
                    atoms=json.loads(x['atom_records_json'])
                    if len(atoms)!=1:raise ValueError('Ambiguous experimental CA')
                    observed[int(x['label_seq_id'])]=np.asarray([float(atoms[0][k]) for k in ['Cartn_x','Cartn_y','Cartn_z']])
            for link in targets[entry+'_'+entity]:
                digest=hashlib.sha256(''.join(x['target_aa'] for x in records).encode()).hexdigest()
                if digest!=link['sequence_sha256']:raise ValueError('Experimental sequence differs')
                seen_links.add((link['entity_id'],link['alphafold_model_id'],link['esmfold_model_id']));af=model(link,'alphafold');esm=model(link,'esmfold')
                for cutoff in [0,70,90]:
                    key=(entry,entity,number,chain,link['alphafold_model_id'],link['esmfold_model_id'],str(cutoff));accepted,row=stored[key]
                    if key in seen:raise ValueError('Duplicate source grid')
                    seen.add(key);positions=[i for i in sorted(observed) if af[i][1]>=cutoff and esm[i][1]>=cutoff];n=len(positions)
                    if json.loads(row['residue_positions_json'])!=positions or int(row['matched_residues'])!=n or int(row['experimental_unambiguous_CA'])!=len(observed) or int(row['sequence_length'])!=len(records) or row['sequence_sha256']!=digest:raise ValueError('Residue mask/count differs')
                    close(row['fraction_full_sequence'],n/len(records))
                    for field,source in [('experimental_method','methods'),('resolution_combined_json','resolution_combined_json'),('initial_release_date','initial_release_date')]:
                        if row[field]!=refs[entry][source]:raise ValueError('Reference metadata changed')
                    if accepted!=(n>=50 and n>=len(records)/2):raise ValueError('Eligibility mismatch')
                    if not accepted:
                        if row['reason']!='fewer_than_50_CA_or_below_half_full_sequence':raise ValueError('Exclusion reason differs')
                        continue
                    xyz=[np.asarray([af[i][0] for i in positions]),np.asarray([esm[i][0] for i in positions]),np.asarray([observed[i] for i in positions])];distances=[cdist(x,x) for x in xyz];p=np.asarray(positions)
                    mask=np.triu((np.minimum.reduce(distances)<=15)&(np.abs(p[:,None]-p[None,:])>=3),1);pairs=int(mask.sum())
                    if pairs!=int(row['common_local_distance_pairs']):raise ValueError('Common distance-pair mask differs')
                    pair_masks+=1
                    for label,i,j in [('af_experiment',0,2),('esm_experiment',1,2),('af_esm',0,1)]:
                        sup=SVDSuperimposer();sup.set(xyz[j],xyz[i]);sup.run();close(row[label+'_rmsd_angstrom'],sup.get_rms());geometries+=1
                        delta=(distances[i]-distances[j])[mask]
                        for suffix,value in [('local_mean_absolute_difference_angstrom',float(np.mean(abs(delta))) if pairs else None),('local_rms_difference_angstrom',float(np.sqrt(np.mean(delta**2))) if pairs else None)]:
                            if value is None:
                                if row[label+'_'+suffix]!='':raise ValueError('Unexpected empty-pair statistic')
                            else:close(row[label+'_'+suffix],value)
        print(entry,'all control masks and geometries verified',flush=True)
    if seen!=set(stored) or geometries!=r['independent_superposition_checks']:raise ValueError('Incomplete comparison audit')
    missing={(x['entity_id'],x['alphafold_model_id'],x['esmfold_model_id']) for x in links}-seen_links
    p=a.comparisons/'unmapped_links.tsv';recorded=list(rows(p)) if p.exists() else []
    if len(recorded)!=len(missing) or {(x['entity_id'],x['alphafold_model_id'],x['esmfold_model_id']) for x in recorded}!=missing:raise ValueError('Missing-link disposition differs')
    result={'status':'passed_full_experimental_control_geometry_readback','accepted_rows':counts[True],'excluded_rows':counts[False],'residue_masks_checked':len(seen),'common_pair_masks_checked':pair_masks,'independent_geometry_triplets':geometries//3,'independent_superpositions':geometries,'links_without_chain_grids':len(missing),'comparison_receipt_sha256':sha(a.comparisons/'receipt.json'),'script_sha256':sha(Path(__file__)),'scope':'Every source chain/model/threshold grid, residue mask, eligibility outcome, common pair count, RMSD and local-distance statistic independently recomputed using full cdist matrices and Bio.SVDSuperimposer. Structure-parser sequence/numbering checks plus float64 mmCIF fields avoid float32 cutoff changes. Protein summaries and biological independence remain separate.'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
