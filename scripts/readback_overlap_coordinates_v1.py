#!/usr/bin/env python3
"""Replay every raw model, original mask and geometry row without producer helpers."""
import argparse
from collections import Counter
import csv
from datetime import datetime,timezone
import hashlib
import io
import json
import math
from pathlib import Path

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Data.PDBData import protein_letters_3to1
import numpy as np
from scipy.spatial.distance import pdist
from scipy.spatial.transform import Rotation

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind,verify


METRICS=['ca_superposition_rmsd_angstrom','local_distance_pairs',
    'local_distance_mean_absolute_change_angstrom','local_distance_rms_change_angstrom',
    'all_distance_pairs','all_distance_mean_absolute_change_angstrom','all_distance_rms_change_angstrom',
    'sequence_local_distance_pairs','sequence_local_distance_mean_absolute_change_angstrom',
    'sequence_local_distance_rms_change_angstrom']


def raw_model(row,prefix):
    path=Path(row[prefix+'_model_path']);data=path.read_bytes()
    if hashlib.sha256(data).hexdigest()!=row[prefix+'_model_sha256']:raise ValueError('Coordinate checksum')
    # Shared standard mmCIF tokenization only; residue/atom/confidence validation is rebuilt here.
    cif=MMCIF2Dict(io.StringIO(data.decode()))
    seqs=cif.get('_entity_poly.pdbx_seq_one_letter_code_can',[])
    if len(seqs)!=1:raise ValueError('Polymer count')
    sequence=''.join(seqs[0].split());n=int(row[prefix+'_length'])
    if len(sequence)!=n or hashlib.sha256(sequence.encode()).hexdigest()!=row[prefix+'_sequence_sha256']:
        raise ValueError('Complete polymer identity')
    keys=['group_PDB','label_seq_id','label_comp_id','label_atom_id','label_asym_id','label_alt_id',
          'pdbx_PDB_model_num','Cartn_x','Cartn_y','Cartn_z','occupancy','B_iso_or_equiv','type_symbol']
    columns=[cif['_atom_site.'+k] for k in keys]
    if len({len(c) for c in columns})!=1:raise ValueError('Column lengths')
    atoms=set();chains=set();coords={};confidence={}
    for record in zip(*columns):
        group,number,residue,atom,chain,alt,model,x,y,z,occupancy,b,element=record
        if group!='ATOM' or alt not in ['.','?'] or model!='1':raise ValueError('Unsupported atom record')
        position=int(number)
        if not 1<=position<=n or protein_letters_3to1.get(residue)!=sequence[position-1]:
            raise ValueError('Atom sequence position')
        if (position,atom) in atoms:raise ValueError('Duplicate atom')
        atoms.add((position,atom));chains.add(chain)
        values=[float(v) for v in [x,y,z,occupancy,b]]
        if not all(math.isfinite(v) for v in values) or not 0<=values[3]<=1 or not 0<=values[4]<=100:
            raise ValueError('Invalid numeric atom record')
        if len(atom)>4 or len(element)>2:raise ValueError('Invalid atom identity')
        if atom=='CA':coords[position]=values[:3];confidence[position]=values[4]
    if len(chains)!=1 or sorted(coords)!=list(range(1,n+1)):raise ValueError('Complete single-chain CA coverage')
    ep=Path(row[prefix+'_encoding_path'])
    if sha(ep)!=row[prefix+'_encoding_sha256']:raise ValueError('Encoding checksum')
    with np.load(ep,allow_pickle=False) as f:d={k:f[k].copy() for k in f.files}
    if str(d['sequence'])!=sequence:raise ValueError('Encoded sequence identity')
    if not np.array_equal(d['ca_plddt'],np.array([confidence[i] for i in range(1,n+1)])):
        raise ValueError('Raw/encoded confidence disagreement')
    for key in ['valid','feature_min_plddt','feature_max_pae','partner_residue_1based']:
        if d[key].shape!=(n,):raise ValueError('Feature dimensions')
    if d['valid'].dtype!=bool or len(str(d['states']))!=n:raise ValueError('State/valid dimensions')
    for key in ['feature_min_plddt','feature_max_pae']:
        if not np.isfinite(d[key][d['valid']]).all():raise ValueError('Nonfinite valid feature confidence')
    d['letters']=list(str(d['states']))
    return d,np.array([coords[i] for i in range(1,n+1)],dtype=np.float64)


def mask_and_counts(a,b,cutoff,pae):
    positions=[];mismatches=0;changed=0;same=0;same_mismatch=0;changed_mismatch=0
    n=len(a['valid'])
    for i in range(n):
        keep=bool(a['valid'][i]) and bool(b['valid'][i])
        keep=keep and float(a['feature_min_plddt'][i])>=cutoff and float(b['feature_min_plddt'][i])>=cutoff
        if pae is not None:keep=keep and float(a['feature_max_pae'][i])<=pae and float(b['feature_max_pae'][i])<=pae
        if not keep:continue
        positions.append(i+1);mismatch=a['letters'][i]!=b['letters'][i]
        change=int(a['partner_residue_1based'][i])!=int(b['partner_residue_1based'][i])
        mismatches+=mismatch;changed+=change;same+=not change
        same_mismatch+=mismatch and not change;changed_mismatch+=mismatch and change
    counts=dict(protein_length=n,retained_residues=len(positions),state_mismatches=mismatches,
        partner_changes=changed,same_partner_residues=same,same_partner_state_mismatches=same_mismatch,
        changed_partner_residues=changed,changed_partner_state_mismatches=changed_mismatch)
    return np.array(positions,dtype=np.int64),counts


def geometry_replay(x,y,positions):
    before=(x.tobytes(),y.tobytes(),positions.tobytes())
    cx=x-x.mean(axis=0);cy=y-y.mean(axis=0)
    rotation,_=Rotation.align_vectors(cx,cy)
    result={'ca_superposition_rmsd_angstrom':float(np.sqrt(np.mean(np.sum((rotation.apply(cy)-cx)**2,axis=1))))}
    first=pdist(x);second=pdist(y);delta=first-second
    ii,jj=np.triu_indices(len(x),1);gap=positions[jj]-positions[ii]
    masks={'all':np.full(len(delta),True,dtype=bool),
        'local':((first<=15)|(second<=15))&(gap>=3),
        'sequence_local':(gap>=3)&(gap<=10)}
    for name,mask in masks.items():
        values=delta[mask];count=len(values)
        result[name+'_distance_pairs']=count
        result[name+'_distance_mean_absolute_change_angstrom']=math.fsum(abs(float(v)) for v in values)/count if count else ''
        result[name+'_distance_rms_change_angstrom']=math.sqrt(math.fsum(float(v)*float(v) for v in values)/count) if count else ''
    assert before==(x.tobytes(),y.tobytes(),positions.tobytes()),'Reader coordinate memory changed'
    assert set(result)==set(METRICS)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['producer','transport','receipt']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();assert not a.receipt.exists();producer=json.loads(a.producer.read_text());transport=json.loads(a.transport.read_text())
    assert producer['status']=='complete_full_matched_predictor_coordinate_benchmark_pending_independent_readback'
    assert transport['validation_sha256']==sha(a.producer) and transport['original_tool_terminal_exit_code']==0
    pins=dict(producer['source_hashes'])
    for q in [a.producer,a.transport,Path(__file__),Path(producer['table'])]:bind(pins,q)
    assert sha(producer['table'])==producer['table_sha256'];verify(pins)
    links=list(csv.DictReader(Path('results/phylogeny/paired-source-model-context-20260926-v1/source_model_context.tsv').open(),delimiter='\t'))
    pairs={tuple(r[k] for k in ['reference_model_id','reference_model_version','local_model_id','local_model_version']):r for r in links}
    assert len(links)==673 and len(pairs)==643
    table=list(csv.DictReader(Path(producer['table']).open(),delimiter='\t'));assert len(table)==7716
    original={}
    for r in csv.DictReader(Path('results/phylogeny/overlap-model-feature-comparison-20260926-v1/model_pair_comparisons.tsv').open(),delimiter='\t'):
        key=tuple(r[k] for k in ['reference_model_id','reference_version','local_model_id','local_version'])
        original[key,int(r['plddt_cutoff']),r['pae_cutoff']]=r
    seen=set();statuses=Counter();max_abs=0.;compared=0;load_errors=0
    for index,(key,link) in enumerate(sorted(pairs.items())):
        loaded=[];issue=False
        try:loaded=[raw_model(link,prefix) for prefix in ['reference','local']]
        except (ValueError,KeyError,AssertionError):issue=True;load_errors+=1
        for ci,cutoff in enumerate([0,70,90]):
            for pi,pae in enumerate([None,5,10,15]):
                r=table[index*12+ci*4+pi];conf=str(pae) if pae is not None else 'unfiltered'
                assert tuple(r[k] for k in ['reference_model_id','reference_version','local_model_id','local_version'])==key
                assert int(r['plddt_cutoff'])==cutoff and r['pae_cutoff']==conf
                tag=(key,cutoff,conf);assert tag not in seen;seen.add(tag)
                identifier=hashlib.sha256(json.dumps(key,separators=(',',':')).encode()).hexdigest()
                assert r['model_pair_id']==identifier and r['sequence_sha256']==link['reference_sequence_sha256']
                assert r['scientific_eligibility']=='False'
                old=original[tag]
                for name in ['protein_length','retained_residues','state_mismatches','partner_changes','same_partner_residues',
                    'same_partner_state_mismatches','changed_partner_residues','changed_partner_state_mismatches']:
                    assert r[name]==old[name]
                assert r['feature_coverage_status']==old['status']
                if issue:
                    assert r['status']=='coordinate_validation_rejected' and r['coordinate_validation_issue']
                    assert all(r[name]=='' for name in METRICS)
                else:
                    (ea,x),(eb,y)=loaded;positions,counts=mask_and_counts(ea,eb,cutoff,pae)
                    assert all(int(r[k])==v for k,v in counts.items())
                    enough=len(positions)>=50 and len(positions)*2>=len(x)
                    assert r['feature_coverage_status']==('compared' if enough else 'insufficient_coverage')
                    assert r['coordinate_validation_issue']==''
                    if not enough:
                        assert r['status']=='insufficient_common_coverage' and all(r[name]=='' for name in METRICS)
                    else:
                        assert r['status']=='coordinates_compared'
                        actual=geometry_replay(x[positions-1],y[positions-1],positions)
                        for name,value in actual.items():
                            if value=='':assert r[name]==''
                            elif name.endswith('_pairs'):assert int(r[name])==value
                            else:
                                saved=float(r[name]);assert math.isfinite(saved) and saved>=0
                                error=abs(saved-value);max_abs=max(max_abs,error)
                                assert math.isclose(saved,value,rel_tol=1e-10,abs_tol=1e-10),(tag,name,saved,value)
                        compared+=1
                statuses[r['status']]+=1
        print('independent_coordinate_pairs',index+1,'/643',flush=True)
    assert len(seen)==7716 and dict(statuses)==producer['status_counts'];verify(pins)
    result=dict(status='passed_full_independent_matched_predictor_coordinate_readback',checked_utc=datetime.now(timezone.utc).isoformat(),
        producer_receipt=str(a.producer),producer_receipt_sha256=sha(a.producer),model_pairs=643,raw_models_attempted=1286,
        confidence_rows=7716,geometry_rows_compared=compared,coordinate_validation_rejected_pairs=load_errors,
        status_counts=dict(statuses),maximum_geometry_absolute_difference=max_abs,geometry_relative_tolerance=1e-10,
        geometry_absolute_tolerance=1e-10,source_hashes=pins,scientific_eligibility=False,all_eight_aims_incomplete=True,
        scope='All original pairs, raw model sequence/all-atom/CA confidence contracts, scalar confidence masks, state/partner '
              'counts and proper-rotation/all-pair/spatial-contact/sequence-local geometry independently rebuilt. Standard '
              'Biopython mmCIF lexical parser is shared; producer atom/geometry/mask helpers are not used. Rotation.align_vectors, '
              'SciPy pdist and math.fsum supply independent geometric arithmetic. Predeclared geometry roundoff tolerance '
              'does not alter the earlier weighted numerical guard. Coverage/rejections explicitly retained. These conditional '
              'predictor differences do not establish experimental accuracy, additive distances, evolutionary acceleration, '
              'branch calibration, independent model pairs, posterior adequacy or biological acceptance.')
    with a.receipt.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
