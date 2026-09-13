#!/usr/bin/env python3
"""Compare completed, independently audited ESMFold controls with exact AF references."""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import numpy as np
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1
from assess_pae_sensitivity import checked_receipt, confidence_mask
from compare_marker_structures import ROOT, sha, geometry
from retrieve_matched_models import polymer_sequences
from retrieve_marker_pae import validate_pae
from prepare_paired_phylogenetic_inputs import write_table


def compare_arrays(a, b, pae_a, pae_b, cutoff):
    if a.shape != b.shape or a.ndim != 2 or a.shape[1] != 4 or len(a) == 0 or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('Invalid matched coordinates/confidence')
    if (a[:,3]<0).any() or (a[:,3]>100).any() or (b[:,3]<0).any() or (b[:,3]>100).any():
        raise ValueError('Confidence outside 0..100')
    for pae in (pae_a,pae_b):
        if pae.shape!=(len(a),len(a)) or not np.isfinite(pae).all() or (pae<0).any():raise ValueError('Invalid PAE')
    mask=(a[:,3]>=cutoff)&(b[:,3]>=cutoff)
    positions=np.flatnonzero(mask)+1
    row={'plddt_cutoff':cutoff,'protein_length':len(a),'matched_residues':len(positions),
         'retained_fraction':len(positions)/len(a),
         'status':'compared' if len(positions)>=50 and len(positions)>=.5*len(a) else 'insufficient_coverage'}
    fields=['ca_superposition_rmsd_angstrom','local_distance_pairs','local_distance_mean_absolute_change_angstrom','local_distance_rms_change_angstrom']
    for p in (5,10,15):fields += [f'pae{p}_local_pairs',f'pae{p}_local_retained_fraction',f'pae{p}_local_mean_absolute_change_angstrom']
    row.update(dict.fromkeys(fields,''))
    if row['status']!='compared':return row
    x,y=a[mask,:3],b[mask,:3];row.update(geometry(x,y,positions,positions))
    dx=np.linalg.norm(x[:,None]-x[None,:],axis=2);dy=np.linalg.norm(y[:,None]-y[None,:],axis=2)
    local=np.triu(np.ones(dx.shape,dtype=bool),1)&((dx<=15)|(dy<=15))&(np.abs(positions[:,None]-positions[None,:])>=3)
    if int(local.sum())!=row['local_distance_pairs']:raise ValueError('Local distance accounting differs')
    for cutoff_pae in (5,10,15):
        retained=local&confidence_mask(pae_a,pae_b,positions,positions,cutoff_pae)
        row.update({f'pae{cutoff_pae}_local_pairs':int(retained.sum()),
            f'pae{cutoff_pae}_local_retained_fraction':float(retained.sum()/local.sum()) if local.any() else '',
            f'pae{cutoff_pae}_local_mean_absolute_change_angstrom':float(np.abs(dx[retained]-dy[retained]).mean()) if retained.any() else ''})
    return row


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['inputs','predictions','audit','pae','output']:p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():raise FileExistsError('Use a new immutable comparison output')
    inputs=checked_receipt(args.inputs);audit=checked_receipt(args.audit);pae_receipt=checked_receipt(args.pae)
    chunk=json.loads((args.predictions/'last_chunk.json').read_text());config=json.loads((args.predictions/'config.json').read_text())
    if chunk['interrupted'] or chunk['oom_deferred'] or chunk['remaining_eligible'] or chunk['length_or_alphabet_deferred']:raise ValueError('Complete control predictions required')
    if config['input_receipt_sha256']!=sha(args.inputs/'receipt.json') or audit['config_sha256']!=sha(args.predictions/'config.json') or audit['chunk_receipt_sha256']!=sha(args.predictions/'last_chunk.json'):raise ValueError('Prediction/audit provenance differs')
    if audit['status']!='complete_artifact_readback' or audit['predictions']!=inputs['prediction_candidates'] or audit['links_sha256']!=sha(args.inputs/'reference_links.tsv'):raise ValueError('Incomplete control readback')
    if pae_receipt['mapping_receipt_sha256']!=inputs['source_mapping_receipt_sha256']:raise ValueError('Reference PAE mapping differs')
    sequences={r.id:str(r.seq) for r in SeqIO.parse(args.inputs/'candidates.faa','fasta')}
    models=json.loads((args.inputs/'reference_models.json').read_text())
    paes={(r['model_id'],str(r['version'])):r for r in json.loads((args.pae/'pae_manifest.json').read_text())}
    if {'S'+m['sequence_sha256'] for m in models}!=set(sequences):raise ValueError('Reference universe differs')
    parser=PDBParser(QUIET=True);rows=[];provenance=[]
    for m in models:
        sid='S'+m['sequence_sha256'];sequence=sequences[sid]
        pr=args.predictions/(sid+'.json');pred=json.loads(pr.read_text())
        if pred['status']!='verified_prediction' or pred['config_sha256']!=sha(args.predictions/'config.json') or pred['sequence_sha256']!=m['sequence_sha256']:raise ValueError('Prediction identity differs')
        for name,digest in pred['artifacts'].items():
            if sha(args.predictions/name)!=digest:raise ValueError('Changed prediction')
        path=ROOT/m['path']
        if sha(path)!=m['sha256']:raise ValueError('Changed reference')
        cif=MMCIF2Dict(io.StringIO(path.read_text()))
        if polymer_sequences(cif)!=[sequence]:raise ValueError('Reference polymer sequence differs')
        ca={}
        for atom,res,x,y,z,c in zip(*[cif['_atom_site.'+k] for k in ['label_atom_id','label_seq_id','Cartn_x','Cartn_y','Cartn_z','B_iso_or_equiv']]):
            if atom=='CA':
                if int(res) in ca:raise ValueError('Duplicate reference CA')
                ca[int(res)]=[float(x),float(y),float(z),float(c)]
        if set(ca)!=set(range(1,len(sequence)+1)):raise ValueError('Reference CA sequence positions differ')
        af=np.array([ca[i] for i in sorted(ca)])
        structure=parser.get_structure(sid,args.predictions/(sid+'.pdb'));residues=list(structure.get_residues())
        if len(list(structure.get_models()))!=1 or len(list(structure.get_chains()))!=1 or ''.join(seq1(r.resname) for r in residues)!=sequence or [r.id[1] for r in residues]!=list(range(1,len(sequence)+1)):raise ValueError('ESMFold residue identity differs')
        with np.load(args.predictions/(sid+'.npz'),allow_pickle=False) as d:
            esm_pae=d['pae'].copy();confidence=d['ca_plddt'].copy()
            if str(d['sequence'])!=sequence or not np.allclose([r['CA'].bfactor for r in residues],confidence,atol=.0051,rtol=0):raise ValueError('ESMFold confidence correspondence differs')
        esm=np.column_stack(([r['CA'].coord for r in residues],confidence))
        pa=paes[(m['model_id'],str(m['version']))];pa_path=ROOT/pa['path']
        if pa['status']!='verified' or pa['sequence_sha256']!=m['sequence_sha256'] or sha(pa_path)!=pa['gzip_sha256']:raise ValueError('Reference PAE differs')
        af_pae=validate_pae(gzip.decompress(pa_path.read_bytes()),len(sequence))
        for cutoff in (0,70,90):rows.append({'sequence_id':sid,'af_model_id':m['model_id'],'af_version':m['version'],**compare_arrays(af,esm,af_pae,esm_pae,cutoff)})
        provenance.append({'sequence_id':sid,'af_model_id':m['model_id'],'af_version':m['version'],'prediction_receipt_sha256':sha(pr),'af_coordinate_sha256':m['sha256'],'af_pae_gzip_sha256':pa['gzip_sha256']})
    args.output.mkdir(parents=True);write_table(args.output/'comparisons.tsv',rows);write_table(args.output/'model_provenance.tsv',provenance)
    result={'status':'complete_audited_esmfold_af_reference_comparison','reference_models':len(models),'unique_sequences':len(sequences),'comparison_rows':len(rows),'script_sha256':sha(Path(__file__)),
            'source_receipts':{n:sha(getattr(args,n)/'receipt.json') for n in ['inputs','audit','pae']},
            'interpretation':'Exact-sequence predictor differences, conditional on jointly selected pLDDT sites and optional bidirectional PAE filters. Direct CA geometry, not additive branch length or experimental accuracy. Local distance change is a custom metric, not lDDT. Selection and shared sequence/training biases prevent causal or independent-replicate interpretations.',
            'artifacts':{x.name:sha(x) for x in args.output.iterdir() if x.is_file()}}
    (args.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
