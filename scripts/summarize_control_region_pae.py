#!/usr/bin/env python3
"""Summarize ESMFold PAE within and between all annotated control hit spans."""
import argparse
import itertools
import json
from collections import defaultdict
from pathlib import Path
import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from audit_joint_path_uncertainty import checked,rows,sha
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['controls','crosswalk','pfam','output']:
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    checked(a.controls);cr=checked(a.crosswalk);checked(a.pfam)
    if cr['source_receipts']['controls']!=sha(a.controls/'receipt.json'):raise ValueError('Control lineage mismatch')
    models=json.loads((a.controls/'model_provenance.json').read_text());links={}
    for row in rows(a.crosswalk/'crosswalk.tsv'):
        sid=row['sequence_sha256']
        if sid in links and links[sid]['alphafold_model_id']!=row['alphafold_model_id']:raise ValueError('Ambiguous paired AF model')
        links[sid]=row
    hits=defaultdict(list)
    for h in rows(a.pfam/'raw_annotated_hits.tsv'):
        if h['sequence_id'][1:] in links:hits[h['sequence_id'][1:]].append(h)
    root=Path(__file__).resolve().parents[1];output=[];npz_pins={}
    for model in models:
        sid=model['sequence_sha256'];link=links[sid];file=root/model['local_pae_npz_path']
        if sha(file)!=model['local_pae_npz_sha256']:raise ValueError('Changed PAE artifact')
        npz_pins[str(file)]=sha(file)
        with np.load(file,allow_pickle=False) as d:pae=d['pae'].copy();esm=d['ca_plddt'].copy()
        n=int(model['length'])
        if pae.shape!=(n,n) or esm.shape!=(n,) or not np.isfinite(pae).all():raise ValueError('Invalid PAE dimensions/values')
        afpath=root/link['alphafold_path']
        if sha(afpath)!=link['alphafold_sha256']:raise ValueError('Changed AF prediction')
        d=MMCIF2Dict(str(afpath));af={}
        for atom,pos,value in zip(d['_atom_site.label_atom_id'],d['_atom_site.label_seq_id'],d['_atom_site.B_iso_or_equiv']):
            if atom=='CA':
                if int(pos) in af:raise ValueError('Duplicate AF CA')
                af[int(pos)]=float(value)
        if sorted(af)!=list(range(1,n+1)):raise ValueError('Incomplete AF sequence positions')
        confidence=np.minimum(esm,np.array([af[i] for i in range(1,n+1)]))
        spans=[]
        for h in hits[sid]:
            lo,hi=int(h['alignment_start']),int(h['alignment_end'])
            if int(h['protein_length'])!=n or not 1<=lo<=hi<=n:raise ValueError('Annotation coordinates differ')
            spans.append((h,np.arange(lo-1,hi)))
        for cutoff in [0,70,90]:
            for (left,li),(right,ri) in itertools.combinations_with_replacement(spans,2):
                within=left['hit_id']==right['hit_id'];overlap=bool(np.intersect1d(li,ri).size) and not within
                i=li[confidence[li]>=cutoff];j=ri[confidence[ri]>=cutoff]
                status=('excluded_overlapping_hits' if overlap else 'excluded_coverage' if min(len(i),len(j))<20 or len(i)*2<len(li) or len(j)*2<len(ri) else 'included')
                row=dict(sequence_sha256=sid,alphafold_model_id=link['alphafold_model_id'],esmfold_model_id=model['model_id'],joint_predicted_plddt_cutoff=cutoff,kind='within_hit' if within else 'between_hits',left_hit_id=left['hit_id'],left_name=left['pfam_name'],right_hit_id=right['hit_id'],right_name=right['pfam_name'],left_residues=len(i),right_residues=len(j),status=status)
                for direction,first,second in [('row_left_column_right',i,j),('row_right_column_left',j,i)]:
                    values=pae[np.ix_(first,second)]
                    if within:values=values[~np.eye(len(first),dtype=bool)]
                    else:values=values.ravel()
                    row[direction+'_pairs']=len(values) if status=='included' else 0
                    for label,q in [('median',.5),('q90',.9)]:
                        row[direction+'_'+label+'_angstrom']=float(np.quantile(values,q)) if status=='included' else ''
                output.append(row)
    a.output.mkdir(parents=True);write_table(a.output/'region_pae.tsv',output)
    r=dict(status='complete_control_region_pae_summary',proteins=len(models),pfam_hits=sum(map(len,hits.values())),rows=len(output),included_rows=sum(x['status']=='included' for x in output),
        source_receipts={name:sha(getattr(a,name)/'receipt.json') for name in ['controls','crosswalk','pfam']},npz_sha256=npz_pins,script_sha256=sha(Path(__file__)),
        interpretation='ESMFold PAE, retaining both matrix directions; within-hit diagonal excluded. All full canonical positions passing joint AF/ESM focal confidence considered, without experimental-coverage restriction. Minimum 20 residues and half span per region. Overlapping different hits excluded explicitly. Raw Pfam hits are not independent structural domains. Predicted uncertainty is not measured error, biological motion or evidence of calibration.',
        artifacts={'region_pae.tsv':sha(a.output/'region_pae.tsv')})
    (a.output/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k!='npz_sha256'},indent=2))


if __name__=='__main__':main()
