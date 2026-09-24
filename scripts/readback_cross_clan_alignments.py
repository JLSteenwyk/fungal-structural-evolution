#!/usr/bin/env python3
"""Independently verify full alignment mappings, RMSD, identity and confidence."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from Bio.Data.PDBData import protein_letters_3to1


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()


def rmsd(x,y):
    a=x-x.mean(axis=0);b=y-y.mean(axis=0)
    u,_,vt=np.linalg.svd(a.T@b)
    correction=np.eye(3);correction[-1,-1]=1 if np.linalg.det(u@vt)>=0 else -1
    return float(np.sqrt(np.mean(np.sum((a@u@correction@vt-b)**2,axis=1))))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        if sha(args.plan)!=ph:raise ValueError('Plan changed')
        for p,h in plan['pins'].items():
            if sha(p)!=h:raise ValueError('Source changed: '+p)
    verify();root=Path(plan['alignments']);coords=Path(plan['coordinates'])
    receipt=json.loads((root/'receipt.json').read_text());cp=json.loads((coords/'receipt.json').read_text())
    producer=json.loads(Path(plan['producer_plan']).read_text())
    if receipt['status']!='complete_cross_clan_alignments_pending_independent_readback':raise ValueError('Completed alignment run required')
    if receipt['plan_sha256']!=plan['pins'][plan['producer_plan']]:raise ValueError('Producer plan differs')
    for p,h in producer['pins'].items():
        if sha(p)!=h:raise ValueError('Producer input changed')
    models={}
    for name,h in cp['pdb_hashes'].items():
        p=coords/'pdb'/name
        if sha(p)!=h:raise ValueError('Coordinate checksum differs')
        positions=[];sequence=[];xyz=[];confidence=[]
        for line in p.read_text().splitlines():
            if line.startswith('ATOM  ') and line[12:16].strip()=='CA':
                positions.append(int(line[22:26]));sequence.append(protein_letters_3to1[line[17:20]])
                xyz.append([float(line[a:b]) for a,b in [(30,38),(38,46),(46,54)]])
                confidence.append(float(line[60:66]))
        if positions!=list(range(1,len(positions)+1)):raise ValueError('Invalid residue grid')
        models[p.stem]=(''.join(sequence),np.asarray(xyz),np.asarray(confidence))
    # Relational reconstruction is separate from the producer's nested set loops.
    members=pd.read_csv(producer['members'],sep='\t');members=members[members.model_exclusive_within_pair_cluster==1]
    keys=['boundary','representative','pfam_left','pfam_right']
    joined=members[members.side=='left'].merge(members[members.side=='right'],on=keys,suffixes=('_a','_b'))
    expected=set(zip(joined.interval_id_a,joined.interval_id_b))|set(zip(joined.interval_id_b,joined.interval_id_a))
    if len(expected)!=receipt['directed_alignments']:raise ValueError('Directed scope differs')
    seen=set();out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    fields=['interval_left','interval_right','length_left','length_right','aligned_length',
            'sequence_identity_exact','rmsd_recomputed','rmsd_native','tm_left_native','tm_right_native',
            'coverage_left','coverage_right','joint_plddt70_pairs','joint_plddt70_fraction']
    max_error=0.;table=out/'alignment_readback.tsv'
    with table.open('w') as handle:
        writer=csv.writer(handle,delimiter='\t',lineterminator='\n');writer.writerow(fields)
        for name,h in sorted(receipt['artifacts'].items()):
            p=root/name
            if sha(p)!=h:raise ValueError('Alignment checkpoint changed')
            record=json.loads(p.read_text());a,b=record['intervals'];key=(a,b)
            if key not in expected or key in seen or record['plan_sha256']!=receipt['plan_sha256']:raise ValueError('Wrong pair identity')
            seen.add(key);m=record['metrics'];lines=record['stdout'].splitlines()
            summary=next(l for l in lines if l.startswith('Aligned length='))
            segments=summary.split(',');n=int(segments[0].split('=')[-1]);native_rmsd=float(segments[1].split('=')[-1]);native_identity=float(segments[2].split('=')[-1])
            scorelines=[l for l in lines if l.startswith('TM-score=')]
            if len(scorelines)!=2:raise ValueError('Unexpected score output')
            scores=[float(l.split()[1]) for l in scorelines]
            marker=next(i for i,l in enumerate(lines) if 'denotes residue pairs of' in l)
            x,marks,y=lines[marker+1:marker+4]
            if (x,marks,y)!=(m['alignment_left'],m['alignment_marks'],m['alignment_right']):raise ValueError('Serialized alignment differs')
            sx,cx,px=models[a];sy,cy,py=models[b]
            if x.replace('-','')!=sx or y.replace('-','')!=sy or len(x)!=len(y) or len(marks)!=len(x):raise ValueError('Sequence mapping differs')
            ix=[];iy=[];i=j=0;identical=0
            for left,right,mark in zip(x,y,marks):
                if left!='-' and right!='-':
                    if mark not in ':.':raise ValueError('Missing matched-residue marker')
                    ix.append(i);iy.append(j);identical+=left==right
                elif mark!=' ' or left==right=='-':raise ValueError('Invalid gap marker')
                i+=left!='-';j+=right!='-'
            if len(ix)!=n or n!=m['aligned_length']:raise ValueError('Aligned residue count differs')
            calculated=rmsd(cx[ix],cy[iy]);error=abs(calculated-native_rmsd);max_error=max(max_error,error)
            if error>.00501:raise ValueError('RMSD differs beyond printed rounding: '+str((key,error)))
            identity=identical/n
            if abs(identity-native_identity)>.000501:raise ValueError('Identity differs beyond rounding')
            if native_rmsd!=m['rmsd'] or native_identity!=m['sequence_identity'] or scores!=[m['tm_left'],m['tm_right']]:raise ValueError('Stored metrics differ from stdout')
            if [len(sx),len(sy)]!=[m['length_left'],m['length_right']]:raise ValueError('Stored lengths differ')
            high=int(np.count_nonzero((px[ix]>=70)&(py[iy]>=70)))
            writer.writerow([a,b,len(sx),len(sy),n,identity,calculated,native_rmsd,*scores,n/len(sx),n/len(sy),high,high/n])
    if seen!=expected:raise ValueError('Incomplete pair grid')
    verify()
    result=dict(status='passed_full_cross_clan_mapping_rmsd_identity_readback',plan_sha256=ph,
                directed_alignments=len(seen),unordered_pairs=len(seen)//2,max_rmsd_rounding_error=max_error,
                artifacts={table.name:sha(table)},
                scope='All source coordinate hashes, directed pair identities, residue correspondences, aligned lengths, '
                'sequence identities and least-squares RMSDs verified independently of producer parser. '
                'TM-scores checked against raw native text only, not independently reoptimized. '
                'Matched confidence uses rounded PDB C-alpha pLDDT; no PAE qualification or homology conclusions.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
