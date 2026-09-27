"""Independent native-text, residue-map, identity and least-squares RMSD checks."""
import math
from pathlib import Path
import numpy as np
from Bio.Data.PDBData import protein_letters_3to1
from readback_cross_clan_alignments import rmsd
from run_ortholog_pair_guide_comparison import sha


def load_pdb(row):
    path=Path(row['path'])
    if sha(path)!=row['sha256']:raise ValueError('Changed PDB bytes')
    positions=[];letters=[];xyz=[];confidence=[]
    for line in path.read_text().splitlines():
        if not line.startswith('ATOM  '):continue
        if line[12:16].strip()!='CA' or line[21]!='A':raise ValueError('Unexpected atom or chain')
        positions.append(int(line[22:26]));letters.append(protein_letters_3to1[line[17:20]])
        xyz.append([float(line[a:b]) for a,b in [(30,38),(38,46),(46,54)]])
        confidence.append(float(line[60:66]))
    if positions!=row['original_positions'] or ''.join(letters)!=row['sequence'] or len(positions)!=row['retained_residues']:
        raise ValueError('PDB residue mapping differs')
    if len(positions)!=len(set(positions)) or positions!=sorted(positions):raise ValueError('Invalid residue order')
    if not np.isfinite(xyz).all() or not all(math.isfinite(x) and 0<=x<=100 for x in confidence):raise ValueError('Invalid PDB values')
    return ''.join(letters),np.asarray(xyz),np.asarray(confidence)


def check_alignment(record,left,right):
    m=record['metrics'];lines=record['stdout'].splitlines()
    summaries=[l for l in lines if l.startswith('Aligned length=')]
    if len(summaries)!=1:raise ValueError('Ambiguous summary')
    segments=summaries[0].split(',')
    n=int(segments[0].split('=')[-1]);native_rmsd=float(segments[1].split('=')[-1]);native_identity=float(segments[2].split('=')[-1])
    scorelines=[l for l in lines if l.startswith('TM-score=')]
    if len(scorelines)!=2:raise ValueError('Unexpected TM-score count')
    scores=[float(l.split()[1]) for l in scorelines]
    for index,line in enumerate(scorelines,1):
        if 'normalized by length of Structure_'+str(index)+':' not in line:raise ValueError('TM-score normalization differs')
    markers=[i for i,l in enumerate(lines) if 'denotes residue pairs of' in l]
    if len(markers)!=1:raise ValueError('Ambiguous alignment block')
    x,marks,y=lines[markers[0]+1:markers[0]+4]
    if (x,marks,y)!=(m['alignment_left'],m['alignment_marks'],m['alignment_right']):raise ValueError('Stored alignment differs')
    sx,cx,px=left;sy,cy,py=right
    if x.replace('-','')!=sx or y.replace('-','')!=sy or len(x)!=len(y) or len(marks)!=len(x):raise ValueError('Alignment sequence differs')
    lengths=[]
    for index in [1,2]:
        match=[l for l in lines if l.startswith('Length of Structure_'+str(index)+':')]
        if len(match)!=1:raise ValueError('Ambiguous native length')
        lengths.append(int(match[0].split(':',1)[1].split()[0]))
    if lengths!=[len(sx),len(sy)] or lengths!=[m['length_left'],m['length_right']]:raise ValueError('Native/stored lengths differ')
    ix=[];iy=[];i=j=identical=0
    for a,b,mark in zip(x,y,marks):
        if a!='-' and b!='-':
            if mark not in ':.':raise ValueError('Missing paired-residue marker')
            ix.append(i);iy.append(j);identical+=a==b
        elif mark!=' ' or a==b=='-':raise ValueError('Invalid gap marker')
        i+=a!='-';j+=b!='-'
    if len(ix)!=n or n!=m['aligned_length'] or n<1:raise ValueError('Aligned length differs')
    if not all(math.isfinite(v) for v in [native_rmsd,native_identity,*scores]) or native_rmsd<0 or not all(0<=v<=1 for v in [native_identity,*scores]):raise ValueError('Invalid metric value')
    calculated=rmsd(cx[ix],cy[iy]);error=abs(calculated-native_rmsd);identity=identical/n
    if error>.00501:raise ValueError('RMSD exceeds printed rounding')
    if abs(identity-native_identity)>.000501:raise ValueError('Identity exceeds printed rounding')
    if native_rmsd!=m['rmsd'] or native_identity!=m['sequence_identity'] or scores!=[m['tm_left'],m['tm_right']]:raise ValueError('Stored numeric metrics differ')
    high=int(np.count_nonzero((px[ix]>=70)&(py[iy]>=70)))
    return dict(aligned_length=n,rmsd_recomputed=calculated,rmsd_native=native_rmsd,
                rmsd_rounding_error=error,sequence_identity_exact=identity,
                tm_left_native=scores[0],tm_right_native=scores[1],
                coverage_left=n/len(sx),coverage_right=n/len(sy),
                joint_plddt70_pairs=high,joint_plddt70_fraction=high/n)
