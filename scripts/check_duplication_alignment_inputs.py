#!/usr/bin/env python3
"""Verify full/masked residue mapping and native US-align behavior on a rigid transform."""
import copy,json,math,subprocess,tempfile
from pathlib import Path
from Bio.PDB import PDBParser
from duplication_alignment_inputs import render_ca
from run_cross_clan_alignments import parse_output

sequence='ACDEFGHIKLMNPQRSTVWY'
xyz=[[3*math.cos(i),3*math.sin(i),1.5*i] for i in range(len(sequence))]
record={'status':'validated','sequence':sequence,'ca_xyz':xyz,'ca_plddt':[70 if i%3 else 69.9 for i in range(len(sequence))]}
with tempfile.TemporaryDirectory() as td:
    p=Path(td);checked=0
    for mask in [None,70]:
        blob,seq,positions=render_ca(record,mask);expected=[i+1 for i,v in enumerate(record['ca_plddt']) if mask is None or v>=mask];assert positions==expected
        assert seq==''.join(sequence[i-1] for i in positions)
        a=p/'a.pdb';a.write_bytes(blob);structure=PDBParser(QUIET=True).get_structure('a',str(a));atoms=list(structure.get_atoms())
        assert [x.get_parent().id[1] for x in atoms]==positions and len(atoms)==len(seq)
        transformed=copy.deepcopy(record);transformed['ca_xyz']=[[-y+10,x-5,z+3] for x,y,z in xyz];blob2,seq2,pos2=render_ca(transformed,mask);assert seq2==seq and pos2==positions;b=p/'b.pdb';b.write_bytes(blob2)
        exe='/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/US-align/USalign'
        for left,right in [(a,b),(b,a)]:
            result=subprocess.run([exe,str(left),str(right),'-mol','prot','-mm','0','-outfmt','0','-ter','2'],check=True,capture_output=True,text=True)
            m=parse_output(result.stdout,[seq,seq]);assert m['rmsd']<.01 and m['tm_left']>.999 and m['tm_right']>.999 and m['aligned_length']==len(seq);checked+=1
    empty=copy.deepcopy(record);empty['ca_plddt']=[0]*len(sequence);assert render_ca(empty,70)[1:]==('',[])
    bad=copy.deepcopy(record);bad['ca_xyz'][0][0]=float('nan')
    try:render_ca(bad)
    except ValueError:pass
    else:raise AssertionError('Nonfinite coordinate accepted')
print(json.dumps({'native_directed_transform_checks':checked,'full_and_noncontiguous_mask_mapping':'passed','empty_mask_retained':'passed','nonfinite_rejected':'passed'}))
