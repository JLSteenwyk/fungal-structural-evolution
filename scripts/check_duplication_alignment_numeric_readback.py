#!/usr/bin/env python3
"""Native structural readback fixtures and corrupted-metric/mapping rejection."""
import copy
import json
import math
import subprocess
import tempfile
from pathlib import Path
from duplication_alignment_inputs import render_ca
from duplication_alignment_numeric_readback import load_pdb,check_alignment
from run_cross_clan_alignments import parse_output
from run_ortholog_pair_guide_comparison import sha

native='/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/US-align/USalign'
seq='ACDEFGHIKLMNPQRSTVWY'*2
with tempfile.TemporaryDirectory() as directory:
    root=Path(directory)
    for threshold in [None,70]:
        inputs=[]
        for side in [0,1]:
            xyz=[]
            for i in range(len(seq)):
                x,y,z=3*math.cos(i),3*math.sin(i),float(i)
                xyz.append([x,y,z] if side==0 else [-y+5,x-2,z+7+0.1*math.sin(i)])
            blob,sequence,positions=render_ca(dict(status='validated',sequence=seq,ca_xyz=xyz,
                               ca_plddt=[90 if i%3 else 20 for i in range(len(seq))]),threshold)
            path=root/(str(side)+'.pdb');path.write_bytes(blob)
            inputs.append(dict(path=str(path),sha256=sha(path),sequence=sequence,
                               original_positions=positions,retained_residues=len(sequence)))
        for a,b in [inputs,inputs[::-1]]:
            raw=subprocess.check_output([native,a['path'],b['path'],'-mol','prot','-mm','0','-outfmt','0','-ter','2'],text=True)
            r=dict(stdout=raw,metrics=parse_output(raw,[a['sequence'],b['sequence']]))
            left,right=load_pdb(a),load_pdb(b)
            result=check_alignment(r,left,right)
            assert result['rmsd_recomputed']>0 and result['rmsd_rounding_error']<=.00501
            for field,value in [('rmsd',9),('tm_left',.123),('length_left',999),('alignment_marks',' '*len(r['metrics']['alignment_marks']))]:
                bad=copy.deepcopy(r);bad['metrics'][field]=value
                try:check_alignment(bad,left,right)
                except ValueError:pass
                else:raise AssertionError('Accepted altered '+field)
            bad=copy.deepcopy(a);bad['original_positions'][0]+=1
            try:load_pdb(bad)
            except ValueError:pass
            else:raise AssertionError('Accepted altered residue mapping')
print('Passed full/masked both-order deformed native alignments, independent RMSD/identity and metric/mapping corruption rejection.')
