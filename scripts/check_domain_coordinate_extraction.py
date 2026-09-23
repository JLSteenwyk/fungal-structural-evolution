#!/usr/bin/env python3
"""Verify PDB boundary, identity, rounding and confidence behavior."""
import hashlib
import io
import json
from pathlib import Path
import tarfile
import tempfile
from Bio.PDB import PDBParser
from extract_domain_coordinates import load_atoms,pdb_bytes,shard


def main():
    model=json.loads(next(Path('results/structures/whole-proteome-afdb-catalog-20260922-v1/models.jsonl').open()))
    sequence,residues,ca=load_atoms(model)
    start,end=2,min(40,len(sequence));data,missing=pdb_bytes(sequence,residues,start,end)
    structure=PDBParser(QUIET=True).get_structure('test',io.StringIO(data.decode()))
    exported=list(structure.get_residues());assert len(exported)==end-start+1
    for offset,residue in enumerate(exported):
        assert residue.id[1]==offset+1
        original=start+offset
        atoms={a[0]:a for a in residues[original]}
        for atom in residue:
            src=atoms[atom.name]
            assert max(abs(float(x)-float(y)) for x,y in zip(atom.coord,src[2:5]))<0.00051
            assert abs(atom.bfactor-src[6])<=0.00501
    assert all(1<=x<=len(sequence) for x in missing)
    bad=dict(model,sha256='0'*64)
    try:load_atoms(bad)
    except ValueError:pass
    else:raise AssertionError('Changed coordinate bytes accepted')
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp);key=Path(model['path']).stem;job=root/'job.jsonl';job.write_text(json.dumps({'model_key':key,'model':model,'intervals':[['test',start,end]]})+'\n')
        output=root/'output';output.mkdir();result=shard(job,output,0);assert result['counts']['exported']==1
        with tarfile.open(output/'job.tar') as archive:
            assert archive.getnames()==['test.pdb'] and archive.extractfile('test.pdb').read()==data
        assert shard(job,output,0)==result
    print('Passed independent PDB parser readback, source-hash rejection, tar member and completed-shard checks')


if __name__=='__main__':main()
