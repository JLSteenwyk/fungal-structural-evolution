#!/usr/bin/env python3
"""Exercise native tar input and verify all fixture AA and CA values."""
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile
import tempfile
import numpy as np
from build_whole_proteome_foldseek_database import readback,index
from catalog_whole_proteome_structures import sha

FOLDSEEK='/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/foldseek/bin/foldseek'


def main():
    root=Path('results/domains/domain-coordinates-20260923-v1/shards');source=root/'shard_00000.tar'
    records={}
    with (root/'shard_00000.jsonl').open() as f:
        for line in f:
            row=json.loads(line)
            if row['status']=='exported':records[row['member']]=row
    with tempfile.TemporaryDirectory() as temporary:
        tmp=Path(temporary);archive=tmp/'fixture.tar';models=[];coordinates={}
        with tarfile.open(source,'r') as original,tarfile.open(archive,'w') as output:
            for member in original.getmembers()[:2]:
                data=original.extractfile(member).read();r=records[member.name]
                assert hashlib.sha256(data).hexdigest()==r['pdb_sha256'];output.addfile(member,io.BytesIO(data))
                models.append({'path':member.name,'length':r['residues'],'sequence_sha256':r['fragment_sequence_sha256']})
                coordinates[Path(member.name).stem]=np.array([[float(line[a:b]) for a,b in [(30,38),(38,46),(46,54)]] for line in data.decode().splitlines() if line.startswith('ATOM  ') and line[12:16].strip()=='CA'],dtype=np.float32)
        paths=tmp/'archives.tsv';paths.write_text(str(archive)+'\n');prefix=tmp/'domains'
        command=[FOLDSEEK,'createdb',str(paths),str(prefix),'--threads','1','--gpu','0','--mask-bfactor-threshold','70','--coord-store-mode','1']
        subprocess.run(command,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE,check=True)
        checked=readback(prefix,models);idx=index(Path(str(prefix)+'_ca.index'))
        with Path(str(prefix)+'_ca').open('rb') as handle,Path(str(prefix)+'.lookup').open() as lookup:
            for line in lookup:
                key,name,_=line.rstrip('\n').split('\t');offset,length=idx[int(key)];handle.seek(offset);raw=handle.read(length)
                observed=np.frombuffer(raw[:-1],dtype=np.float32).reshape(3,-1).T
                if not np.array_equal(observed,coordinates[name]):raise ValueError('Foldseek C-alpha coordinates differ from exported PDB')
    result={'status':'passed_native_domain_tar_input_sequence_and_coordinate_fixture','models':checked['models'],'residues':checked['residues'],'foldseek_sha256':sha(FOLDSEEK),'script_sha256':sha(__file__),'scope':'Two serialized domain records from a completed source-atom-checked archive; verifies TSV-of-tar input, lookup names, exact AA hashes, 3Di dimensions/alphabet and every float32 C-alpha coordinate. Does not validate full domain database or 3Di feature accuracy.'}
    Path('metadata/domain_foldseek_archive_fixture_checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
