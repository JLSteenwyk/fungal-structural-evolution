#!/usr/bin/env python3
"""Exercise coordinate validation and checkpoint rejection on an isolated source copy."""
import argparse,json,tempfile
from pathlib import Path
from validate_duplication_coordinates import validate_model,run_shard

p=argparse.ArgumentParser();p.add_argument('--models',type=Path,required=True);a=p.parse_args()
with a.models.open() as f:m=json.loads(next(f))
with tempfile.TemporaryDirectory() as td:
    root=Path(td);source=root/'source.cif';source.write_bytes(Path(m['path']).read_bytes());m=dict(m,path=str(source))
    r=validate_model(m);assert len(r['ca_xyz'])==len(r['ca_plddt'])==r['length']==m['length']
    result=run_shard((0,[m],str(root),'fixture',0));assert result['counts']=={'validated':1}
    assert run_shard((0,[m],str(root),'fixture',0))==result
    bad=dict(m,mean_ca_plddt=m['mean_ca_plddt']+1)
    result=run_shard((1,[bad],str(root),'fixture',0));assert result['counts']=={'rejected_content':1}
    source.write_bytes(source.read_bytes()+b'\n')
    for index in [0,2]:
        try:run_shard((index,[m],str(root),'fixture',0))
        except ValueError:pass
        else:raise AssertionError('Changed source accepted')
print('Passed coordinate array/count check, checkpoint reuse, content rejection, changed-source checkpoint/new-run rejection.')
