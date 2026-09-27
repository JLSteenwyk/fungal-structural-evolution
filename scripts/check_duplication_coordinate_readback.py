#!/usr/bin/env python3
"""Test independent readback against altered exported values and source identities."""
import argparse,copy,json
from pathlib import Path
from validate_duplication_coordinates import validate_model
from readback_duplication_coordinates import check_record
p=argparse.ArgumentParser();p.add_argument('--models',type=Path,required=True);a=p.parse_args()
with a.models.open() as f:m=json.loads(next(f))
blob=Path(m['path']).read_bytes();r=dict(model_id=m['model_id'],version=m['version'],source_path=m['path'],source_sha256=m['sha256'],sequence_sha256=m['sequence_sha256'],**validate_model(m))
assert check_record(r,m,blob)==m['length']
changes=[lambda x:x['ca_xyz'][0].__setitem__(0,x['ca_xyz'][0][0]+.1),lambda x:x['ca_plddt'].__setitem__(0,x['ca_plddt'][0]+.1),lambda x:x.update(sequence='X'+x['sequence'][1:]),lambda x:x.update(residues_ge70=x['residues_ge70']+1),lambda x:x.update(model_id='wrong'),lambda x:x.update(mean_ca_plddt=x['mean_ca_plddt']+1)]
for mutate in changes:
    altered=copy.deepcopy(r);mutate(altered)
    try:check_record(altered,m,blob)
    except ValueError:pass
    else:raise AssertionError('Changed export accepted')
try:check_record(r,m,blob+b'\n')
except ValueError:pass
else:raise AssertionError('Changed raw bytes accepted')
print('Exact source readback passed; rejected six export corruptions and changed raw bytes.')
