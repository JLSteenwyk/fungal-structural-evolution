#!/usr/bin/env python3
"""Independently verify platform RMSDs with BioPython and inspect confidence masks."""
import json
from pathlib import Path
import numpy as np
from Bio.SVDSuperimposer import SVDSuperimposer
from compare_spark_qualification import coords, sha

source=Path('metadata/dgx_spark_qualification_comparison.json')
rows=json.loads(source.read_text())['rows']
inputs=json.loads(Path('metadata/dgx_spark_qualification_input_receipt.json').read_text())
checks=[]
for item,row in zip(inputs['sequences'],rows,strict=True):
    assert item['sequence_id']==row['sequence_id']
    local=Path(item['local_receipt']).with_suffix('.pdb')
    remote=Path('results/predictions/dgx-spark-qualification-v1')/local.name
    _,x=coords(local);_,y=coords(remote)
    fit=SVDSuperimposer();fit.set(x,y);fit.run()
    error=abs(fit.get_rms()-row['ca_superposition_rmsd_angstrom']);assert error<1e-10
    with np.load(local.with_suffix('.npz')) as a,np.load(remote.with_suffix('.npz')) as b:
        mask=(a['ca_plddt']>=70)&(b['ca_plddt']>=70)
    qualified=None
    if mask.sum()>=3:
        fit.set(x[mask],y[mask]);fit.run();qualified=float(fit.get_rms())
    checks.append({'sequence_id':row['sequence_id'],'full_rmsd_independent_error':error,
                   'joint_plddt70_ca_count':int(mask.sum()),
                   'joint_plddt70_ca_superposition_rmsd_angstrom':qualified})
result={'status':'independent_BioPython_SVD_check_passed','source_sha256':sha(source),
        'script_sha256':sha(Path(__file__)),'maximum_full_rmsd_error':max(r['full_rmsd_independent_error'] for r in checks),
        'checks':checks,'interpretation':'Confidence filtering is descriptive, not domain segmentation or proof of hardware/software causality.'}
Path('metadata/dgx_spark_qualification_independent_readback.json').write_text(json.dumps(result,indent=2)+'\n')
print(result['maximum_full_rmsd_error'])
