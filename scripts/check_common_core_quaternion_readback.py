#!/usr/bin/env python3
"""Independent quaternion fit checks, including output corruption rejection."""
import numpy as np
from fit_domain_triad_common_residues import fit_triplet
from readback_domain_triad_common_fits import core_metrics,compare_row
rng=np.random.default_rng(1981)
for n in [3,7,50,100]:
    reference=rng.normal(size=(n,3));a=reference+rng.normal(size=(n,3))*.1;b=reference+rng.normal(size=(n,3))*.5
    for reflected in [False,True]:
        coords=[a,b,reference*np.array([-1,1,1]) if reflected else reference]
        letters=['A'*n,'C'*n,'A'*n];confidence=[np.full(n,v) for v in [90.,60.,80.]]
        observed=fit_triplet(coords,letters,confidence);wanted=core_metrics(coords,letters,confidence)
        compare_row({k:str(v) for k,v in observed.items()},wanted)
        for field in ['rmsd_ar','rmsd_ar_minus_br','sequence_identity_ab','joint_plddt70_fraction']:
            bad={k:str(v) for k,v in observed.items()};bad[field]=str(float(bad[field])+.1)
            try:compare_row(bad,wanted)
            except ValueError:pass
            else:raise AssertionError('Corruption accepted: '+field)
print('Quaternion and SVD fits agree for small/large/noisy/reflected cores; altered RMSD, contrast, identity and confidence rejected.')
