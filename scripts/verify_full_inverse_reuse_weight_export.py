#!/usr/bin/env python3
"""Check complete exported bytes and control redundancy; not a count readback."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import psutil

from ancestral_chain_attempt import sha
from full_exact_covariance_sources import runtime_caps
from inverse_reuse_weight_controls import ARRAYS
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal
from reference_measurement_union_sources import verify
from run_full_inverse_reuse_weights import PRODUCER


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();assert not a.output.exists()
    caps=runtime_caps();proc=psutil.Process()
    identity=dict(pid=proc.pid,created=proc.create_time(),cmdline=proc.cmdline())
    print(json.dumps(dict(original_audit_process=identity)),flush=True)
    pp=Path('metadata/full_inverse_reuse_weight_plan_20261003_v1.json');plan=json.loads(pp.read_text())
    root=Path(plan['output']);rp=root/'receipt.json';receipt=json.loads(rp.read_text())
    assert receipt['status']==PRODUCER and receipt['plan_sha256']==sha(pp)
    assert receipt['scientific_eligibility'] is receipt['raw_reml_basis_qualification_complete'] is False
    assert receipt['residual_diagonal_prepared'] is True
    lp=Path('metadata/full_inverse_reuse_weights_v1_launch_20261003.json');launch=json.loads(lp.read_text());launch['launch']=str(lp)
    assert fingerprint(launch) is None;terminal=journal_terminal(launch)
    verify(receipt['source_hashes'])
    artifacts={str(root/name):h for name,h in receipt['artifacts'].items()};verify(artifacts)
    manifest=json.loads((root/'control_manifest.json').read_text())
    assert len(manifest)==plan['expected']['cohorts']==4340
    assert [r['cohort_id'] for r in manifest]==sorted({r['cohort_id'] for r in manifest})
    equal=0;unequal=0;occurrences=0;maximum_weight_difference=0.;maximum_diagonal_difference=0.
    for r in manifest:
        with np.load(root/r['array_path'],allow_pickle=False) as z:
            assert z.files==ARRAYS
            n=r['records'];assert z['reuse_counts'].shape==(3,n)
            w,d=z['weights'],z['reciprocal_diagonals']
            assert w.shape==d.shape==(4,n) and np.isfinite(w).all() and np.isfinite(d).all()
            same=np.array_equal(w[1],w[2]) and np.array_equal(d[1],d[2])
            equal+=int(same);unequal+=int(not same);occurrences+=n
            maximum_weight_difference=max(maximum_weight_difference,float(np.max(np.abs(w[1]-w[2]))))
            maximum_diagonal_difference=max(maximum_diagonal_difference,float(np.max(np.abs(d[1]-d[2]))))
    assert occurrences==receipt['case_row_occurrences']==34110120
    bindings={**receipt['source_hashes'],**artifacts,str(rp):sha(rp),str(lp):sha(lp),str(Path(__file__)):sha(__file__)}
    verify(bindings)
    result=dict(status='verified_complete_export_hashes_and_all_cohort_control_redundancy_pending_sql_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(),process=identity,actual_cgroup_limits=caps,
        original_producer_terminal=terminal,cohorts=4340,case_row_occurrences=occurrences,
        source_bindings_rehashed=len(receipt['source_hashes']),export_artifacts_rehashed=len(artifacts),
        background_node_pair_exact_equal_cohorts=equal,background_node_pair_different_cohorts=unequal,
        maximum_background_node_pair_weight_difference=maximum_weight_difference,
        maximum_background_node_pair_diagonal_difference=maximum_diagonal_difference,
        source_hashes=bindings,sql_fraction_count_readback_complete=False,
        raw_reml_basis_qualification_complete=False,scientific_eligibility=False,
        scope='Full original producer terminal journal and fresh4384source/8681artifact hashes, with exact exported node/pair weight and diagonal comparisons over every cohort. No independent source-group counting, weighted numerical qualification, fitted effect, biological independence, final process peak or stage closure claim. Original policies and settings remain separate even if numerically identical.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','original_producer_terminal','process','scope']}),flush=True)


if __name__=='__main__':main()
