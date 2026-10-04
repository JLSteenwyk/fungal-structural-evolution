#!/usr/bin/env python3
"""Independent extended-precision raw-kernel reference for the exact failed cohort.

Consumes diagnostic exports only after their producer closes. Groups repeated
species-factor rows and uses wider arithmetic. Disagreement stays diagnostic;
this reader neither relaxes production guards nor qualifies biological fits.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import time

import numpy as np
from scipy import sparse

from ancestral_chain_attempt import sha


def compare(a, b, names):
    a = np.asarray(a, dtype=np.longdouble); b = np.asarray(b, dtype=np.longdouble)
    d = abs(a-b); threshold = np.longdouble(2e-8) + np.longdouble(3e-9)*abs(b)
    return dict(passes_original_allclose=not bool(np.any(d > threshold)),
                maximum_absolute_difference=float(d.max()),
                mismatches=[dict(row=int(i),column=int(j),row_name=names[i],column_name=names[j],
                                 reference=float(a[i,j]),compared=float(b[i,j]),
                                 absolute_difference=float(d[i,j]),tolerance=float(threshold[i,j]))
                            for i,j in np.argwhere(d > threshold)])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--producer', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); assert not a.output.exists(); started = time.monotonic()
    producer = json.loads(a.producer.read_text())
    assert producer['status'] == 'completed_exact_failed_cohort_covariance_diagnostic_v1'
    assert producer['records'] == 22881 and producer['selected_groups'] == 40
    assert producer['candidate_census_rows'] == 1200 and producer['biological_fits'] == 0
    assert producer['production_tolerance_changed'] is False
    assert np.finfo(np.longdouble).eps < np.finfo(float).eps
    for path, checksum in producer['artifacts'].items(): assert sha(path) == checksum, path
    root = Path(json.loads(Path(producer['plan']).read_text())['output'])
    census = json.loads((root/'candidate_census.json').read_text())
    assert len(census) == len({r['candidate_id'] for r in census}) == 1200
    eligible = Counter(r['group_id'] for r in census if r['group_id'])
    assert len(eligible) == 40
    bank = {}; contexts = {}; references = []; records = []
    for item in producer['groups']:
        identity = item['representative']; mode = identity['loading_mode']; tree = identity['tree']
        group = item['group_id']; names = item['representative']['retained_kernel_names']
        assert eligible[group] == item['eligible_candidates']
        assert all(r['benchmark_candidate_id'] == identity['candidate_id'] for r in census if r['group_id'] == group)
        detail = json.loads((root/'groups'/(group+'.json')).read_text())
        if mode not in bank:
            ops = [sparse.load_npz(root/(mode+'-'+name+'.npz')).astype(np.longdouble) for name in names[1:-1]]
            raw = np.zeros((len(names),len(names)),dtype=np.longdouble); raw[0,0] = producer['records']
            for i,z in enumerate(ops,1):
                raw[0,i] = raw[i,0] = z.multiply(z).sum()
                for j,w in enumerate(ops,1):
                    cross = z.T@w
                    raw[i,j] = cross.multiply(cross).sum()
            bank[mode] = (ops,raw,names)
        ops,raw,names = bank[mode]
        key = mode,tree
        if key not in contexts:
            factor = np.load(root/(tree+'-factor.npy'))
            assert factor.dtype == np.dtype('float64') and factor.shape == (producer['records'],301)
            wide = factor.astype(np.longdouble)
            unique, counts = np.unique(factor, axis=0, return_counts=True)
            unique = unique.astype(np.longdouble); counts = counts.astype(np.longdouble)
            core = (unique.T * counts) @ unique
            reference = raw.copy()
            reference[0,-1] = reference[-1,0] = np.sum(counts * np.sum(unique*unique,axis=1))
            reference[-1,-1] = np.sum(core*core)
            for i,z in enumerate(ops,1):
                cross = z.T @ wide
                reference[i,-1] = reference[-1,i] = np.sum(cross*cross)
            contexts[key] = reference
            references.append(dict(loading_mode=mode,tree=tree,names=names,
                                   factor_unique_rows=len(unique),factor_shape=list(factor.shape),
                                   reference_raw_gram=reference.astype(float).tolist()))
            print('extended_precision_reference',mode,tree,'unique_factor_rows',len(unique),flush=True)
        reference = contexts[key]
        results = dict(group_id=group,
                       reference_vs_fresh=compare(reference,detail['fresh']['raw_gram'],names),
                       reference_vs_inherited=compare(reference,detail['inherited']['raw_gram'],names),
                       reference_vs_latent=compare(reference,detail['latent']['raw'],names))
        records.append(results)
    assert len(records) == len({r['group_id'] for r in records}) == 40
    for path, checksum in producer['artifacts'].items(): assert sha(path) == checksum, path
    receipt = dict(status='completed_extended_precision_failed_cohort_raw_reference_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(),producer=str(a.producer),producer_sha256=sha(a.producer),
        source_hashes={str(Path(__file__)):sha(__file__),str(a.producer):sha(a.producer)},
        producer_artifact_bindings_verified=len(producer['artifacts']),
        selected_groups=40,raw_reference_contexts=len(contexts),
        arithmetic_mantissa_bits=np.finfo(np.longdouble).nmant,
        references=references,groups=records,elapsed_seconds=time.monotonic()-started,
        projected_gram_independently_recomputed=False,full_source_archive_rehashed=False,
        scientific_eligibility=False,production_tolerance_changed=False,biological_fits=0,gpu=False,
        scope='Independent raw-kernel arithmetic on original completed diagnostic exports: '
              'sparse latent entity overlaps, longdouble entity/species products and counted '
              'unique factor-row longdouble species Gram. Same allclose tolerances retained. '
              'A wider floating-point reference is not an interval proof, source repair, '
              'projected-kernel acceptance, full original source replay or biological result.')
    with a.output.open('x') as f:json.dump(receipt,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['groups','references','source_hashes']},indent=2))


if __name__ == '__main__': main()
