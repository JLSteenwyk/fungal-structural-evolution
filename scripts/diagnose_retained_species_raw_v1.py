#!/usr/bin/env python3
"""Scoped wider-precision species Gram reference while whole-source loading runs.

Verifies consumed factor/row files against saved receipts and compact closure
links. Does not replay the complete original source archive or production guard.
"""
import argparse
import csv
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
import time

import numpy as np

from ancestral_chain_attempt import sha
from full_expanded_model_design_sources import array_digest
from reference_measurement_union_sources import bind, verify


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True)
    a=p.parse_args();plan=json.loads(a.plan.read_text());verify(plan['pins'])
    root=Path(plan['output']);root.mkdir(exist_ok=False);started=time.monotonic();pins=dict(plan['pins'])
    bind(pins,a.plan)
    qpath=Path('metadata/full_uniform_covariance_qualification_plan_20261002.json')
    q=json.loads(qpath.read_text());bind(pins,qpath)
    paths={name:Path(q[name+'_plan']) for name in ['operator','design']}
    configs={name:json.loads(path.read_text()) for name,path in paths.items()}
    for path in paths.values():bind(pins,path)
    designroot=Path(configs['design']['output']);oproot=Path(configs['operator']['output'])
    cp=Path(configs['operator']['covariance_plan']);c=json.loads(cp.read_text());covroot=Path(c['output']);bind(pins,cp)
    # Consume separately closed factor/design/operator producer receipts, without
    # pretending that this diagnostic repeats their complete source closures.
    for name,path in [('design',q['design_completion']),('operator',q['operator_completion']),
                      ('covariance','metadata/full_expanded_covariance_completed_20261002.json')]:
        closure=json.loads(Path(path).read_text());bind(pins,path)
        rp=Path(closure['producer_receipt']);assert sha(rp)==closure['producer_receipt_sha256'];bind(pins,rp)
        receipt=json.loads(rp.read_text())
        if name=='design':dr=receipt
        elif name=='operator':opr=receipt
        else:cr=receipt
    manifestpath=designroot/'cohort_manifest.json';bind(pins,manifestpath,dr['artifacts']['cohort_manifest.json'])
    manifest=json.loads(manifestpath.read_text());cohort=next(c for c in manifest if c['cohort_id']==plan['cohort_id'])
    rowsfile=designroot/cohort['path'];assert sha(rowsfile)==cohort['sha256'];bind(pins,rowsfile)
    with np.load(rowsfile) as f:rows=f['case_rows']
    assert rows.dtype==np.dtype('int64') and len(rows)==22881
    assert array_digest(rows,'<i8')==cohort['case_rows_sha256']
    idfile=oproot/'case_ids.json';assert sha(idfile)==opr['artifacts']['case_ids.json'];bind(pins,idfile)
    ids=json.loads(idfile.read_text())
    assert array_digest([ids[i] for i in rows],'S64')==cohort['ordered_case_ids_sha256']
    indexfile=covroot/'case_covariance_index.tsv.gz';assert sha(indexfile)==cr['artifacts'][indexfile.name];bind(pins,indexfile)
    with gzip.open(indexfile,'rt') as f:index=list(csv.DictReader(f,delimiter='\t'))
    assert [r['case_id'] for r in index]==ids
    patterns=np.asarray([int(r['species_pattern_row']) for r in index],dtype=np.int64)[rows]
    auditfile=Path(q['output'])/'design_covariance_audits.jsonl.gz';bind(pins,auditfile)
    audits=[]
    with gzip.open(auditfile,'rt') as f:
        for line in f:
            r=json.loads(line)
            if r['cohort_id']==cohort['cohort_id']:audits.append(r)
            elif audits:break
    assert len(audits)==300 and all(r['records']==22881 for r in audits)
    (root/'original_cohort_audits.json').write_text(json.dumps(audits,indent=2)+'\n')
    assert np.finfo(np.longdouble).eps < np.finfo(float).eps
    results=[]
    for tree in q['trees']:
        factorfile=covroot/(tree+'.npz');assert sha(factorfile)==cr['artifacts'][factorfile.name];bind(pins,factorfile)
        with np.load(factorfile) as f:factor=f['factor'][patterns]
        assert factor.dtype==np.dtype('float64') and factor.shape==(22881,301)
        unique,counts=np.unique(factor,axis=0,return_counts=True)
        wide=unique.astype(np.longdouble);weights=counts.astype(np.longdouble)
        core=(wide.T*weights)@wide
        reference=np.asarray([np.sum(weights*np.sum(wide*wide,axis=1)),np.sum(core*core)],dtype=np.longdouble)
        direct=np.asarray([np.sum(factor*factor),np.sum((factor.T@factor)**2)],dtype=float)
        old=[]
        for audit in audits:
            if audit['tree']!=tree or audit['numerical_audit'] is None:continue
            numerical=audit['numerical_audit'];species=numerical['kernel_names'].index('species')
            old.append(np.asarray(numerical['raw_gram'],dtype=float)[[0,species],species].tolist())
        assert len(old)==60
        checks=[]
        for values in old:
            delta=abs(reference-np.asarray(values,dtype=np.longdouble))
            threshold=np.longdouble(2e-8)+np.longdouble(3e-9)*abs(np.asarray(values,dtype=np.longdouble))
            checks.append(dict(original_values=values,absolute_difference=delta.astype(float).tolist(),
                               passes_original_allclose=not bool(np.any(delta>threshold))))
        result=dict(tree=tree,factor_shape=list(factor.shape),factor_dtype=str(factor.dtype),
                    unique_factor_rows=len(unique),extended_reference=reference.astype(float).tolist(),
                    fresh_float64=direct.tolist(),reference_vs_fresh_passes_original_allclose=bool(np.all(abs(reference-direct)<=2e-8+3e-9*abs(direct))),
                    original_audits_compared=len(checks),original_raw_species_pairs_distinct=len({tuple(r) for r in old}),
                    inherited_disagreements=sum(not r['passes_original_allclose'] for r in checks),checks=checks)
        results.append(result);print('species_raw_reference',tree,'disagreeing_original_audits',result['inherited_disagreements'],flush=True)
    verify(pins)
    proof=dict(status='completed_scoped_failed_cohort_species_raw_reference_v1',checked_utc=datetime.now(timezone.utc).isoformat(),
        cohort_id=cohort['cohort_id'],records=len(rows),trees_checked=5,original_audits_checked=300,
        raw_entries_per_audit=2,arithmetic_mantissa_bits=int(np.finfo(np.longdouble).nmant),results=results,
        source_hashes=pins,artifacts={str(f):sha(f) for f in root.iterdir() if f.is_file()},
        elapsed_seconds=time.monotonic()-started,full_original_source_archive_replayed=False,
        original_guard_reexecuted=False,projected_gram_recomputed=False,production_tolerance_changed=False,
        scientific_eligibility=False,biological_fits=0,gpu=False,
        scope='All300original audits of exact failed cohort11 inspected for residual/species and '
              'species/species raw entries. Counted unique factor rows and wider floating-point '
              'arithmetic compared with direct float64 and unchanged3e-9/2e-8criterion. Consumed '
              'factor/row inputs match separately closed producer receipts; whole source and '
              'original production guard are not replayed. No source acceptance or original restart.')
    with a.receipt.open('x') as f:json.dump(proof,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in proof.items() if k not in ['results','source_hashes','artifacts']},indent=2))


if __name__=='__main__':main()
