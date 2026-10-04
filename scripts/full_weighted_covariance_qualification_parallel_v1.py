#!/usr/bin/env python3
"""Fresh complete-grid weighted raw/REML qualification and latent readback.

Every original setting is retained. Segments bound individual output files;
they do not select cohorts, controls, trees, settings or numerical cases.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import multiprocessing
import hashlib
import fcntl
import resource
import psutil
from scipy import sparse
from collections import Counter
import csv
from datetime import datetime, timezone
import gzip
import itertools
import json
import os
from pathlib import Path

import numpy as np

from ancestral_chain_attempt import sha
from full_exact_covariance_sources import closed_subset
from full_expanded_model_design_sources import AXES, DEGREES, SETTING_FIELDS, array_digest, digest
from full_expanded_model_input_sources import ORDERS
from full_weighted_covariance_sources_v2 import MODES, POLICIES, basis, cohort_controls, design_matrix, jsonl, load, retained_operators
from independent_positive_diagonal_basis_context import IndependentPositiveDiagonalBasisContext
from positive_diagonal_basis_context import PositiveDiagonalBasisContext
from readback_full_covariance_qualification import numeric
from reference_measurement_union_sources import bind, verify


from full_weighted_covariance_qualification import (SCHEMA, SUMMARY, LINK_FIELDS, source_inputs, runtime, settings_by_cohort, expected_record, disposition, failure_capture, write, PRODUCER, READER)

CACHE = None


def arrays(value, prefix='source'):
    if isinstance(value, np.ndarray):
        yield prefix, value
    elif sparse.issparse(value):
        for key in ['data', 'indices', 'indptr']:
            yield prefix+'.'+key, getattr(value, key)
    elif isinstance(value, dict):
        for key, child in value.items():
            yield from arrays(child, prefix+'/'+str(key))
    elif isinstance(value, (list, tuple)):
        for index, child in enumerate(value):
            yield from arrays(child, prefix+'/'+str(index))


def array_bytes(value):
    h = hashlib.sha256()
    h.update(json.dumps([str(value.dtype), list(value.shape)]).encode())
    if value.dtype.hasobject:
        h.update(json.dumps(value.tolist(), sort_keys=True).encode())
    else:
        storage = np.ascontiguousarray(value)
        if storage.size: h.update(memoryview(storage).cast('B'))
    return h.hexdigest()


def worker_init(source, plan, root, reader):
    global CACHE
    limit = plan['resources']['worker_address_space_gib']*2**30
    resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    refs = list(arrays(source)); assert refs
    before = {name:array_bytes(value) for name,value in refs}
    CACHE = (source, plan, Path(root), reader, refs, before)


def cohort_job(task):
    source, plan, root, reader, refs, before = CACHE
    ci, cohort, control, cohort_designs, settings = task
    assert {name:array_bytes(value) for name,value in refs} == before, 'Worker cached numeric source changed before cohort'
    setting_groups = {cohort['cohort_id']:settings}
    counts = Counter(); link_counts = Counter(); policy_counts = Counter(); artifacts = {}
    audit_n = link_n = design_n = conservative_reviews = 0; maximum_error = 0.
    selected, diagonals, control_record = cohort_controls(source, cohort, control)
    assert {(d['order_contrast'], d['sequence_axis'], d['degree']) for d in cohort_designs} == set(itertools.product(ORDERS, AXES, DEGREES))
    matrices = {d['design_id']: design_matrix(source, cohort, selected, d) for d in cohort_designs}
    design_n += len(cohort_designs); statuses = {}; cohort_counts = Counter(); cohort_link_counts = Counter()
    ap = root / 'cohorts' / (cohort['cohort_id'] + '.audits.jsonl.gz')
    lp = root / 'cohorts' / (cohort['cohort_id'] + '.links.tsv.gz')
    exports = jsonl(ap) if reader else None
    out = None if reader else gzip.open(ap, 'wt')
    try:
        for mode, tree in itertools.product(MODES, plan['trees']):
            routes = [basis(source, cohort, mode, d) for d in diagonals]
            contexts = {}; operators_by_basis = {}
            factor = source['factors'][tree][selected]
            for route in routes:
                key = tuple(route['names'])
                if key not in contexts:
                    ops = retained_operators(source, selected, mode, route['names']); operators_by_basis[key] = ops
                    cls = IndependentPositiveDiagonalBasisContext if reader else PositiveDiagonalBasisContext
                    contexts[key] = cls(source['labels'][selected], ops, factor)
            for d in cohort_designs:
                x = matrices[d['design_id']]; projected = {}; context_errors = {}
                if d['disposition'] == 'full_rank_design':
                    for key, context in contexts.items():
                        try: projected[key] = context.design(x)
                        except (ValueError, ArithmeticError) as e: context_errors[key] = e
                for j, policy in enumerate(POLICIES):
                    route = routes[j]; key = tuple(route['names']); diagonal = diagonals[j]
                    expected = expected_record(source, cohort, d, route, diagonal, control, mode, tree, policy)
                    if reader:
                        saved = next(exports); reference = None
                        try:
                            assert all(saved[k] == v for k,v in expected.items())
                            if d['disposition'] != 'full_rank_design':
                                assert saved['disposition'] == d['disposition'] and saved['numerical_audit'] is None
                            elif saved['numerical_audit'] is None:
                                assert saved['disposition'] == 'numerical_covariance_qualification_requires_review'
                                primary = PositiveDiagonalBasisContext(source['labels'][selected], operators_by_basis[key], factor)
                                try: primary.design(x).audit(diagonal)
                                except (ValueError, ArithmeticError) as e:
                                    assert saved['error_type'] == type(e).__name__ and saved['error_message'] == str(e)
                                else: raise AssertionError('False primary numerical failure')
                            else:
                                assert key in projected, 'Independent design requires review'
                                reference = projected[key].audit(diagonal); value = saved['numerical_audit']
                                assert value['fixed_effect_columns'] == x.shape[1] and value['residual_dimension'] == len(selected) - x.shape[1]
                                assert value['family_components'] == len(np.unique(source['labels'][selected]))
                                assert value['exactly_zero_incidence_names'] == [k for k,z in operators_by_basis[key].items() if not z.nnz]
                                # Frozen verifier uses supplied kernels and fresh D bounds;
                                # its historical uniform disposition name is never exported.
                                review, error = numeric(value, reference, route['names'], len(selected))
                                conservative_reviews += int(review); maximum_error = max(maximum_error, error)
                                assert saved['disposition'] == disposition(value)
                        except (AssertionError, ValueError, ArithmeticError, KeyError) as e:
                            failure_capture(root, expected, saved, reference, selected, diagonal, x, factor, e)
                            raise
                    else:
                        saved = dict(expected, numerical_audit=None)
                        if d['disposition'] != 'full_rank_design': saved['disposition'] = d['disposition']
                        else:
                            try:
                                if key in context_errors: raise context_errors[key]
                                saved['numerical_audit'] = projected[key].audit(diagonal)
                                saved['disposition'] = disposition(saved['numerical_audit'])
                            except (ValueError, ArithmeticError) as e:
                                saved.update(disposition='numerical_covariance_qualification_requires_review',
                                    error_type=type(e).__name__, error_message=str(e))
                        out.write(json.dumps(saved, sort_keys=True, allow_nan=False) + '\n')
                    status = saved['disposition']; aid = expected['audit_id']
                    assert (d['design_id'], mode, tree, policy) not in statuses
                    statuses[d['design_id'], mode, tree, policy] = aid, status
                    counts[status] += 1; cohort_counts[status] += 1; policy_counts[policy + ':' + status] += 1; audit_n += 1
    finally:
        if out is not None: out.close()
    if reader: assert next(exports, None) is None
    assert len(statuses) == 1200
    def expected_links():
        for ordinal, row in setting_groups[cohort['cohort_id']]:
            setting = dict(zip(SETTING_FIELDS, row)); assert setting['cohort_id'] == cohort['cohort_id']
            for mode, tree, policy in itertools.product(MODES, plan['trees'], POLICIES):
                aid, status = statuses[setting['design_id'], mode, tree, policy]
                combined = setting['disposition'] if setting['disposition'] != 'ready_for_working_covariance_fit' else status
                yield [str(ordinal), digest(row), *row, mode, tree, policy, aid, status, combined]
    with gzip.open(lp, 'rt' if reader else 'wt') as f:
        if reader:
            links = csv.reader(f, delimiter='\t'); assert next(links) == LINK_FIELDS
        else:
            links = csv.writer(f, delimiter='\t', lineterminator='\n'); links.writerow(LINK_FIELDS)
        for row in expected_links():
            if reader: assert next(links) == row
            else: links.writerow(row)
            link_counts[row[-1]] += 1; cohort_link_counts[row[-1]] += 1; link_n += 1
        if reader: assert next(links, None) is None
    entry = dict(cohort_id=cohort['cohort_id'], audits=1200, links=len(setting_groups[cohort['cohort_id']]) * 40,
        audit_file=str(ap.relative_to(root)), audit_sha256=sha(ap), link_file=str(lp.relative_to(root)), link_sha256=sha(lp),
        audit_status_counts=dict(cohort_counts), link_status_counts=dict(cohort_link_counts))
    after = {name:array_bytes(value) for name,value in refs}
    if after != before:
        failure = root/'failures'/(cohort['cohort_id']+'-source-mutation')
        failure.mkdir(parents=True, exist_ok=False)
        write(failure/'failure.json', dict(before=before, after=after, cohort_id=cohort['cohort_id'], scientific_eligibility=False))
        for name,value in refs:
            if before[name] != after[name]:
                with (failure/(hashlib.sha256(name.encode()).hexdigest()+'.npy')).open('xb') as f: np.save(f,value,allow_pickle=False)
        raise AssertionError('Worker cached numeric source changed during cohort')
    process = psutil.Process()
    checkpoint = dict(cohort_index=ci, cohort_id=cohort['cohort_id'], entry=entry,
        counts=dict(counts), link_counts=dict(link_counts), policy_counts=dict(policy_counts),
        audits=audit_n, links=link_n, designs=design_n,
        conservative_reviews=conservative_reviews, maximum_error=maximum_error,
        worker=dict(pid=process.pid, created=process.create_time(), cmdline=process.cmdline()),
        worker_address_space_limit_bytes=resource.getrlimit(resource.RLIMIT_AS)[0],
        cached_input_array_bindings_checked=len(refs), cached_numeric_inputs_preserved=True,
        scientific_eligibility=False)
    write(root/'checkpoints'/((str(ci).zfill(5))+('.reader.json' if reader else '.producer.json')),checkpoint)
    return checkpoint


def run(path, reader=False):
    """Hold a stage lock, including during source loading and worker teardown."""
    path = Path(path)
    plan = json.loads(path.read_text())
    lock_path = Path(str(plan['output']) + '.parallel-stage.lock')
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        return run_locked(path, reader)


def run_locked(path, reader=False):
    path=Path(path);plan=json.loads(path.read_text());runtime(plan)
    r=plan['resources'];assert 1 <= r['workers'] <= 16
    assert r['worker_address_space_gib'] == 12 and r['worker_address_space_gib']*r['workers'] <= r['reservation_capacity_gib']
    assert r['reservation_capacity_gib']+r['parent_headroom_gib'] <= r['memory_gib']
    assert r['workers'] <= r['cpus']
    assert psutil.disk_usage('.').free >= r['minimum_free_disk_gib']*2**30
    source,bindings=source_inputs(plan,path);original=dict(bindings)
    root=Path(plan['output']);receipt_name='readback.json' if reader else 'receipt.json'
    assert not (root/receipt_name).exists()
    if reader:
        receipt=json.loads((root/'receipt.json').read_text())
        assert receipt['status']==PRODUCER and receipt['source_contract']==source['contract']
        assert receipt['plan_sha256']==sha(path) and receipt['source_hashes']==original
        assert not receipt['scientific_eligibility']
        bind(bindings,root/'receipt.json')
        for n,h in receipt['artifacts'].items():bind(bindings,root/n,h)
        manifest=json.loads((root/'cohort_manifest.json').read_text())
        assert [row['cohort_id'] for row in manifest]==[c['cohort_id'] for c in source['cohorts']]
        assert json.loads((root/'stage_plan.json').read_text())==dict(schema=SCHEMA,plan_sha256=sha(path),source_contract=source['contract'])
    else:
        assert psutil.virtual_memory().available >= r['memory_gib']*2**30
        root.mkdir(parents=True,exist_ok=False);(root/'cohorts').mkdir();(root/'checkpoints').mkdir()
        write(root/'stage_plan.json',dict(schema=SCHEMA,plan_sha256=sha(path),source_contract=source['contract']))
    setting_groups,settings_n=settings_by_cohort(source);assert settings_n==plan['expected']['settings']
    design_stream=jsonl(source['root']/'unique_designs.jsonl')
    tasks=[(ci,c,ctl,[next(design_stream) for _ in range(30)],setting_groups[c['cohort_id']])
        for ci,(c,ctl) in enumerate(zip(source['cohorts'],source['controls']))]
    assert next(design_stream,None) is None and len(tasks)==plan['expected']['cohorts']
    assert len({t[1]['cohort_id'] for t in tasks})==len(tasks)
    outcomes={};remaining=iter(tasks)
    # A bounded pending queue never admits the whole grid before a failure.
    with ProcessPoolExecutor(max_workers=r['workers'],mp_context=multiprocessing.get_context('fork'),
        initializer=worker_init,initargs=(source,plan,str(root),reader)) as pool:
        pending={}
        def submit():
            task=next(remaining,None)
            if task is not None:pending[pool.submit(cohort_job,task)]=task[0]
        for _ in range(r['workers']):submit()
        try:
            while pending:
                done,_=wait(pending,return_when=FIRST_COMPLETED)
                for future in done:
                    ci=pending.pop(future);value=future.result();assert value['cohort_index']==ci and ci not in outcomes
                    outcomes[ci]=value
                    print('parallel_weighted_covariance_cohorts',len(outcomes),'/',len(tasks),'reader',reader,flush=True)
                # Consume every completion before submitting another cohort.
                for _ in done: submit()
        except BaseException:
            for future in pending:future.cancel()
            # Already running tasks can finish; no new cohort is submitted.
            raise
    assert sorted(outcomes)==list(range(len(tasks)))
    counts=Counter();link_counts=Counter();policy_counts=Counter();artifacts={}
    audit_n=link_n=design_n=conservative_reviews=0;maximum_error=0.;entries=[]
    for ci in sorted(outcomes):
        value=outcomes[ci];entry=value['entry'];entries.append(entry)
        if reader:assert manifest[ci]==entry
        counts.update(value['counts']);link_counts.update(value['link_counts']);policy_counts.update(value['policy_counts'])
        audit_n+=value['audits'];link_n+=value['links'];design_n+=value['designs']
        conservative_reviews+=value['conservative_reviews'];maximum_error=max(maximum_error,value['maximum_error'])
        artifacts[entry['audit_file']]=entry['audit_sha256'];artifacts[entry['link_file']]=entry['link_sha256']
        producer_checkpoint=root/'checkpoints'/(str(ci).zfill(5)+'.producer.json')
        saved_checkpoint=json.loads(producer_checkpoint.read_text())
        assert saved_checkpoint['entry']==entry and saved_checkpoint['cached_numeric_inputs_preserved']
        assert saved_checkpoint['cohort_index']==ci and saved_checkpoint['worker_address_space_limit_bytes']==12*2**30
        for field in ['counts','link_counts','policy_counts','audits','links','designs']:
            assert saved_checkpoint[field]==value[field]
        artifacts[str(producer_checkpoint.relative_to(root))]=sha(producer_checkpoint)
    mp=root/'cohort_manifest.json'
    if not reader:write(mp,entries)
    artifacts['cohort_manifest.json']=sha(mp);artifacts['stage_plan.json']=sha(root/'stage_plan.json')
    summary=dict(logical_cases=len(source['ids']),cohorts=len(source['cohorts']),designs=design_n,settings=settings_n,
        numerical_audit_rows=audit_n,setting_audit_links=link_n,audit_status_counts=dict(counts),
        link_status_counts=dict(link_counts),policy_audit_status_counts=dict(policy_counts),policies=POLICIES,trees=plan['trees'],loading_modes=MODES)
    assert design_n==plan['expected']['designs'] and audit_n==plan['expected']['numerical_audit_rows']
    assert link_n==plan['expected']['setting_audit_links']
    if reader:assert all(receipt[k]==v for k,v in summary.items()) and receipt['artifacts']==artifacts
    verify(bindings)
    result=dict(status=READER if reader else PRODUCER,checked_utc=datetime.now(timezone.utc).isoformat(),
        plan_sha256=sha(path),source_contract=source['contract'],**summary,source_hashes=bindings,artifacts=artifacts,
        parallel_workers=r['workers'],worker_address_space_gib=r['worker_address_space_gib'],
        cached_numeric_inputs_preserved_for_all_cohorts=True,
        working_model_fits_computed=0,scientific_eligibility=False,nonuniform_weighting_accepted=False,scope=plan['scope'])
    if reader:result.update(producer_receipt_sha256=sha(root/'receipt.json'),conservative_diagnostic_disagreements=conservative_reviews,maximum_projected_gram_difference=maximum_error,
        reader_checkpoint_artifacts={str(p.relative_to(root)):sha(p) for p in sorted((root/'checkpoints').glob('*.reader.json'))})
    write(root/receipt_name,result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts','scope']}),flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--reader',action='store_true')
    a=p.parse_args();run(a.plan,a.reader)
