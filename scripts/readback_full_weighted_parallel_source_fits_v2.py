#!/usr/bin/env python3
"""Read every original weighted candidate, independent audit and setting link."""
import argparse
from collections import Counter
import fcntl
import gzip
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from full_expanded_model_design_sources import digest
from full_weighted_covariance_qualification import settings_by_cohort
from full_weighted_fit_exports import SCHEMA, PRODUCER, READER, SUMMARY, INDEPENDENT_STATUSES, atomic, finish_gzip, validate_export, key, link_rows, process_links, summary, failure_capture
from full_weighted_shared_entity_fit_sources_parallel_v1 import load, cohorts, cases, operators_for
from readback_weighted_shared_entity_candidate import numeric
from weighted_fit_numeric_memory_guard_v1 import ArrayGuard, guarded_call
from reference_measurement_union_sources import bind, verify


def run(path, output, stop_after_cohorts=None):
    path = Path(path); output = Path(output); plan = json.loads(path.read_text()); source, bindings = load(plan, path)
    original = dict(bindings); root = Path(plan['output']); rp = root / 'receipt.json'; receipt = json.loads(rp.read_text())
    lock = (root / 'readback.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    completed = root / 'independent_readback_completed.json'
    if completed.exists():
        prior = json.loads(completed.read_text()); assert sha(prior['output']) == prior['output_sha256']
        raise AssertionError('Completed independent weighted fitting reader cannot restart')
    assert not output.exists()
    assert receipt['status'] == PRODUCER and receipt['plan_sha256'] == sha(path)
    assert receipt['source_contract'] == source['fit_contract'] and receipt['source_hashes'] == original
    assert receipt['scientific_eligibility'] is receipt['nonuniform_weighting_accepted'] is receipt['component_variance_attribution_accepted'] is False
    assert receipt['cached_numeric_inputs_preserved_for_all_cohorts'] is True and receipt['native_call_numeric_inputs_checked'] is True
    bind(bindings, rp)
    manifest = json.loads((root / 'cohort_manifest.json').read_text())
    assert [m['cohort_id'] for m in manifest] == [c['cohort_id'] for c in source['cohorts']]
    names = {'stage_plan.json', 'cohort_manifest.json'} | {m[k] for m in manifest for k in ['candidates', 'links', 'receipt']}
    assert set(receipt['artifacts']) == names
    for n, h in receipt['artifacts'].items(): bind(bindings, root / n, h)
    verify(bindings)
    cached_guard = ArrayGuard(source, 'cached-source')
    state = dict(schema='full-four-control-guarded-parallel-source-fit-v2', plan_sha256=sha(path), source_contract=source['fit_contract'])
    assert json.loads((root / 'stage_plan.json').read_text()) == state
    independent = root / 'independent'; independent.mkdir(exist_ok=True)
    read_state = dict(schema='full-four-control-fit-spectral-reader-v1', plan_sha256=sha(path), producer_receipt_sha256=sha(rp))
    marker = independent / 'stage_plan.json'
    if marker.exists(): assert json.loads(marker.read_text()) == read_state
    else: atomic(marker, read_state)
    settings, setting_count = settings_by_cohort(source); assert setting_count == source['numerical_completion']['settings']
    counts = Counter(); link_counts = Counter(); audited = Counter(); link_audited = Counter(); total = links = 0; artifacts = {}
    for number, ((cohort, rows, entries), part) in enumerate(zip(cohorts(source, plan), manifest), 1):
        cached_guard.check(root, {'cohort_id':cohort['cohort_id']}, 'cached-before')
        cohort_guard = ArrayGuard((rows, entries), 'cohort-input')
        cid = cohort['cohort_id']
        assert part['candidates'] == 'cohorts/' + cid + '.candidates.jsonl.gz'
        assert part['links'] == 'cohorts/' + cid + '.links.tsv.gz' and part['receipt'] == 'cohorts/' + cid + '.receipt.json'
        cp = root / part['receipt']; checkpoint = json.loads(cp.read_text())
        assert checkpoint['stage'] == state and checkpoint['cohort_id'] == cid
        assert checkpoint['cached_numeric_inputs_preserved'] is True and checkpoint['cohort_numeric_inputs_preserved'] is True
        assert checkpoint['native_call_numeric_inputs_checked'] is True
        assert checkpoint['cached_input_array_bindings_checked'] == len(cached_guard.refs)
        assert checkpoint['cohort_input_array_bindings_checked'] == len(cohort_guard.refs)
        assert checkpoint['candidate_sha256'] == part['candidates_sha256'] == receipt['artifacts'][part['candidates']]
        assert checkpoint['links_sha256'] == part['links_sha256'] == receipt['artifacts'][part['links']]
        assert part['receipt_sha256'] == receipt['artifacts'][part['receipt']]
        ap = independent / (cid + '.audits.jsonl.gz'); lp = independent / (cid + '.links.tsv.gz'); ip = independent / (cid + '.receipt.json')
        prior = json.loads(ip.read_text()) if ip.exists() else None
        if prior:
            assert prior['stage'] == read_state and prior['audits_sha256'] == sha(ap) and prior['links_sha256'] == sha(lp)
            sink = gzip.open(ap, 'rt'); previous = (json.loads(line) for line in sink)
        else:
            assert not ap.exists() and not lp.exists(), 'Preserve incomplete original reader files'
            temporary = ap.with_name(ap.name + '.partial'); sink = gzip.open(temporary, 'xt')
        index = {}; local = Counter(); checked = Counter(); seen = 0; operators = {}
        try:
            with gzip.open(root / part['candidates'], 'rt') as f:
                exports = (json.loads(line) for line in f)
                for expected, matrix, response, source_audit, route, diagonal in cases(source, plan, cohort, entries):
                    value = None
                    try:
                        value = next(exports); validate_export(expected, value)
                        basis_key = expected['loading_mode'], tuple(source_audit['retained_kernel_names'])
                        if basis_key not in operators: operators[basis_key] = operators_for(source, rows, source_audit)
                        # Even reused closed chunks are freshly numerically replayed:
                        # their own newly hashed checkpoint cannot promote a review.
                        audit = guarded_call(numeric, source, plan, rows, expected, matrix, response,
                            operators[basis_key], source_audit, route, diagonal, root=root, exported=value, reader=True)
                        assert audit['disposition'] in INDEPENDENT_STATUSES and audit['scientific_eligibility'] is False
                        record = dict(candidate_id=expected['candidate_id'], source_candidate_sha256=digest(value), **audit)
                        if prior: assert next(previous) == record
                        else: sink.write(json.dumps(record, sort_keys=True, allow_nan=False) + '\n')
                    except (AssertionError, ValueError, ArithmeticError, KeyError, StopIteration) as error:
                        failure_capture(independent / 'failures', expected, value, rows, diagonal, matrix, response,
                            source['factors'][expected['tree']][rows], error)
                        raise
                    k = key(expected); assert k not in index; index[k] = expected, value['disposition'], audit['disposition']
                    local[value['disposition']] += 1; checked[audit['disposition']] += 1; seen += 1
                assert next(exports, None) is None
            if prior: assert next(previous, None) is None
        finally:
            sink.close()
            cached_guard.check(root, {'cohort_id':cid}, 'cached-after')
            cohort_guard.check(root, {'cohort_id':cid}, 'cohort-after')
        assert seen == checkpoint['candidate_rows'] == 4800 and checkpoint['candidate_status_counts'] == dict(local)
        if not prior: finish_gzip(temporary, ap)
        n, setting_counts, setting_audits = process_links(root / part['links'], link_rows(settings[cid], plan, index),
            reader=True, independent_path=lp, index=index)
        assert n == checkpoint['setting_fit_links'] and checkpoint['setting_status_counts'] == dict(setting_counts)
        new = dict(stage=read_state, cohort_id=cid, candidate_rows=seen, setting_fit_links=n,
            independent_candidate_status_counts=dict(checked), independent_setting_status_counts=dict(setting_audits),
            audits_sha256=sha(ap), links_sha256=sha(lp),
            cached_numeric_inputs_preserved=True, cached_input_array_bindings_checked=len(cached_guard.refs),
            cohort_numeric_inputs_preserved=True, cohort_input_array_bindings_checked=len(cohort_guard.refs),
            native_call_numeric_inputs_checked=True)
        if prior: assert new == prior
        else: atomic(ip, new)
        for p in [ap, lp, ip]: artifacts[str(p.relative_to(root))] = sha(p)
        counts.update(local); link_counts.update(setting_counts); audited.update(checked); link_audited.update(setting_audits)
        total += seen; links += n
        print('independent_full_weighted_fit_cohort', number, '/', len(manifest), 'candidates', seen, flush=True)
        if stop_after_cohorts == number: raise InterruptedError('Software weighted readback checkpoint contract')
    result_summary = summary(source, plan, total, links, counts, link_counts)
    assert all(receipt[k] == result_summary[k] for k in SUMMARY)
    artifacts[str(marker.relative_to(root))] = sha(marker)
    for n, h in artifacts.items(): bind(bindings, root / n, h)
    verify(bindings)
    result = dict(status=READER, plan_sha256=sha(path), producer_receipt_sha256=sha(rp), source_contract=source['fit_contract'],
        **result_summary, independent_candidate_status_counts=dict(audited), independent_setting_status_counts=dict(link_audited),
        independent_backend=plan['independent_backend'], artifacts=artifacts, source_hashes=bindings,
        cached_numeric_inputs_preserved_for_all_cohorts=True, native_call_numeric_inputs_checked=True,
        scientific_eligibility=False, nonuniform_weighting_accepted=False, component_variance_attribution_accepted=False, scope=plan['scope'])
    atomic(output, result); atomic(completed, dict(output=str(output), output_sha256=sha(output),
        plan_sha256=sha(path), producer_receipt_sha256=sha(rp)))
    print(json.dumps(result_summary), flush=True); return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True); a = p.parse_args(); run(a.plan, a.output)
