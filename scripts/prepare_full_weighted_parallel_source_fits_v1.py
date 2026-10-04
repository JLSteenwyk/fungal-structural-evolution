#!/usr/bin/env python3
"""Checkpoint every original four-control candidate and original setting link."""
import argparse
from collections import Counter
import fcntl
import gzip
import json
from pathlib import Path
import shutil

from ancestral_chain_attempt import sha
from full_weighted_covariance_qualification import settings_by_cohort
from full_weighted_fit_exports import SCHEMA, PRODUCER, atomic, finish_gzip, validate_export, key, link_rows, process_links, summary, failure_capture
from full_weighted_shared_entity_fit_sources_parallel_v1 import load, cohorts, cases, operators_for
from reference_measurement_union_sources import verify
from weighted_shared_entity_candidate import candidate


def run(path, stop_after_cohorts=None):
    path = Path(path); plan = json.loads(path.read_text()); source, bindings = load(plan, path)
    root = Path(plan['output']); assert shutil.disk_usage(root.parent).free >= plan['resources']['minimum_free_disk_gib'] * 2**30
    root.mkdir(exist_ok=True); lock = (root / 'stage.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    assert not (root / 'receipt.json').exists(), 'Completed full weighted fitting cannot restart'
    state = dict(schema=SCHEMA, plan_sha256=sha(path), source_contract=source['fit_contract'])
    marker = root / 'stage_plan.json'
    if marker.exists(): assert json.loads(marker.read_text()) == state
    else: atomic(marker, state)
    folder = root / 'cohorts'; folder.mkdir(exist_ok=True)
    settings, setting_count = settings_by_cohort(source); assert setting_count == source['numerical_completion']['settings']
    manifest = []; statuses = Counter(); link_statuses = Counter(); total = links = 0
    for number, (cohort, rows, entries) in enumerate(cohorts(source, plan), 1):
        assert shutil.disk_usage(root).free >= plan['resources']['minimum_free_disk_gib'] * 2**30
        cid = cohort['cohort_id']; fp = folder / (cid + '.candidates.jsonl.gz')
        lp = folder / (cid + '.links.tsv.gz'); cp = folder / (cid + '.receipt.json')
        saved = json.loads(cp.read_text()) if cp.exists() else None
        if saved:
            assert saved['stage'] == state and saved['cohort_id'] == cid
            assert saved['candidate_sha256'] == sha(fp) and saved['links_sha256'] == sha(lp)
            stream = gzip.open(fp, 'rt'); exports = (json.loads(line) for line in stream)
        else:
            assert not fp.exists() and not lp.exists(), 'Preserve incomplete original fitting files'
            temporary = fp.with_name(fp.name + '.partial'); stream = gzip.open(temporary, 'xt')
        index = {}; local = Counter(); seen = 0; operators = {}
        try:
            for expected, matrix, y, audit, route, diagonal in cases(source, plan, cohort, entries):
                basis_key = expected['loading_mode'], tuple(audit['retained_kernel_names'])
                if basis_key not in operators: operators[basis_key] = operators_for(source, rows, audit)
                if saved: value = next(exports)
                else:
                    value = None
                    try:
                        value = candidate(source, plan, rows, expected, matrix, y, operators[basis_key], audit, route, diagonal)
                        validate_export(expected, value)
                    except (AssertionError, ValueError, ArithmeticError, KeyError) as error:
                        failure_capture(root / 'failures', expected, value, rows, diagonal, matrix, y,
                            source['factors'][expected['tree']][rows], error)
                        raise
                    stream.write(json.dumps(value, sort_keys=True, allow_nan=False) + '\n')
                validate_export(expected, value)
                k = key(expected); assert k not in index; index[k] = expected, value['disposition']
                local[value['disposition']] += 1; seen += 1
            if saved: assert next(exports, None) is None
        finally: stream.close()
        assert seen == 4800
        if not saved: finish_gzip(temporary, fp)
        n, link_counts, _ = process_links(lp, link_rows(settings[cid], plan, index), reader=bool(saved))
        checkpoint = dict(stage=state, cohort_id=cid, candidate_rows=seen, setting_fit_links=n,
            candidate_status_counts=dict(local), setting_status_counts=dict(link_counts), candidate_sha256=sha(fp), links_sha256=sha(lp))
        if saved: assert saved == checkpoint
        else: atomic(cp, checkpoint)
        manifest.append(dict(cohort_id=cid, candidates=str(fp.relative_to(root)), candidates_sha256=sha(fp),
            links=str(lp.relative_to(root)), links_sha256=sha(lp), receipt=str(cp.relative_to(root)), receipt_sha256=sha(cp)))
        statuses.update(local); link_statuses.update(link_counts); total += seen; links += n
        print('full_weighted_fit_cohort', number, '/', len(source['cohorts']), 'candidates', seen, flush=True)
        if stop_after_cohorts == number: raise InterruptedError('Software weighted fitting checkpoint contract')
    atomic(root / 'cohort_manifest.json', manifest)
    result_summary = summary(source, plan, total, links, statuses, link_statuses); verify(bindings)
    names = ['stage_plan.json', 'cohort_manifest.json'] + [m[k] for m in manifest for k in ['candidates', 'links', 'receipt']]
    result = dict(status=PRODUCER, plan_sha256=sha(path), source_contract=source['fit_contract'], **result_summary,
        artifacts={n:sha(root / n) for n in names}, source_hashes=bindings,
        scientific_eligibility=False, nonuniform_weighting_accepted=False, component_variance_attribution_accepted=False, scope=plan['scope'])
    atomic(root / 'receipt.json', result); print(json.dumps(result_summary), flush=True); return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True); run(p.parse_args().plan)
