"""Closed recovery overlay and raw TSV sources for separate marginal checks."""
from collections import Counter
import csv
import json
from pathlib import Path

import numpy as np
from background_measurement_union_sources import closed_source
from reference_measurement_union_sources import bind, verify


CUTS = ['250', '500']
KINDS = ['scalar', 'length']
SUMMARY_FIELDS = ['full_quartets', 'complete_quartets', 'unresolved_quartets',
    'variable_cutoff_rows', 'source_status_counts', 'numeric_comparison_status_counts',
    'maximum_absolute_errors', 'quartets_passing_every_scalar', 'quartets_passing_every_length_scalar']


def numeric_status(row):
    if row['unresolved_numeric_metrics']:
        return 'singular_metric_requires_review'
    return ('all_defined_metrics_compared' if 'rhat' in row['original']
            else 'original_disposition_only_no_defined_metrics')


def load(plan, path):
    bindings = dict(plan['pins']); bind(bindings, path)
    complete = closed_source(plan['diagnostic_completion'],
        'complete_verified_full_baliphy_recovery_diagnostics',
        'complete_verified_full_baliphy_recovery_diagnostics_archive', 2, bindings)
    lookup = {str(Path(p).resolve()): d for p, d in bindings.items()}
    def doc(p, expected=None):
        p = Path(p); d = lookup[str(p.resolve())]
        assert expected is None or expected == d
        bind(bindings, p, d)
        return json.loads(p.read_text())
    producer = doc(complete['producer_receipt'], complete['producer_receipt_sha256'])
    recovery = doc(plan['recovery_completion'])
    overlay = doc(recovery['overlay'], recovery['overlay_sha256'])
    native = {r['chain']['chain_id']: r for r in overlay['rows']}
    assert len(native) == len(overlay['rows']) == 1620
    assert len(producer['groups']) == complete['full_quartets'] == 405
    assert complete['complete_quartets'] == 403 and complete['unresolved_quartets'] == 2
    assert overlay['original_samples_concatenated'] is False
    assert overlay['scientific_eligibility'] is False
    seen = set()
    for group, info in producer['groups'].items():
        assert len(info['chain_ids']) == 4 and len(set(info['chain_ids'])) == 4
        assert not seen.intersection(info['chain_ids']); seen.update(info['chain_ids'])
        members = [native[c] for c in info['chain_ids']]
        assert all(m['model_input_identity'] == group for m in members)
        assert len({m['chain']['seed'] for m in members}) == 4
        ready = all(m['selected_disposition']['status'] == 'all_saved_alignments_and_candidate_nodes_checked' for m in members)
        assert info['status'] == ('complete_scalar_length_category_screens_not_posterior_qualification'
            if ready else 'unresolved_failed_native_chain_retained')
    assert seen == set(native)
    verify(bindings)
    return dict(groups=producer['groups'], native=native, complete=complete, doc=doc), bindings


def trace_input(source, bindings, group, info, kind, cutoff, cache):
    doc = source['doc']
    record = doc(info[kind]['receipt'], info[kind]['receipt_sha256'])
    output = record['outputs'][cutoff]; report = doc(output['path'], output['sha256'])
    assert report['status'] == 'scalar_diagnostics_complete_not_posterior_qualification'
    assert report['arviz_version'] == '0.22.0' and report['numpy_version'] == '2.2.6'
    assert report['discard_through_iteration'] == int(cutoff)
    assert report['thresholds'] == dict(rhat_strict_upper=1.01, bulk_ess_minimum=400, tail_ess_minimum=400)
    manifests = [p for p in report['pins'] if Path(p).name == 'manifest-' + cutoff + '.json']
    assert len(manifests) == 1
    manifest = doc(manifests[0], report['pins'][manifests[0]])
    expected = list(range(1001)) if kind == 'scalar' else list(range(0, 1001, 10))
    assert manifest['expected_iterations'] == expected and manifest['discard_through_iteration'] == int(cutoff)
    variables = manifest['variables']; assert len(set(variables)) == len(variables) == (37 if kind == 'scalar' else 4)
    assert set(report['variables']) == set(variables)
    chains = manifest['chains']; assert len(chains) == 4
    assert {r['chain_id'] for r in chains} == set(info['chain_ids'])
    assert len({str(Path(r['log']).resolve()) for r in chains}) == 4
    retained = []
    for chain in chains:
        original = source['native'][chain['chain_id']]
        assert chain['model_input_identity'] == group and chain['seed'] == original['chain']['seed']
        assert report['pins'][chain['log']] == chain['log_sha256']
        bind(bindings, chain['log'], chain['log_sha256'])
        if kind == 'scalar':
            selected = original['selected_disposition']
            audit = doc(selected['sample_audit'], selected['sample_audit_sha256'])
            assert Path(audit['scalar_log']).resolve() == Path(chain['log']).resolve()
            assert audit['scalar_log_sha256'] == chain['log_sha256']
        key = (chain['log'], chain['log_sha256'], tuple(variables), tuple(expected))
        if key not in cache:
            with Path(chain['log']).open() as f:
                reader = csv.reader(f, delimiter='\t'); header = next(reader)
                assert len(header) == len(set(header)) and all(v in header for v in variables)
                indices = [header.index(v) for v in variables]; iteration = header.index('iter')
                records = list(reader)
            assert all(len(r) == len(header) for r in records)
            assert [int(r[iteration]) for r in records] == expected
            cache[key] = np.asarray([[float(r[i]) for i in indices] for r in records])
        retained.append(cache[key][np.asarray(expected) > int(cutoff)])
    values = np.asarray(retained)
    assert values.shape == (4, (750 if cutoff == '250' else 500) if kind == 'scalar' else (75 if cutoff == '250' else 50), len(variables))
    return variables, values, report, output


def summarize(groups, source):
    assert set(groups) == set(source['groups'])
    statuses = Counter(); numeric = Counter(); maximum = Counter(); total = complete = failed = 0
    passes = {kind: {cutoff: 0 for cutoff in CUTS} for kind in KINDS}
    for group, result in groups.items():
        original = source['groups'][group]
        assert result['chain_ids'] == original['chain_ids'] and result['scientific_eligibility'] is False
        if original['status'] == 'unresolved_failed_native_chain_retained':
            assert result['status'] == original['status'] and not result['rows']
            failed += 1; continue
        complete += 1
        rows = result['rows']; assert len(rows) == 82
        assert len({(r['kind'], r['cutoff'], r['variable']) for r in rows}) == 82
        for row in rows:
            assert row['kind'] in KINDS and row['cutoff'] in CUTS
            assert row['original']['status'] == row['independent']['status']
            assert row['scientific_eligibility'] is False
            assert row['numeric_status'] == numeric_status(row)
            statuses[row['kind'] + ':' + row['cutoff'] + ':' + row['original']['status']] += 1
            numeric[row['numeric_status']] += 1
            for name, value in row['absolute_errors'].items():
                assert np.isfinite(value) and value >= 0
                maximum[name] = max(maximum[name], value)
        for kind in KINDS:
            for cutoff in CUTS:
                subset = [r for r in rows if r['kind'] == kind and r['cutoff'] == cutoff]
                assert len(subset) == (37 if kind == 'scalar' else 4)
                passes[kind][cutoff] += all(r['original']['status'] == 'passes_scalar_screen_only' for r in subset)
        total += len(rows)
    closed = source['complete']
    assert complete == 403 and failed == 2 and total == 33046
    for kind in KINDS:
        for key, count in closed[kind + '_status_counts'].items():
            assert statuses[kind + ':' + key] == count
        observed = {k.split(':', 1)[1]: v for k, v in statuses.items() if k.startswith(kind + ':')}
        assert observed == closed[kind + '_status_counts']
    assert passes['scalar'] == closed['quartets_passing_every_scalar']
    assert passes['length'] == closed['quartets_passing_every_length_scalar']
    return dict(full_quartets=405, complete_quartets=complete, unresolved_quartets=failed,
        variable_cutoff_rows=total, source_status_counts=dict(statuses),
        numeric_comparison_status_counts=dict(numeric), maximum_absolute_errors=dict(maximum),
        quartets_passing_every_scalar=passes['scalar'], quartets_passing_every_length_scalar=passes['length'])
