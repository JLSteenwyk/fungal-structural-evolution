"""Complete closed recovery sources for separate categorical numerical checks."""
from collections import Counter
import json
from pathlib import Path

import numpy as np

from independent_baliphy_scalar_sources import load as scalar_sources
from reference_measurement_union_sources import bind


CUTS = ['250', '500']
ALPHABET = list('ACDEFGHIKLMNPQRSTVWYX-')
SUMMARY_FIELDS = ['full_quartets', 'complete_quartets', 'unresolved_quartets',
    'pattern_cutoff_rows', 'declared_indicator_rows', 'coordinate_counts',
    'pattern_status_counts', 'coordinate_status_counts', 'numeric_comparison_status_counts',
    'maximum_absolute_errors', 'singular_indicator_rows', 'unanchored_count_rows']


def load(plan, path):
    source, bindings = scalar_sources(plan, path)
    producer = source['doc'](source['complete']['producer_receipt'],
        source['complete']['producer_receipt_sha256'])
    source['states'] = producer['states']
    assert set(source['states']) == set(source['native'])
    assert plan['comparison_tolerances'] == dict(atol=1e-8, rtol=1e-8)
    if 'expected' in plan:
        closed = source['complete']
        actual = dict(full_quartets=closed['full_quartets'], complete_quartets=closed['complete_quartets'],
            unresolved_quartets=closed['unresolved_quartets'],
            pattern_cutoff_rows=sum(closed['categorical_pattern_counts'].values()),
            declared_indicator_rows=len(ALPHABET) * sum(closed['categorical_pattern_counts'].values()),
            unanchored_count_rows=8 * closed['complete_quartets'])
        assert actual == plan['expected']
    return source, bindings


def category_input(source, bindings, group, info):
    doc = source['doc']
    category = doc(info['categorical']['receipt'], info['categorical']['receipt_sha256'])
    assert category['status'] == 'verified_quartet_categorical_reports_complete_not_posterior_qualification'
    assert category['group'] == group
    manifest = doc(category['manifest'], category['manifest_sha256'])
    report = doc(category['report'], category['report_sha256'])
    assert manifest['status'] == 'provenance_checked_state_quartet'
    assert report['status'] == 'both_cutoff_categorical_reports_complete_not_posterior_qualification'
    assert report['arviz_version'] == '0.22.0' and report['numpy_version'] == '2.2.6'
    assert report['pins'][category['manifest']] == category['manifest_sha256']
    assert report['pins'][manifest['arrays']] == manifest['arrays_sha256']
    assert set(report['outputs']) == set(CUTS)
    coords = manifest['coordinates']; alphabet = list(coords['alphabet'])
    assert alphabet == ALPHABET
    assert len(coords['nodes']) == len(set(coords['nodes'])) == 4
    assert coords['nodes'] == sorted(coords['nodes'])
    assert len(coords['tips']) == len({r['tip'] for r in coords['tips']})
    assert all(isinstance(r['length'], int) and r['length'] >= 0 for r in coords['tips'])
    total_length = sum(r['length'] for r in coords['tips'])
    assert total_length > 0
    assert manifest['expected_iterations'] == list(range(0, 1001, 10))
    chains = manifest['chains']
    assert len(chains) == len({c['chain_id'] for c in chains}) == len({c['seed'] for c in chains}) == 4
    assert {c['chain_id'] for c in chains} == set(info['chain_ids'])
    bind(bindings, manifest['arrays'], manifest['arrays_sha256'])
    with np.load(manifest['arrays'], allow_pickle=False) as saved:
        assert set(saved.files) == {'values', 'iterations', 'unanchored_residue_counts'}
        values = saved['values']; iterations = saved['iterations']; free = saved['unanchored_residue_counts']
    assert values.dtype == np.uint8 and values.shape == (4, 101, 4, total_length)
    assert np.all(values < len(alphabet))
    assert np.array_equal(iterations, manifest['expected_iterations'])
    assert np.issubdtype(iterations.dtype, np.integer)
    assert free.shape == (4, 101, 4) and np.issubdtype(free.dtype, np.integer) and np.all(free >= 0)
    for number, chain in enumerate(chains):
        cid = chain['chain_id']; native = source['native'][cid]
        assert chain['seed'] == native['chain']['seed'] and chain['model_input_identity'] == group
        audit = doc(native['selected_disposition']['sample_audit'], native['selected_disposition']['sample_audit_sha256'])
        assert Path(chain['log']).resolve() == Path(audit['scalar_log']).resolve()
        assert chain['log_sha256'] == audit['scalar_log_sha256']
        assert coords['input_alignment_sha256'] == native['chain']['alignment_sha256']
        assert coords['nodes'] == sorted({r['source_node'] for r in audit['candidate_samples']})
        state = source['states'][cid]
        receipt = doc(state['receipt'], state['receipt_sha256'])
        assert len(receipt['summaries']) == 1
        summary = receipt['summaries'][0]
        assert summary['source_audit_sha256'] == native['selected_disposition']['sample_audit_sha256']
        assert summary['candidate_anchor_coordinates'] == 4 * total_length
        arrays = [p for p in summary['artifacts'] if Path(p).name == 'states.npz']
        coordinate_docs = [p for p in summary['artifacts'] if Path(p).name == 'coordinates.json']
        assert len(arrays) == len(coordinate_docs) == 1
        assert doc(coordinate_docs[0], summary['artifacts'][coordinate_docs[0]]) == coords
        bind(bindings, arrays[0], summary['artifacts'][arrays[0]])
        with np.load(arrays[0], allow_pickle=False) as saved:
            assert saved['states'].dtype == np.uint8
            assert np.array_equal(saved['states'], values[number])
            assert np.array_equal(saved['iterations'], iterations)
            assert np.array_equal(saved['unanchored_residue_counts'], free[number])
    return dict(category=category, manifest=manifest, report=report, coordinates=coords,
        values=values, iterations=iterations, free=free)


def cutoff_input(source, bindings, inputs, cutoff):
    category, report = inputs['category'], inputs['report']
    root = Path(category['report']).parent
    folder = root / ('discard-' + cutoff)
    summary = source['doc'](folder / 'summary.json')
    retained = inputs['values'][:, inputs['iterations'] > int(cutoff)]
    free = inputs['free'][:, inputs['iterations'] > int(cutoff)]
    output = report['outputs'][cutoff]
    assert summary['discard_through'] == int(cutoff)
    assert summary['retained_samples_per_chain'] == output['retained_samples_per_chain'] == retained.shape[1] == (75 if cutoff == '250' else 50)
    assert summary['coordinates'] == output['coordinates'] == np.prod(retained.shape[2:])
    assert summary['patterns'] == output['patterns']
    files = {}
    for name in ['patterns.npz', 'diagnostics.jsonl.gz']:
        p = folder / name; relative = str(p.relative_to(root)); digest = report['artifacts'][relative]
        bind(bindings, p, digest); files[name] = dict(path=str(p), sha256=digest)
    return retained, free, summary, files


def summarize(results, source):
    assert set(results) == set(source['groups'])
    numeric = Counter(); maximum = Counter(); pattern_counts = Counter(); coords = Counter()
    pattern_status = Counter(); coordinate_status = Counter()
    complete = failed = indicators = singular = unanchored = 0
    for group, result in results.items():
        info = source['groups'][group]
        assert result['chain_ids'] == info['chain_ids'] and result['scientific_eligibility'] is False
        if info['status'] == 'unresolved_failed_native_chain_retained':
            assert result['status'] == info['status'] and result['cutoffs'] == {}
            failed += 1; continue
        assert result['status'] == 'complete_categorical_comparison_not_posterior_qualification'
        assert set(result['cutoffs']) == set(CUTS)
        complete += 1
        for cutoff, row in result['cutoffs'].items():
            assert row['declared_indicator_rows'] == row['patterns'] * len(ALPHABET)
            assert row['scientific_eligibility'] is False
            assert sum(row['pattern_status_counts'].values()) == row['patterns']
            assert sum(row['coordinate_status_counts'].values()) == row['coordinates']
            assert sum(row['numeric_comparison_status_counts'].values()) == row['declared_indicator_rows']
            assert len(row['unanchored_count_rows']) == 4
            pattern_counts[cutoff] += row['patterns']; coords[cutoff] += row['coordinates']
            for counts, key in [(pattern_status, 'pattern_status_counts'), (coordinate_status, 'coordinate_status_counts')]:
                counts.update({cutoff + ':' + k:v for k,v in row[key].items()})
            numeric.update(row['numeric_comparison_status_counts'])
            indicators += row['declared_indicator_rows']; singular += row['singular_indicator_rows']
            unanchored += len(row['unanchored_count_rows'])
            for k,v in row['maximum_absolute_errors'].items(): maximum[k] = max(maximum[k], v)
    closed = source['complete']
    assert complete == 403 and failed == 2 and unanchored == 3224
    assert dict(pattern_counts) == closed['categorical_pattern_counts']
    assert dict(coords) == closed['categorical_coordinate_counts']
    assert dict(pattern_status) == closed['pattern_status_counts']
    assert dict(coordinate_status) == closed['coordinate_status_counts']
    return dict(full_quartets=405, complete_quartets=complete, unresolved_quartets=failed,
        pattern_cutoff_rows=sum(pattern_counts.values()), declared_indicator_rows=indicators,
        coordinate_counts=dict(coords), pattern_status_counts=dict(pattern_status),
        coordinate_status_counts=dict(coordinate_status), numeric_comparison_status_counts=dict(numeric),
        maximum_absolute_errors=dict(maximum), singular_indicator_rows=singular,
        unanchored_count_rows=unanchored)
