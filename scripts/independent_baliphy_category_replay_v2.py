"""Full temporal-pattern reconstruction and streamed marginal numerical replay."""
from collections import Counter
from contextlib import ExitStack
import gzip
import hashlib
import json
import os
from pathlib import Path
import time

import numpy as np

from independent_ancestral_categorical_diagnostics_v2 import diagnose_states, compare_states, group_patterns
from independent_ancestral_scalar_diagnostics_v2 import diagnose
from binary_ess_boundary_certificate import compare_with_binary_certificate as compare
from independent_baliphy_category_sources_v2 import category_input, cutoff_input, CUTS, ALPHABET
from run_ortholog_pair_guide_comparison import sha


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def pattern_row(original, pattern, number, multiplicity, plan):
    assert set(original) == {'pattern_id', 'coordinate_multiplicity', 'diagnostic'}
    assert original['pattern_id'] == number and original['coordinate_multiplicity'] == int(multiplicity)
    independent = diagnose_states(pattern, ALPHABET)
    errors, unresolved = compare_states(original['diagnostic'], independent, pattern,
        ALPHABET, **plan['comparison_tolerances'])
    return dict(pattern_id=number, coordinate_multiplicity=int(multiplicity),
        source_diagnostic_sha256=digest(original['diagnostic']), independent=independent,
        maximum_absolute_errors=errors, unresolved_numeric_metrics=unresolved, scientific_eligibility=False)


def replay_cutoff(source, bindings, group, inputs, cutoff, target, plan, verify_output):
    retained, free, summary, files = cutoff_input(source, bindings, inputs, cutoff)
    patterns, mapping = group_patterns(retained)
    assert len(patterns) == summary['patterns']
    with np.load(files['patterns.npz']['path'], allow_pickle=False) as saved:
        assert set(saved.files) == {'patterns', 'coordinate_pattern_ids'}
        assert saved['patterns'].dtype == np.uint8
        assert np.array_equal(saved['patterns'], patterns)
        assert saved['coordinate_pattern_ids'].dtype == np.int32
        assert np.array_equal(saved['coordinate_pattern_ids'], mapping)
    multiplicities = np.bincount(mapping.ravel(), minlength=len(patterns))
    assert len(multiplicities) == len(patterns) and np.all(multiplicities > 0)
    assert int(multiplicities.sum()) == summary['coordinates']
    del retained, mapping
    temporary = target.with_suffix(target.suffix + '.partial')
    maximum = Counter(); statuses = Counter(); coordinate_statuses = Counter(); numeric = Counter()
    singular = boundary = 0; started = time.perf_counter()
    with ExitStack() as stack:
        native = stack.enter_context(gzip.open(files['diagnostics.jsonl.gz']['path'], 'rt'))
        if verify_output:
            exported = stack.enter_context(gzip.open(target, 'rt')); raw = None
        else:
            assert not target.exists()
            raw = stack.enter_context(temporary.open('wb'))
            exported = stack.enter_context(gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0))
        for number, (pattern, multiplicity) in enumerate(zip(patterns, multiplicities)):
            line = native.readline(); assert line, ('Missing original pattern', group, cutoff, number)
            original = json.loads(line)
            row = pattern_row(original, pattern, number, multiplicity, plan)
            if verify_output:
                line = exported.readline(); assert line, ('Missing serialized pattern', group, cutoff, number)
                assert json.loads(line) == row, ('Changed serialized pattern', group, cutoff, number)
            else:
                exported.write(canonical(row) + b'\n')
            state = row['independent']['status']; statuses[state] += 1
            coordinate_statuses[state] += int(multiplicity)
            singular_states = {r['state'] for r in row['unresolved_numeric_metrics'] if r['reason'] == 'rhat_zero_split_within_variance_requires_review'}
            boundary_states = {r['state'] for r in row['unresolved_numeric_metrics'] if r['reason'] == 'exact_binary_zero_pair_truncation_requires_review'}
            boundary += len(boundary_states)
            singular += len(singular_states)
            for label, indicator in row['independent']['indicators'].items():
                screen = indicator['screen']
                if label in boundary_states:
                    numeric['binary_zero_pair_metric_requires_review'] += 1
                elif label in singular_states:
                    numeric['singular_metric_requires_review'] += 1
                elif any(screen.get(k) is not None for k in ['rhat', 'bulk_ess', 'tail_ess', 'mean_mcse']):
                    numeric['all_defined_metrics_compared'] += 1
                else:
                    numeric['original_disposition_only_no_defined_metrics'] += 1
            for name, error in row['maximum_absolute_errors'].items(): maximum[name] = max(maximum[name], error)
            if (number + 1) % 10000 == 0:
                print('independent_categorical_patterns', group, cutoff, number + 1, '/', len(patterns),
                    'seconds', round(time.perf_counter() - started, 3), flush=True)
        assert native.readline() == '', ('Extra original pattern', group, cutoff)
        if verify_output:
            assert exported.readline() == '', ('Extra serialized pattern', group, cutoff)
        else:
            exported.close(); raw.flush(); os.fsync(raw.fileno())
    if not verify_output:
        os.replace(temporary, target)
    assert dict(statuses) == summary['pattern_status_counts']
    assert dict(coordinate_statuses) == summary['coordinate_status_counts']
    nodes = inputs['coordinates']['nodes']
    assert set(summary['unanchored_count_screens']) == set(nodes)
    unanchored = {}
    for index, node in enumerate(nodes):
        original = summary['unanchored_count_screens'][node]; independent = diagnose(free[:, :, index])
        errors, unresolved = compare(original, independent, values=free[:, :, index], **plan['comparison_tolerances'])
        unanchored[node] = dict(original=original, independent=independent,
            absolute_errors=errors, unresolved_numeric_metrics=unresolved, scientific_eligibility=False)
        for name, error in errors.items(): maximum[name] = max(maximum[name], error)
    return dict(patterns=len(patterns), coordinates=summary['coordinates'],
        declared_indicator_rows=len(patterns) * len(ALPHABET),
        original_patterns=files['patterns.npz'], original_diagnostics=files['diagnostics.jsonl.gz'],
        serialized_comparison_sha256=sha(target), pattern_status_counts=dict(statuses),
        coordinate_status_counts=dict(coordinate_statuses), numeric_comparison_status_counts=dict(numeric),
        maximum_absolute_errors=dict(maximum), singular_indicator_rows=singular, binary_zero_pair_indicator_rows=boundary,
        unanchored_count_rows=unanchored, scientific_eligibility=False)


def replay_group(source, bindings, group, info, root, plan, readback=False):
    cutoffs = {}
    if info['status'] != 'unresolved_failed_native_chain_retained':
        inputs = category_input(source, bindings, group, info)
        folder = root / 'groups' / group
        if not readback: folder.mkdir(exist_ok=True)
        for cutoff in CUTS:
            target = folder / ('discard-' + cutoff + '.jsonl.gz')
            assert not readback or target.exists()
            result = replay_cutoff(source, bindings, group, inputs, cutoff, target, plan,
                verify_output=readback or target.exists())
            result['serialized_comparison'] = str(target.relative_to(root)); cutoffs[cutoff] = result
    return dict(status=('complete_categorical_comparison_not_posterior_qualification'
        if cutoffs else 'unresolved_failed_native_chain_retained'),
        chain_ids=info['chain_ids'], cutoffs=cutoffs, scientific_eligibility=False)
