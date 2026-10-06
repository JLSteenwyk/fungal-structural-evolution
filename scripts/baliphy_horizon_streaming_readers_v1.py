"""Horizon-aware streaming readers; integrity never implies posterior admission.

Qualified short-run sources remain unchanged. These readers hold at most one
scalar row or joint frame, plus its arrays. Native custody, resource caps and
long-run production admission remain controller responsibilities.
"""
from collections import Counter
import csv
from decimal import Decimal
from itertools import zip_longest
import json
import math
import os
from pathlib import Path

from ancestral_chain_attempt import sha
from baliphy_scalar_json_logger_v6c import SCHEMA, validate_record
from independent_joint_ancestral_frames import decode, verify_arrays, write_arrays
from independent_native_ancestral_alignment import fasta_records, native_alignments
from independent_native_ancestral_topology import match_trees
from independent_short_sampler_outputs_v2 import NATIVE_ALPHABET, mapping_rows, strict_json
from read_baliphy_scalar_json_v6b import load, read_record
from readback_baliphy_log_alpha_v10 import HEADER, check as latent_check, load as decimal_load


def horizon(iterations, interval=10):
    assert type(iterations) is int and iterations >= 0
    assert type(interval) is int and interval > 0
    return range(0, iterations + 1, interval)


def scalar_records(path, iterations, parser=load):
    horizon(iterations)
    with Path(path).open() as handle:
        header = next(handle, None)
        assert header is not None and parser(header) == HEADER
        count = 0
        for count, line in enumerate(handle, 1):
            assert count <= iterations + 1, 'Extra scalar row'
            record = parser(line)
            assert type(record['iter']) is int and record['iter'] == count - 1, 'Missing/duplicate/unordered scalar iteration'
            yield record
        assert count == iterations + 1, 'Missing scalar row'


def flattened(row):
    flat = {'iter': row['iter']}

    def visit(value, prefix):
        for key, item in value.items():
            if isinstance(item, dict):
                if key.endswith('/'):
                    visit(item, prefix + key)
                else:
                    for subkey, scalar in item.items():
                        name = prefix + key + '[' + subkey + ']'
                        assert name not in flat
                        flat[name] = scalar
            else:
                assert prefix + key not in flat
                flat[prefix + key] = item

    visit({k: v for k, v in row['statistics//'].items() if k != '__project_scalar_v6_quality__'}, '')
    visit(row['parameters//'], '')
    return flat


def scalar_trace(directory, iterations, review_sink=None):
    """First pass retains tags; TSV mapping is checked only for a finite role.

    A sink receives each review separately, avoiding a horizon-sized review list.
    Original 2e-13 relative/zero absolute mapping and 1e-7 score checks remain.
    """
    directory = Path(directory)
    path = directory / 'C1.log.json'
    counts = Counter()
    literal_rows = 0
    for row in scalar_records(path, iterations):
        primary, independent = validate_record(row), read_record(row)
        assert primary['iteration'] == independent['iteration'] == row['iter']
        assert primary['finite_record'] == independent['finite_record']
        assert primary['context']['numeric_leaves'] + primary['parameters']['numeric_leaves'] == independent['numeric_leaves']
        for review in independent['nonfinite_reviews']:
            counts[review['kind']] += 1
            if review_sink is not None:
                review_sink(dict(iteration=row['iter'], **review))
        literals = primary['context']['literal_null_paths'] or primary['parameters']['literal_null_paths']
        if literals:
            literal_rows += 1
            if review_sink is not None:
                review_sink(dict(iteration=row['iter'], kind='literal_null_review'))
    finite = not counts and literal_rows == 0
    mapped = 0
    with (directory / 'C1.log').open() as handle:
        tsv = csv.DictReader(handle, delimiter='\t')
        assert tsv.fieldnames and len(tsv.fieldnames) == len(set(tsv.fieldnames))
        mapping = load((directory / 'C1.log.column-map.json').read_text())
        assert type(mapping) is dict and mapping
        assert all(type(k) is str and type(v) is str for k, v in mapping.items())
        assert set(mapping.values()) == set(tsv.fieldnames)
        missing = object()
        for index, (row, reference) in enumerate(zip_longest(scalar_records(path, iterations), tsv, fillvalue=missing)):
            assert row is not missing and reference is not missing, 'JSON/TSV row census mismatch'
            assert set(reference) == set(tsv.fieldnames) and all(v is not None for v in reference.values())
            assert reference['iter'] == str(index)
            scores = [float(reference[k]) for k in ('prior', 'likelihood', 'posterior')]
            assert all(math.isfinite(v) for v in scores)
            assert abs(scores[0] + scores[1] - scores[2]) < 1e-7
            if finite:
                flat = flattened(row)
                assert set(flat) == set(mapping)
                for field, column in mapping.items():
                    value, expected = flat[field], float(reference[column])
                    assert type(value) in (int, float) and math.isfinite(expected)
                    assert math.isclose(value, expected, rel_tol=2e-13, abs_tol=0), field
                    mapped += 1
    return dict(schema=SCHEMA, rows=iterations + 1, nonfinite_review_counts=dict(counts),
        nonfinite_reviews=sum(counts.values()), literal_null_rows=literal_rows, mapped_values_compared=mapped,
        fully_finite_mapped_integrity_checked=finite, scientific_eligibility=False, posterior_qualified=False)


def latent_trace(directory, iterations, prior, overflow_sink=None):
    directory = Path(directory)
    mu, scale = {'broad': (Decimal(0), Decimal(2)), 'centered': (Decimal(0), Decimal(1)),
                 'package': (Decimal(6), Decimal(2))}[prior]
    missing, count, overflows = object(), 0, 0
    paired = zip_longest(scalar_records(directory / 'C1.P1.log-alpha-samples.jsonl', iterations, decimal_load),
                        scalar_records(directory / 'C1.log.json', iterations, decimal_load), fillvalue=missing)
    for row, reference in paired:
        assert row is not missing and reference is not missing
        overflow = latent_check(row, mu, scale, reference)
        count += 1
        overflows += int(overflow)
        if overflow and overflow_sink is not None:
            state = row['parameters//']['S1/']
            overflow_sink(dict(iteration=row['iter'], finite_latent_log_alpha=str(state['latentLogAlpha']),
                latent_log_density=str(state['latentLogDensity']), derived_alpha=state['derivedAlpha']))
    assert count == iterations + 1
    return dict(rows=count, overflow_rows=overflows, decimal_precision=90,
                scientific_eligibility=False, posterior_qualified=False)


def frame_context(chain, directory, mapping):
    for name in ('alignment', 'tree'):
        assert sha(chain[name]) == chain[name + '_sha256']
    matched = match_trees(Path(chain['tree']).read_text(), (Path(directory) / 'runtime-tree.nwk').read_text())
    observed = {name: seq.replace('-', '') for name, seq in fasta_records(Path(chain['alignment']).read_text().splitlines()).items()}
    assert sorted(observed) == matched['tips'] and len(observed) == chain['proteins']
    bits = {tip: 1 << i for i, tip in enumerate(matched['tips'])}
    candidates = {}
    for row in mapping_rows(mapping, chain):
        tips = strict_json(row['retained_set_json'])
        node = row['source_node']
        assert tips and len(tips) == len(set(tips)) and set(tips) <= set(bits)
        mask = sum(bits[tip] for tip in tips)
        assert node not in candidates and matched['source_labels'][node]['mask'] == mask
        assert (int(row['level']) == 3) == (mask == (1 << len(bits)) - 1)
        if 'retained_descendants' in row:
            assert int(row['retained_descendants']) == len(tips)
        candidates[node] = matched['runtime_index'][mask]['label']
    assert len(candidates) == len(set(candidates.values())) == 4
    return dict(observed=observed, runtime_labels=matched['runtime_labels'], candidates=candidates,
                negative_branch_pairs_require_review=matched['negative_branch_pairs_require_review'])


def joint_frames(directory, context, iterations, interval=10):
    """Yield one decoded frame; exhaustion is required to establish completeness."""
    expected = horizon(iterations, interval)
    count = 0
    with (Path(directory) / 'C1.P1.site-property-samples.jsonl').open() as handle:
        for index, line in enumerate(handle):
            assert index < len(expected), 'Extra joint frame'
            frame = strict_json(line)
            yield decode(frame, expected[index], context['observed'], context['runtime_labels'], context['candidates'])
            count += 1
    assert count == len(expected), 'Missing joint frame'


def fasta_trace(directory, context, iterations, interval=10):
    expected = horizon(iterations, interval)
    count = residues = 0
    with (Path(directory) / 'C1.P1.fastas').open() as handle:
        for index, (iteration, sequences) in enumerate(native_alignments(handle)):
            assert index < len(expected) and iteration == expected[index], 'FASTA frame horizon mismatch'
            assert set(sequences) == set(context['runtime_labels'])
            assert len({len(seq) for seq in sequences.values()}) == 1
            for tip, source in context['observed'].items():
                observed = sequences[tip].replace('-', '')
                assert len(observed) == len(source)
                assert all(a == b or a == 'X' for a, b in zip(source, observed))
            for label in context['candidates'].values():
                residues += len(sequences[label].replace('-', ''))
            count += 1
    assert count == len(expected), 'Missing FASTA frame'
    return dict(saved_alignments=count, candidate_frames=count * len(context['candidates']),
                candidate_residues=residues, scientific_eligibility=False, posterior_qualified=False)


def verify_existing_exports(directory, context, iterations, records, interval=10, summary_sink=None):
    assert len(records) == len(horizon(iterations, interval))
    count = tip_pairs = ancestral_pairs = 0
    for index, (summary, arrays) in enumerate(joint_frames(directory, context, iterations, interval)):
        record = records[index]
        assert all(record[k] == v for k, v in summary.items()), 'Changed native frame summary'
        path = Path(record['projection_array'])
        assert sha(path) == record['projection_array_sha256']
        verify_arrays(path, arrays)
        if summary_sink is not None:
            summary_sink(summary)
        count += 1
        tip_pairs += summary['tip_pairs']
        ancestral_pairs += summary['ancestral_pairs']
    assert count == len(records)
    return dict(saved_joint_frames=count, joint_tip_pairs=tip_pairs, joint_ancestral_pairs=ancestral_pairs,
                scientific_eligibility=False, posterior_qualified=False)


def stage_exports(directory, context, iterations, destination, prior, interval=10, require_latent=True):
    """Transactionally stage software-checked frames; no biological admission.

    A failing frame leaves its pending namespace and failure record untouched.
    No final directory can appear until scalar, FASTA and every joint frame pass.
    Native custody and resource/source-change guards are required separately.
    """
    destination = Path(destination)
    pending = destination.with_name(destination.name + '.pending')
    assert not destination.exists() and not pending.exists()
    names = ['C1.log', 'C1.log.json', 'C1.log.column-map.json', 'C1.P1.fastas',
             'C1.P1.site-property-samples.jsonl', 'runtime-tree.nwk']
    if require_latent:
        names.append('C1.P1.log-alpha-samples.jsonl')
    source_hashes = {str(Path(directory) / name): sha(Path(directory) / name) for name in names}
    scalar = scalar_trace(directory, iterations)
    latent = latent_trace(directory, iterations, prior) if require_latent else None
    if not scalar['fully_finite_mapped_integrity_checked']:
        assert all(sha(p) == digest for p, digest in source_hashes.items())
        return dict(status='numeric_review_no_exports_created', scalar=scalar, latent=latent,
                    source_hashes=source_hashes, scientific_eligibility=False, posterior_qualified=False)
    fasta = fasta_trace(directory, context, iterations, interval)
    pending.mkdir(parents=True, exist_ok=False)
    count = 0
    try:
        with (pending / 'frames.jsonl').open('x') as ledger:
            for summary, arrays in joint_frames(directory, context, iterations, interval):
                name = 'frame-' + str(summary['iteration']) + '.npz'
                path = pending / name
                write_arrays(path, arrays)
                verify_arrays(path, arrays)
                ledger.write(json.dumps(dict(summary, projection_array=name, projection_array_sha256=sha(path))) + '\n')
                count += 1
        assert count == len(horizon(iterations, interval))
        assert all(sha(p) == digest for p, digest in source_hashes.items()), 'Source changed during staging'
        result = dict(status='complete_horizon_output_integrity_not_posterior', saved_joint_frames=count,
            scalar=scalar, latent=latent, fasta=fasta, frame_ledger_sha256=sha(pending / 'frames.jsonl'),
            source_hashes=source_hashes, scientific_eligibility=False, posterior_qualified=False,
            native_custody_qualified_here=False)
        with (pending / 'complete.json').open('x') as handle:
            json.dump(result, handle, indent=2)
            handle.write('\n')
        assert not destination.exists()
        os.rename(pending, destination)
        return result
    except Exception as error:
        with (pending / 'failure.json').open('x') as handle:
            json.dump(dict(status='failed_pending_exports_retained', complete_joint_frames_written=count,
                error_type=type(error).__name__, error=str(error), scientific_eligibility=False,
                posterior_qualified=False), handle, indent=2)
            handle.write('\n')
        raise
