#!/usr/bin/env python3
"""Replay every expanded pair/event coverage screen with decimal ceilings and SQL joins."""
import argparse
import csv
import gzip
import json
import math
import sqlite3
from collections import Counter, defaultdict
from decimal import Decimal, ROUND_CEILING
from itertools import zip_longest
from pathlib import Path

from screen_duplication_alignment_reuse import sha


def exclusions(row, lengths, spec):
    why, aligned = [], []
    for order in [0, 1]:
        if row[f'order{order}_status'] != 'aligned':
            why.append(f'order{order}_not_numerically_usable')
        else:
            n = Decimal(row[f'order{order}_aligned_length'])
            assert n == n.to_integral_value() and 0 < n <= min(lengths)
            aligned.append(int(n))
    if any(n < spec['minimum_aligned_residues'] for n in aligned):
        why.append('short_alignment')
    minimum = max(int((Decimal(length) * Decimal(str(spec['minimum_original_coverage']))).to_integral_value(rounding=ROUND_CEILING))
                  for length in lengths)
    if any(n < minimum for n in aligned):
        why.append('low_original_protein_coverage')
    return why


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    plan = json.loads(args.plan.read_text())
    root, queue = Path(plan['output']), Path(plan['queue'])
    r = json.loads((root / 'receipt.json').read_text())
    assert r['status'] == 'complete_full_expanded_duplication_coverage_pending_independent_readback'
    assert r['plan_sha256'] == sha(args.plan)
    bindings = {str(args.plan): sha(args.plan), **plan['pins'], str(root / 'receipt.json'): sha(root / 'receipt.json')}
    for name, digest in r['artifacts'].items():
        bindings[str(root / name)] = digest

    def verify():
        for path, digest in bindings.items():
            if sha(path) != digest:
                raise ValueError('Changed coverage readback source: ' + path)

    verify()
    screens = json.loads(Path(plan['screens_plan']).read_text())['screens']
    assert screens == r['screens'] == plan['screens']
    db = sqlite3.connect(':memory:')
    db.execute('CREATE TABLE models (model TEXT, version INTEGER, length INTEGER, PRIMARY KEY (model,version))')
    for line in (queue / 'models.jsonl').open():
        m = json.loads(line)
        db.execute('INSERT INTO models VALUES (?,?,?)', (m['model_id'], m['version'], m['length']))
    assert db.execute('SELECT COUNT(*) FROM models').fetchone()[0] == plan['expected']['models']
    db.execute('CREATE TABLE pairs (pair TEXT PRIMARY KEY, a TEXT, av INTEGER, b TEXT, bv INTEGER)')
    db.executemany('INSERT INTO pairs VALUES (?,?,?,?,?)',
                   ((x['pair_key'], x['model_a'], int(x['version_a']), x['model_b'], int(x['version_b']))
                    for x in csv.DictReader((queue / 'model_pairs.tsv').open(), delimiter='\t')))
    lengths = {key: (a, av, b, bv, al, bl) for key, a, av, b, bv, al, bl in db.execute(
        'SELECT p.*,a.length,b.length FROM pairs p JOIN models a ON p.a=a.model AND p.av=a.version '
        'JOIN models b ON p.b=b.model AND p.bv=b.version')}
    assert len(lengths) == plan['expected']['pairs']
    numeric, flags = {}, {}
    pair_counts, pair_exclusions = Counter(), Counter()
    pair_decisions = 0
    for source, output in zip_longest(
            csv.DictReader((Path(plan['summary']) / 'pair_mask_order_summary.tsv').open(), delimiter='\t'),
            csv.DictReader((root / 'pair_mask_coverage.tsv').open(), delimiter='\t')):
        assert source is not None and output is not None
        key = source['pair_key'], source['mask']
        assert key not in numeric and (output['pair_key'], output['mask']) == key
        numeric[key] = source
        a, av, b, bv, al, bl = lengths[key[0]]
        for field, expected in [('model_a', a), ('version_a', av), ('model_b', b), ('version_b', bv), ('length_a', al), ('length_b', bl)]:
            assert output[field] == str(expected)
        for order in [0, 1]:
            for field in ['status', 'native_status', 'numerical_exclusion_reasons', 'aligned_length']:
                assert source[f'order{order}_{field}'] == output[f'order{order}_{field}']
            for side, length in [('a', al), ('b', bl)]:
                observed = output[f'order{order}_original_coverage_{side}']
                if source[f'order{order}_status'] == 'aligned':
                    expected = float(Decimal(source[f'order{order}_aligned_length']) / Decimal(length))
                    assert math.isclose(float(observed), expected, rel_tol=1e-15, abs_tol=1e-15)
                else:
                    assert observed == ''
        flags[key] = {}
        for spec in screens:
            name = spec['id']
            why = exclusions(source, [al, bl], spec)
            assert output[name + '_exclusions'] == ';'.join(why) and int(output[name + '_pass']) == int(not why)
            flags[key][name] = why
            pair_counts[key[1] + ':' + name] += not why
            for reason in why:
                pair_exclusions[key[1] + ':' + name + ':' + reason] += 1
            pair_decisions += 1
    assert set(numeric) == {(key, mask) for key in lengths for mask in ['full', 'plddt70']}
    assert len(numeric) == r['pair_mask_rows'] == plan['expected']['pair_mask_rows']
    joint = {s['id']: sum(not flags[key, 'full'][s['id']] and not flags[key, 'plddt70'][s['id']]
                         for key in lengths) for s in screens}
    assert dict(pair_counts) == r['pair_pass_counts'] and dict(pair_exclusions) == r['pair_exclusion_counts']
    assert joint == r['pair_both_masks_pass_counts']
    db.execute('CREATE TABLE links (guide TEXT,family TEXT,taxon TEXT,node TEXT,lo TEXT,hi TEXT,payload TEXT, '
               'PRIMARY KEY (guide,family,taxon,node,lo,hi))')
    for x in csv.DictReader((queue / 'event_model_pair_links.tsv').open(), delimiter='\t'):
        db.execute('INSERT INTO links VALUES (?,?,?,?,MIN(?,?),MAX(?,?),?)',
                   (x['guide'], x['family'], x['taxon_id'], x['gene_node'], x['gene_a'], x['gene_b'],
                    x['gene_a'], x['gene_b'], json.dumps(x)))
    assert db.execute('SELECT COUNT(*) FROM links').fetchone()[0] == plan['expected']['modeled_events']
    taxa = {x['taxon_id']: x for x in csv.DictReader(Path(plan['manifest']).open(), delimiter='\t')}
    assert len(taxa) == plan['expected']['taxa']
    counts, excluded, dispositions, both = Counter(), Counter(), Counter(), Counter()
    taxon_counts, taxon_pass = defaultdict(Counter), Counter()
    seen_events, seen_links = set(), set()
    event_decisions = 0
    with gzip.open(plan['events'], 'rt') as source, gzip.open(root / 'event_mask_coverage.tsv.gz', 'rt') as output:
        original = csv.DictReader(source, delimiter='\t')
        actual = csv.DictReader(output, delimiter='\t')
        for ix, event in enumerate(original, 1):
            prefix = tuple(event[k] for k in ['guide', 'family', 'taxon_id', 'gene_node'])
            gene_set = frozenset([event['gene_a'], event['gene_b']])
            identity = prefix, gene_set
            assert identity not in seen_events and event['taxon_id'] in taxa
            seen_events.add(identity)
            q = db.execute('SELECT payload FROM links WHERE guide=? AND family=? AND taxon=? AND node=? '
                           'AND lo=MIN(?,?) AND hi=MAX(?,?)', (*prefix, event['gene_a'], event['gene_b'], event['gene_a'], event['gene_b'])).fetchone()
            covered = event['new_coverage_class'] == 'both_models'
            assert bool(q) == covered
            if q:
                seen_links.add(identity)
                link = json.loads(q[0])
                for field in ['guide', 'family', 'taxon_id', 'gene_node']:
                    assert event[field] == link[field]
                mapped_versions = {}
                for side in ['a', 'b']:
                    matches = [qside for qside in ['a', 'b'] if event['gene_' + side] == link['gene_' + qside]]
                    assert len(matches) == 1
                    qside = matches[0]
                    assert event['new_model_' + side] == link['model_' + qside]
                    mapped_versions[side] = link['version_' + qside]
                assert int(event['new_same_model']) == int(link['same_model'])
                disposition = {'queued_distinct_models': 'new_model_pair',
                               'identical_model_no_alignment': 'identical_model_no_alignment',
                               'unresolved_tree': 'unresolved_tree'}[link['comparison_status']]
            else:
                disposition = {'one_model': 'one_model', 'neither_model': 'no_models'}[event['new_coverage_class']]
                link = None
            guide, taxon = event['guide'], event['taxon_id']
            group = guide, taxon
            taxon_counts[group]['events'] += 1
            taxon_counts[group]['both_models'] += covered
            taxon_counts[group]['distinct_model_pair'] += disposition == 'new_model_pair'
            passed = {}
            for mask in ['full', 'plddt70']:
                row = next(actual)
                assert all(row[k] == value for k, value in event.items())
                assert row['mask'] == mask and row['measurement_disposition'] == disposition
                assert row['pair_key'] == (link['pair_key'] if link else '')
                for side in ['a', 'b']:
                    assert row['version_' + side] == (mapped_versions[side] if link else '')
                dispositions[guide + ':' + mask + ':' + disposition] += 1
                passed[mask] = {}
                for spec in screens:
                    name = spec['id']
                    why = flags[link['pair_key'], mask][name] if disposition == 'new_model_pair' else [disposition]
                    assert row[name + '_exclusions'] == ';'.join(why) and int(row[name + '_pass']) == int(not why)
                    passed[mask][name] = not why
                    counts[guide + ':' + mask + ':' + name] += not why
                    taxon_pass[guide, taxon, mask, name] += not why
                    for reason in why:
                        excluded[guide + ':' + mask + ':' + name + ':' + reason] += 1
                    event_decisions += 1
            for spec in screens:
                both[guide + ':' + spec['id']] += passed['full'][spec['id']] and passed['plddt70'][spec['id']]
            if ix % 100000 == 0:
                print('Independently checked full event rows', ix, '/', plan['expected']['events'], flush=True)
        assert next(actual, None) is None
    assert len(seen_events) == r['events'] == plan['expected']['events'] and len(seen_links) == plan['expected']['modeled_events']
    assert dict(counts) == r['event_pass_counts'] and dict(excluded) == r['event_exclusion_counts']
    assert dict(dispositions) == r['event_dispositions'] and dict(both) == r['event_both_masks_pass_counts']
    seen = set()
    for row in csv.DictReader((root / 'taxon_screen_coverage.tsv').open(), delimiter='\t'):
        key = tuple(row[k] for k in ['guide', 'taxon_id', 'mask', 'screen'])
        assert key not in seen
        seen.add(key)
        guide, taxon, mask, name = key
        assert taxon in taxa
        for field in ['species_name', 'study_role', 'lineage']:
            assert row[field] == taxa[taxon][field]
        for field in ['events', 'both_models', 'distinct_model_pair']:
            assert int(row[field]) == taxon_counts[guide, taxon][field]
        assert int(row['passed_events']) == taxon_pass[key]
    assert seen == {(g, t, m, s['id']) for g in ['profile', 'mafft'] for t in taxa
                    for m in ['full', 'plddt70'] for s in screens}
    assert len(seen) == r['taxon_screen_rows']
    verify()
    result = dict(status='passed_full_expanded_pair_event_and_taxon_coverage_readback',
                  producer_receipt_sha256=sha(root / 'receipt.json'), plan_sha256=sha(args.plan),
                  checker_sha256=sha(__file__), pair_mask_rows=len(numeric), pair_screen_decisions_checked=pair_decisions,
                  event_rows=len(seen_events), event_mask_rows=2 * len(seen_events),
                  event_screen_decisions_checked=event_decisions, taxon_screen_rows=len(seen),
                  pair_pass_counts=dict(pair_counts), pair_both_masks_pass_counts=joint,
                  event_pass_counts=dict(counts), event_both_masks_pass_counts=dict(both), source_hashes=bindings,
                  scientific_eligibility=False,
                  scope='Every pair screen independently recomputed with decimal-ceiling original-length thresholds and SQL model/length joins. '
                        'Every source event field, masked disposition, screen and exclusion retained and checked against SQL queue joins and independently recomputed pair decisions. '
                        'Every aggregate, both-mask intersection and all manifest taxon/guide/mask/screen cells, including zero-event cells, verified. '
                        'No predictor calibration, matched controls, species independence or biological effect acceptance.')
    with args.output.open('x') as handle:
        handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True)


if __name__ == '__main__':
    main()
