#!/usr/bin/env python3
"""Retain full expanded pair and terminal-event denominators under six coverage screens."""
import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from screen_duplication_alignment_reuse import sha
from screen_background_whole_protein_coverage import reasons


def event_key(row):
    return tuple(row[k] for k in ['guide', 'family', 'taxon_id', 'gene_node']) + tuple(sorted([row['gene_a'], row['gene_b']]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    pins = {str(args.plan): sha(args.plan), **plan['pins']}

    def verify():
        for path, digest in pins.items():
            if sha(path) != digest:
                raise ValueError('Changed coverage source: ' + path)

    verify()
    queue = Path(plan['queue'])
    qr = json.loads((queue / 'receipt.json').read_text())
    summary_root = Path(plan['summary'])
    sr = json.loads((summary_root / 'receipt.json').read_text())
    audit = json.loads(Path(plan['summary_audit']).read_text())
    assert audit['status'] == 'passed_full_primary_usable_order_summary_readback'
    assert audit['source_receipt_sha256'] == sha(summary_root / 'receipt.json')
    assert qr['status'] == 'complete_reviewed_duplication_model_pair_queue'
    for root, receipt in [(queue, qr), (summary_root, sr)]:
        for name, digest in receipt['artifacts'].items():
            assert sha(root / name) == digest
    screens = json.loads(Path(plan['screens_plan']).read_text())['screens']
    assert screens == plan['screens']
    models = {}
    for line in (queue / 'models.jsonl').open():
        row = json.loads(line)
        key = row['model_id'], row['version']
        assert key not in models and row['length'] > 0
        models[key] = row['length']
    assert len(models) == qr['unique_models'] == plan['expected']['models']
    pairs = {}
    for row in csv.DictReader((queue / 'model_pairs.tsv').open(), delimiter='\t'):
        endpoints = [(row['model_' + s], int(row['version_' + s])) for s in ['a', 'b']]
        assert endpoints == sorted(endpoints) and endpoints[0] != endpoints[1]
        assert row['pair_key'] == hashlib.sha256(json.dumps(endpoints, separators=(',', ':')).encode()).hexdigest()
        assert row['pair_key'] not in pairs and all(x in models for x in endpoints)
        pairs[row['pair_key']] = dict(row, length_a=models[endpoints[0]], length_b=models[endpoints[1]])
    assert len(pairs) == qr['unique_distinct_model_pairs'] == plan['expected']['pairs']
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    screen_fields = [s['id'] + suffix for s in screens for suffix in ['_pass', '_exclusions']]
    fields = ['pair_key', 'mask', 'model_a', 'version_a', 'model_b', 'version_b', 'length_a', 'length_b']
    for order in [0, 1]:
        fields += [f'order{order}_status', f'order{order}_native_status', f'order{order}_numerical_exclusion_reasons',
                   f'order{order}_aligned_length', f'order{order}_original_coverage_a', f'order{order}_original_coverage_b']
    fields += screen_fields
    flags, pair_counts, pair_exclusions = {}, Counter(), Counter()
    with (out / 'pair_mask_coverage.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for original in csv.DictReader((summary_root / 'pair_mask_order_summary.tsv').open(), delimiter='\t'):
            key = original['pair_key'], original['mask']
            assert key not in flags and key[1] in ['full', 'plddt70']
            pair = pairs[key[0]]
            work = dict(pair, mask=key[1], measurement_disposition='new_model_pair',
                        **{'fit_' + k: v for k, v in original.items()})
            row = {k: work[k] for k in fields if k in work}
            for order in [0, 1]:
                for name in ['status', 'native_status', 'numerical_exclusion_reasons', 'aligned_length']:
                    row[f'order{order}_{name}'] = original[f'order{order}_{name}']
                for side in ['a', 'b']:
                    row[f'order{order}_original_coverage_{side}'] = (
                        int(float(original[f'order{order}_aligned_length'])) / pair['length_' + side]
                        if original[f'order{order}_status'] == 'aligned' else '')
            flags[key] = {}
            for spec in screens:
                name = spec['id']
                why = reasons(work, spec)
                flags[key][name] = why
                row[name + '_pass'] = int(not why)
                row[name + '_exclusions'] = ';'.join(why)
                pair_counts[key[1] + ':' + name] += not why
                for reason in why:
                    pair_exclusions[key[1] + ':' + name + ':' + reason] += 1
            writer.writerow(row)
    assert len(flags) == audit['pair_mask_rows'] == plan['expected']['pair_mask_rows']
    assert set(flags) == {(key, mask) for key in pairs for mask in ['full', 'plddt70']}
    pair_joint = {s['id']: sum(not flags[key, 'full'][s['id']] and not flags[key, 'plddt70'][s['id']]
                             for key in pairs) for s in screens}
    links = {}
    for row in csv.DictReader((queue / 'event_model_pair_links.tsv').open(), delimiter='\t'):
        key = event_key(row)
        assert key not in links
        links[key] = row
    assert len(links) == qr['counts']['event_links'] == plan['expected']['modeled_events']
    taxa = {r['taxon_id']: r for r in csv.DictReader(Path(plan['manifest']).open(), delimiter='\t')}
    assert len(taxa) == plan['expected']['taxa']
    seen_links, seen_events = set(), set()
    event_counts, event_exclusions, dispositions, event_joint = Counter(), Counter(), Counter(), Counter()
    taxon_counts = defaultdict(Counter)
    taxon_pass = Counter()
    with gzip.open(plan['events'], 'rt') as source, gzip.open(out / 'event_mask_coverage.tsv.gz', 'wt', compresslevel=1) as handle:
        reader = csv.DictReader(source, delimiter='\t')
        identity_fields = reader.fieldnames
        fields = identity_fields + ['mask', 'pair_key', 'version_a', 'version_b', 'measurement_disposition'] + screen_fields
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for ix, original in enumerate(reader, 1):
            key = event_key(original)
            assert key not in seen_events and original['taxon_id'] in taxa
            seen_events.add(key)
            link = links.get(key)
            covered = original['new_coverage_class'] == 'both_models'
            assert bool(link) == covered
            if link:
                seen_links.add(key)
                assert all(original[k] == link[k] for k in ['guide', 'family', 'taxon_id', 'gene_node'])
                gene_models = {link['gene_' + side]: (link['model_' + side], link['version_' + side]) for side in ['a', 'b']}
                assert set(gene_models) == {original['gene_a'], original['gene_b']}
                assert all(original['new_model_' + side] == gene_models[original['gene_' + side]][0] for side in ['a', 'b'])
                assert int(original['new_same_model']) == int(link['same_model'])
                disposition = ('new_model_pair' if link['comparison_status'] == 'queued_distinct_models'
                               else link['comparison_status'])
                assert disposition in ['new_model_pair', 'identical_model_no_alignment', 'unresolved_tree']
            else:
                assert original['new_coverage_class'] in ['one_model', 'neither_model']
                disposition = 'one_model' if original['new_coverage_class'] == 'one_model' else 'no_models'
            guide, taxon = original['guide'], original['taxon_id']
            group = guide, taxon
            taxon_counts[group]['events'] += 1
            taxon_counts[group]['both_models'] += covered
            taxon_counts[group]['distinct_model_pair'] += disposition == 'new_model_pair'
            passed = {}
            for mask in ['full', 'plddt70']:
                row = dict(original, mask=mask, pair_key=link['pair_key'] if link else '',
                           version_a=gene_models[original['gene_a']][1] if link else '', version_b=gene_models[original['gene_b']][1] if link else '',
                           measurement_disposition=disposition)
                dispositions[guide + ':' + mask + ':' + disposition] += 1
                passed[mask] = {}
                for spec in screens:
                    name = spec['id']
                    why = flags[link['pair_key'], mask][name] if disposition == 'new_model_pair' else [disposition]
                    row[name + '_pass'], row[name + '_exclusions'] = int(not why), ';'.join(why)
                    passed[mask][name] = not why
                    event_counts[guide + ':' + mask + ':' + name] += not why
                    taxon_pass[guide, taxon, mask, name] += not why
                    for reason in why:
                        event_exclusions[guide + ':' + mask + ':' + name + ':' + reason] += 1
                writer.writerow(row)
            for spec in screens:
                event_joint[guide + ':' + spec['id']] += passed['full'][spec['id']] and passed['plddt70'][spec['id']]
            if ix % 100000 == 0:
                print('Full event coverage rows', ix, '/', plan['expected']['events'], flush=True)
    assert len(seen_events) == plan['expected']['events'] and seen_links == set(links)
    taxon_path = out / 'taxon_screen_coverage.tsv'
    with taxon_path.open('w') as handle:
        fields = ['guide', 'taxon_id', 'species_name', 'study_role', 'lineage', 'mask', 'screen',
                  'events', 'both_models', 'distinct_model_pair', 'passed_events']
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for guide in ['profile', 'mafft']:
            for taxon, record in sorted(taxa.items()):
                for mask in ['full', 'plddt70']:
                    for spec in screens:
                        count = taxon_counts[guide, taxon]
                        writer.writerow(dict(guide=guide, taxon_id=taxon, **{k: record[k] for k in ['species_name', 'study_role', 'lineage']},
                                             mask=mask, screen=spec['id'], **{k: count[k] for k in ['events', 'both_models', 'distinct_model_pair']},
                                             passed_events=taxon_pass[guide, taxon, mask, spec['id']]))
    verify()
    result = dict(status='complete_full_expanded_duplication_coverage_pending_independent_readback',
                  plan_sha256=sha(args.plan), source_hashes=pins, screens=screens,
                  pairs=len(pairs), pair_mask_rows=len(flags), events=len(seen_events), event_mask_rows=2 * len(seen_events),
                  modeled_events=len(links), taxon_screen_rows=2 * len(taxa) * 2 * len(screens),
                  pair_pass_counts=dict(pair_counts), pair_both_masks_pass_counts=pair_joint,
                  pair_exclusion_counts=dict(pair_exclusions), event_pass_counts=dict(event_counts),
                  event_both_masks_pass_counts=dict(event_joint), event_exclusion_counts=dict(event_exclusions),
                  event_dispositions=dict(dispositions), artifacts={p.name: sha(p) for p in out.iterdir()}, scientific_eligibility=False,
                  scope='All distinct model pairs and the complete frozen terminal singleton-side event ledger under both masks and six unchanged coverage screens. '
                        'Both input orders required; confidence-masked coverage uses original full protein lengths. '
                        'All missing/same-model/unresolved and numerical/coverage exclusions retained; all 526 manifest entries receive explicit zero-event rows where needed. '
                        'These are overlapping dependent candidate events under two guides, not all reconciled duplications, confidence calibration, domain orientation validation, matched controls or biological effect inference.')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'event_exclusion_counts', 'pair_exclusion_counts']}, indent=2), flush=True)


if __name__ == '__main__':
    main()
