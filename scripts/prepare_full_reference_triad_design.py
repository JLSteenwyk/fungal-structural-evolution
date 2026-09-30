#!/usr/bin/env python3
"""Preserve every context/tie while designing all three duplicate/reference edges."""
import argparse
import csv
import gzip
import hashlib
import json
import shutil
from collections import Counter, defaultdict
from pathlib import Path
from reference_triad_design_sources import load_sources
from reference_measurement_union_sources import bind, verify, pair_ends
from run_ortholog_pair_guide_comparison import sha


def digest(value):
    return hashlib.sha256(json.dumps(value, separators=(',', ':')).encode()).hexdigest()


def edge(start, end, roles, pairs, kind):
    key = digest(sorted([start, end])); found = pairs.get(key)
    if start == end: disposition = 'identical_model_not_independent'; endpoints = order = None; work = None
    elif found:
        ends = pair_ends(found); assert sorted(ends) == sorted([start, end])
        disposition = 'in_full_' + kind + '_measurement_design'; endpoints = [list(e) for e in ends]
        order = 0 if ends == [start, end] else 1
        work = found.get('work_disposition', 'primary_pair')
    else: disposition = 'outside_full_' + kind + '_measurement_design'; endpoints = order = work = None
    return dict(pair_key=key, measurement_disposition=disposition, desired_role_order=roles,
                current_pair_endpoints=endpoints, desired_direction_order=order, pair_work_disposition=work)


def physical(models, primary, references):
    a, b, r = models
    edges = dict(ab=edge(a, b, ['a', 'b'], primary, 'primary'),
                 ar=edge(r, a, ['reference', 'a'], references, 'reference'),
                 br=edge(r, b, ['reference', 'b'], references, 'reference'))
    versioned = len(set(models)); identities = len({m[0] for m in models}); reasons = []
    for label, left, right in [('duplicate_model_identity', a, b), ('a_reference_model_identity', a, r), ('b_reference_model_identity', b, r)]:
        if left == right: reasons.append(label)
    if identities < 3: reasons.append('shared_model_id_across_roles')
    for label, value in edges.items():
        if value['measurement_disposition'].startswith('outside_'): reasons.append(label + '_' + value['measurement_disposition'])
    all_pairs = edges['ab']['measurement_disposition'] == 'in_full_primary_measurement_design' and all(edges[s]['measurement_disposition'] == 'in_full_reference_measurement_design' for s in ['ar', 'br'])
    return dict(triad_id=digest(models), role_order=['a', 'b', 'reference'], models=[list(m) for m in models],
                distinct_versioned_models=versioned, distinct_model_ids=identities, edges=edges,
                all_three_pairs_in_measurement_design=all_pairs, model_work_exclusion_reasons=sorted(reasons))


def link(triad, reference, context, queued):
    parent = context['parent_context_eligible']; reasons = list(triad['model_work_exclusion_reasons']) if triad else ['reference_model_missing']
    if not parent: reasons.append('parent_context_excluded')
    if queued != 'queued_distinct_models': reasons.append('duplicate_comparison_' + queued)
    ready = bool(triad and triad['distinct_versioned_models'] == 3 and triad['distinct_model_ids'] == 3 and triad['all_three_pairs_in_measurement_design'] and parent and queued == 'queued_distinct_models')
    return dict(reference_gene=reference['reference_gene'], triad_id=triad['triad_id'] if triad else None,
                three_distinct_versioned_models=bool(triad and triad['distinct_versioned_models'] == 3),
                three_distinct_model_ids=bool(triad and triad['distinct_model_ids'] == 3),
                all_three_pairs_in_measurement_design=bool(triad and triad['all_three_pairs_in_measurement_design']),
                source_design_ready_for_correspondence=ready,
                source_design_ready_native_own=ready and reference['native_coorthology'][context['source_guide']] == 'both',
                source_design_ready_native_both_guides=ready and all(reference['native_coorthology'][g] == 'both' for g in ['profile', 'mafft']),
                work_exclusion_reasons=sorted(set(reasons)))


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True)
    args = p.parse_args(); plan = json.loads(args.plan.read_text())
    if 'resources' in plan: assert shutil.disk_usage('.').free >= plan['resources']['minimum_free_disk_gib'] * 2 ** 30
    source, queue, inventory, upstream, bindings = load_sources(plan); bind(bindings, args.plan)
    catalogs = []
    for root, n in [(queue, plan['expected']['primary_pairs']), (inventory, plan['expected']['reference_pairs'])]:
        with (root / 'model_pairs.tsv').open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        pairs = {r['pair_key']: r for r in rows}; assert len(pairs) == len(rows) == n
        for key, value in pairs.items(): assert key == digest(sorted(pair_ends(value))) and len(set(pair_ends(value))) == 2
        catalogs.append(pairs)
    primary, references = catalogs; out = Path(plan['output']); out.mkdir(parents=True, exist_ok=False)
    triads, seen, identities = {}, set(), set(); groups, duplicates, empty, guides = defaultdict(Counter), Counter(), Counter(), Counter()
    contexts = ties = ledger_links = 0
    with source.open() as original, gzip.open(out / 'context_triad_design.jsonl.gz', 'xt', compresslevel=1) as target:
        for line in original:
            work = json.loads(line); context = work['native_context']; guide = context['source_guide']; s = context['source']
            ordinal = guide, context['source_row_number']; identity = guide, s['family'], s['gene_node'], *sorted([s['gene_a'], s['gene_b']])
            assert ordinal not in seen and identity not in identities; seen.add(ordinal); identities.add(identity)
            a, b = [(work['duplicate_models'][side]['model_id'], work['duplicate_models'][side]['version']) for side in ['a', 'b']]
            ab = edge(a, b, ['a', 'b'], primary, 'primary'); assert ab['pair_key'] == work['duplicate_pair_key']
            queued = work['duplicate_comparison_status']; designs = {}
            if queued == 'queued_distinct_models': assert ab['measurement_disposition'] == 'in_full_primary_measurement_design'
            duplicates[f'{guide}|{queued}|{ab["measurement_disposition"]}'] += 1
            for design, refs in work['measurement_designs'].items():
                projected = []; prefix = f'{guide}|{design}|parent={int(context["parent_context_eligible"])}'
                if not refs: empty[prefix] += 1
                for ref in refs:
                    reference = ref['reference']; r = (reference['reference_model'], int(reference['reference_version'])) if reference['reference_model'] else None
                    triad = physical([a, b, r], primary, references) if r else None
                    value = link(triad, reference, context, queued)
                    if triad:
                        for side, label in [('a', 'ar'), ('b', 'br')]:
                            actual = ref['side_work'][side]; e = triad['edges'][label]
                            assert e['pair_key'] == actual['pair_key'] and e['measurement_disposition'] == actual['measurement_disposition']
                            if e['current_pair_endpoints'] is not None:
                                assert e['current_pair_endpoints'][actual['current_pair_focal_endpoint']] == list([a, b][['a', 'b'].index(side)])
                        key = triad['triad_id']
                        if key not in triads: triads[key] = dict(**triad, logical_reference_links=0, source_design_ready_links=0)
                        assert all(triads[key][k] == v for k, v in triad.items())
                        triads[key]['logical_reference_links'] += 1
                        triads[key]['source_design_ready_links'] += int(value['source_design_ready_for_correspondence'])
                    projected.append(value); ties += 1
                    ledger_links += sum(int(x['availability_ledger_link']) for x in ref['side_work'].values())
                    g = groups[prefix]; g['reference_tie_records'] += 1
                    g['modeled_reference_links'] += int(triad is not None); g['unmodeled_reference_links'] += int(triad is None)
                    for field in ['three_distinct_versioned_models', 'three_distinct_model_ids', 'all_three_pairs_in_measurement_design', 'source_design_ready_for_correspondence', 'source_design_ready_native_own', 'source_design_ready_native_both_guides']:
                        g[field] += int(value[field])
                    g['lexical_reference_links'] += int(reference['lexical_choice'])
                    g['lexical_source_design_ready_native_both_guides'] += int(reference['lexical_choice'] and value['source_design_ready_native_both_guides'])
                designs[design] = projected
            target.write(json.dumps(dict(source_design=work, duplicate_pair_design=ab, triad_designs=designs), separators=(',', ':')) + '\n')
            guides[guide] += 1; contexts += 1
            if contexts % 25000 == 0: print('Full three-edge triad design', contexts, '/', plan['expected']['target_contexts'], flush=True)
    assert contexts == plan['expected']['target_contexts'] and dict(guides) == upstream['guide_contexts']
    assert seen == {(g, n) for g, count in guides.items() for n in range(1, count + 1)}
    assert ties == plan['expected']['reference_tie_records'] and 2 * ties == plan['expected']['duplicate_reference_links'] and ledger_links == plan['expected']['availability_side_links']
    with (out / 'ordered_model_triads.jsonl').open('x') as handle:
        for key in sorted(triads): handle.write(json.dumps(triads[key], separators=(',', ':')) + '\n')
    model_counts = Counter(f'versioned={t["distinct_versioned_models"]}|ids={t["distinct_model_ids"]}' for t in triads.values())
    ready = sum(t['source_design_ready_links'] > 0 for t in triads.values()); verify(bindings)
    result = dict(status='complete_full_reference_triad_work_design_pending_independent_readback', plan_sha256=sha(args.plan),
                  target_contexts=contexts, context_design_records=2 * contexts, reference_tie_records=ties, duplicate_reference_links=2 * ties,
                  availability_side_links_preserved=ledger_links, guide_contexts=dict(guides), unique_ordered_model_triads=len(triads), correspondence_work_triads=ready,
                  potential_correspondence_states=16 * ready, triad_model_identity_counts=dict(model_counts), duplicate_design_counts=dict(duplicates),
                  reference_link_counts={k: dict(v) for k, v in groups.items()}, empty_reference_context_counts=dict(empty), source_hashes=bindings,
                  artifacts={name: sha(out / name) for name in ['context_triad_design.jsonl.gz', 'ordered_model_triads.jsonl']}, scientific_eligibility=False,
                  scope=plan['scope'])
    with (out / 'receipt.json').open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'artifacts']}, indent=2), flush=True)


if __name__ == '__main__': main()
