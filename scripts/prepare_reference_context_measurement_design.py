#!/usr/bin/env python3
"""Attach full model/version work states to all native reference contexts and ties."""
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path
from reference_context_measurement_sources import load_sources, load_tables, event_key, pair_hash
from reference_measurement_union_sources import bind, verify, pair_ends
from run_ortholog_pair_guide_comparison import sha


def project(context, duplicates, pairs, links, seen_links):
    source = context['source']; guide = context['source_guide']; identity = event_key(guide, source)
    assert context['source_row_sha256'] == hashlib.sha256(json.dumps(source, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    assert context['parent_context_eligible'] == (source['status'] in ['provisional_reference_available', 'no_modeled_nonfocal_sister'])
    queued = duplicates[identity]
    for field in ['family', 'gene_node', 'taxon_id']: assert source[field] == queued[field]
    gene_models = {queued['gene_' + s]: (queued['model_' + s], int(queued['version_' + s])) for s in ['a', 'b']}
    assert set(gene_models) == {source['gene_a'], source['gene_b']}
    mapped = {side: dict(gene=source['gene_' + side], model_id=gene_models[source['gene_' + side]][0], version=gene_models[source['gene_' + side]][1]) for side in ['a', 'b']}
    assert all(mapped[s]['model_id'] == source['model_' + s] for s in ['a', 'b'])
    assert queued['pair_key'] == pair_hash(list(gene_models.values()))
    assert bool(int(queued['same_model'])) == (gene_models[source['gene_a']] == gene_models[source['gene_b']])
    designs = {}
    for design in ['availability', 'sequence_first']:
        refs = context['designs'][design]
        assert [r['reference_gene'] for r in refs] == sorted({r['reference_gene'] for r in refs})
        assert [r['lexical_choice'] for r in refs] == ([True] + [False] * (len(refs) - 1) if refs else [])
        projected = []
        for reference in refs:
            assert all(value in ['both', 'only_a', 'only_b', 'neither'] for value in reference['native_coorthology'].values())
            work = {}
            for side in ['a', 'b']:
                focal = gene_models[source['gene_' + side]]
                ref = (reference['reference_model'], int(reference['reference_version'])) if reference['reference_model'] else None
                if ref is None: assert reference['reference_version'] == ''
                pair = pair_hash([focal, ref]) if ref else ''
                if ref is None: disposition = 'reference_model_missing'; endpoint = ''; ledger_work = ''
                elif focal == ref: disposition = 'identical_model_not_independent'; endpoint = ''; ledger_work = ''
                elif pair in pairs:
                    ends = pair_ends(pairs[pair]); assert sorted(ends) == sorted([focal, ref])
                    disposition = 'in_full_reference_measurement_design'; endpoint = ends.index(focal); ledger_work = pairs[pair]['work_disposition']
                else: disposition = 'outside_full_reference_measurement_design'; endpoint = ''; ledger_work = ''
                key = (*identity, reference['reference_gene'], source['gene_' + side])
                availability_link = links.get(key) if design == 'availability' else None
                if availability_link is not None:
                    seen_links.add(key)
                    assert availability_link['focal_model'] == focal[0] and int(availability_link['focal_version']) == focal[1]
                    assert ref and availability_link['reference_model'] == ref[0] and int(availability_link['reference_version']) == ref[1]
                    assert availability_link['pair_key'] == pair and int(availability_link['lexical_representative']) == int(reference['lexical_choice'])
                    assert availability_link['work_disposition'] == ('identical_model' if focal == ref else ledger_work)
                if design == 'availability' and context['parent_context_eligible']:
                    assert availability_link is not None and disposition in ['in_full_reference_measurement_design', 'identical_model_not_independent']
                work[side] = dict(duplicate_gene=source['gene_' + side], duplicate_model=focal[0], duplicate_version=focal[1],
                                  reference_model=ref[0] if ref else '', reference_version=ref[1] if ref else '', pair_key=pair,
                                  measurement_disposition=disposition, current_pair_focal_endpoint=endpoint,
                                  pair_work_disposition=ledger_work, availability_ledger_link=bool(availability_link))
            projected.append(dict(reference=reference, side_work=work))
        designs[design] = projected
    return dict(native_context=context, duplicate_models=mapped, duplicate_pair_key=queued['pair_key'],
                duplicate_comparison_status=queued['comparison_status'], measurement_designs=designs)


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True)
    args = p.parse_args(); plan = json.loads(args.plan.read_text())
    source, queue, inventory, upstream, bindings = load_sources(plan); bind(bindings, args.plan)
    duplicates, pairs, links = load_tables(queue, inventory)
    assert len(duplicates) == plan['expected']['target_contexts'] and len(pairs) == plan['expected']['model_pairs'] and len(links) == plan['expected']['availability_side_links']
    out = Path(plan['output']); out.mkdir(parents=True, exist_ok=False)
    seen, seen_links, row_ids = set(), set(), set(); counts, guide_counts = Counter(), Counter(); ties = logical = 0
    with source.open() as original, (out / 'context_measurement_design.jsonl').open('x') as target:
        for line in original:
            context = json.loads(line); key = event_key(context['source_guide'], context['source'])
            ordinal = context['source_guide'], context['source_row_number']
            assert key not in seen and ordinal not in row_ids; seen.add(key); row_ids.add(ordinal)
            row = project(context, duplicates, pairs, links, seen_links); guide_counts[context['source_guide']] += 1
            for design, refs in row['measurement_designs'].items():
                if not refs: counts[f'{context["source_guide"]}|{design}|parent={int(context["parent_context_eligible"])}|no_reference_gene'] += 1
                for ref in refs:
                    ties += 1
                    for work in ref['side_work'].values():
                        logical += 1
                        counts[f'{context["source_guide"]}|{design}|parent={int(context["parent_context_eligible"])}|{work["measurement_disposition"]}'] += 1
            target.write(json.dumps(row, separators=(',', ':')) + '\n')
            if len(seen) % 25000 == 0: print('Full reference context/model design', len(seen), '/', len(duplicates), flush=True)
    assert seen == set(duplicates) and seen_links == set(links)
    assert dict(guide_counts) == upstream['guide_contexts']
    assert row_ids == {(g, n) for g, count in guide_counts.items() for n in range(1, count + 1)}
    assert ties == plan['expected']['reference_tie_records'] and logical == plan['expected']['duplicate_reference_links']
    verify(bindings)
    result = dict(status='complete_full_reference_context_measurement_design_pending_independent_readback',
                  plan_sha256=sha(args.plan), target_contexts=len(seen), context_design_records=2 * len(seen),
                  reference_tie_records=ties, duplicate_reference_links=logical, availability_side_links_checked=len(seen_links),
                  guide_contexts=dict(guide_counts), measurement_disposition_counts=dict(counts),
                  source_hashes=bindings, artifacts={'context_measurement_design.jsonl': sha(out / 'context_measurement_design.jsonl')}, scientific_eligibility=False,
                  scope='All original closed native contexts/fields, both reference designs, all ties, lexical choices, native memberships and parent exclusions retained. Duplicate models/versions joined by actual gene identity rather than table-side position. Every availability-side ledger link checked; every distinct endpoint pair mapped to current full reference ledger or explicitly outside it. Missing reference structures and identical models remain separate states; parent exclusions are never overridden by measurements. This fixed source-only design does not qualify coverage, residue correspondence, prediction uncertainty, biological orthology, phylogeny or asymmetry.')
    with (out / 'receipt.json').open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'artifacts']}, indent=2), flush=True)


if __name__ == '__main__': main()
