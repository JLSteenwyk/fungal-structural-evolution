#!/usr/bin/env python3
"""Link full structural order grids to every original context and tied reference."""
import argparse
import csv
import gzip
import json
import shutil
from collections import Counter
from pathlib import Path
from full_triad_context_geometry_sources import load_sources, DEFINITIONS, MASKS, GATES, POLICIES
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def run(plan_path):
    plan_path = Path(plan_path); plan = json.loads(plan_path.read_text())
    assert shutil.disk_usage('.').free >= plan['resources']['minimum_free_disk_gib'] * 2 ** 30
    contexts, measurements, upstream, bindings = load_sources(plan, plan_path)
    index = {}; groups = 0
    with gzip.open(measurements, 'rt') as handle:
        for line in handle:
            row = json.loads(line); tid, mask, definition = row['triad_id'], row['mask'], row['mapping_definition']
            assert mask in MASKS[:2] and definition in DEFINITIONS
            state = index.setdefault(tid, dict(models=row['models'], masks={}))
            assert state['models'] == row['models']; by_definition = state['masks'].setdefault(mask, {}); assert definition not in by_definition
            by_definition[definition] = {s['id']: row['screens'][s['id']]['order_pass_bits'] for s in plan['screens']}; groups += 1
            assert all(type(v) is int and 0 <= v <= 255 for v in by_definition[definition].values())
    assert len(index) == plan['expected']['measured_triads'] and groups == plan['expected']['measured_robustness_groups']
    assert all(set(s['masks']) == set(MASKS[:2]) and all(set(v) == set(DEFINITIONS) for v in s['masks'].values()) for s in index.values())
    out = Path(plan['output']); out.mkdir(exist_ok=False, parents=True)
    guides, presence, baseline, counts = Counter(), Counter(), Counter(), Counter(); n = ties = 0
    with gzip.open(contexts, 'rt') as source, gzip.open(out / 'context_geometry.jsonl.gz', 'xt', compresslevel=1) as dest:
        for line in source:
            original = json.loads(line); native = original['source_design']['native_context']; guide = native['source_guide']; parent = native['parent_context_eligible']
            projected = {}
            for design in ['availability', 'sequence_first']:
                references = original['triad_designs'][design]; measurement_refs = original['source_design']['measurement_designs'][design]
                assert len(references) == len(measurement_refs); linked = []
                lexical = measurement_refs[0]['reference'] if measurement_refs else None
                for label, value in [('source_contexts', 1), ('parent_eligible_contexts', parent), ('parent_eligible_lexical_gene_present', bool(parent and lexical)), ('parent_eligible_lexical_model_present', bool(parent and lexical and lexical['reference_model']))]:
                    baseline[guide, design, label] += value
                for ix, (ref, source_ref) in enumerate(zip(references, measurement_refs)):
                    chosen = source_ref['reference']; assert chosen['lexical_choice'] == (ix == 0) and chosen['reference_gene'] == ref['reference_gene']
                    flags = [ref[gate] for gate in GATES]; assert all(type(v) is bool for v in flags)
                    assert not any(flags) or parent
                    state = index.get(ref['triad_id']); measured = state is not None
                    if measured:
                        models = [[original['source_design']['duplicate_models'][role]['model_id'], original['source_design']['duplicate_models'][role]['version']] for role in ['a', 'b']]
                        models.append([chosen['reference_model'], int(chosen['reference_version'])]); assert state['models'] == models
                    if flags[0]: assert measured
                    bits = {}; qualified = {}
                    for mask in MASKS:
                        bits[mask] = {}; qualified[mask] = {}
                        for definition in DEFINITIONS:
                            cells = {}
                            for screen in plan['screens']:
                                sid = screen['id']; value = None
                                if measured:
                                    value = state['masks'][mask][definition][sid] if mask != 'both_masks' else state['masks']['full'][definition][sid] & state['masks']['plddt70'][definition][sid]
                                cells[sid] = value
                            bits[mask][definition] = cells
                            qualified[mask][definition] = {sid: [bool(flag and value == 255) for flag in flags] for sid, value in cells.items()}
                    linked.append(dict(reference_gene=ref['reference_gene'], triad_id=ref['triad_id'], physical_measurement_present=measured,
                                       physical_robustness_keys=[[ref['triad_id'], mask, definition] for mask in MASKS[:2] for definition in DEFINITIONS] if measured else None,
                                       source_eligibility_flags=flags, order_pass_bits=bits, qualified_all_order_flags=qualified))
                    presence[f'{guide}|{design}|parent={int(parent)}|measured={int(measured)}|source_ready={int(flags[0])}'] += 1; ties += 1
                policies = {}
                for mask in MASKS:
                    policies[mask] = {}
                    for definition in DEFINITIONS:
                        cells = {}
                        for screen in plan['screens']:
                            sid = screen['id']; q = [ref['qualified_all_order_flags'][mask][definition][sid] for ref in linked]
                            first = q[0] if q else [False] * 3
                            cells[sid] = [*first, any(v[2] for v in q), bool(q and all(v[2] for v in q))]
                            for policy, value in zip(POLICIES, cells[sid]): counts[guide, design, mask, definition, sid, policy] += value
                        policies[mask][definition] = cells
                projected[design] = dict(references=linked, context_policy_flags=policies)
            dest.write(json.dumps(dict(source_work_design=original, structural_designs=projected), separators=(',', ':')) + '\n'); n += 1; guides[guide] += 1
            if n % 25000 == 0: print('Full original context geometric linkage', n, '/', plan['expected']['target_contexts'], flush=True)
    assert n == plan['expected']['target_contexts'] and ties == plan['expected']['reference_tie_records'] and dict(guides) == upstream['guide_contexts']
    fields = ['guide', 'design', 'mask', 'mapping_definition', 'screen', 'policy', 'source_contexts', 'parent_eligible_contexts', 'parent_eligible_lexical_gene_present', 'parent_eligible_lexical_model_present', 'passed_contexts']
    with (out / 'context_geometry_counts.tsv').open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t'); writer.writeheader()
        for key, count in sorted(counts.items()):
            value = dict(zip(fields[:6], key)); value.update({label: baseline[key[0], key[1], label] for label in fields[6:-1]}); value['passed_contexts'] = count; writer.writerow(value)
    summary = dict(target_contexts=n, context_design_records=2 * n, reference_tie_records=ties, duplicate_reference_links=2 * ties,
                   logical_reference_screen_decisions=ties * len(MASKS) * len(DEFINITIONS) * len(plan['screens']),
                   context_screen_states=n * 2 * len(MASKS) * len(DEFINITIONS) * len(plan['screens']),
                   context_policy_decisions=n * 2 * len(MASKS) * len(DEFINITIONS) * len(plan['screens']) * len(POLICIES),
                   summary_rows=len(counts), guide_contexts=dict(guides), reference_measurement_presence_counts=dict(presence), measured_triads=len(index), measured_robustness_groups=groups)
    assert summary['duplicate_reference_links'] == plan['expected']['duplicate_reference_links']; verify(bindings)
    result = dict(status='complete_full_triad_context_geometry_pending_independent_readback', plan_sha256=sha(plan_path), **summary,
                  masks=MASKS, definitions=DEFINITIONS, source_gate_order=GATES, context_policy_order=POLICIES, screens=plan['screens'], source_hashes=bindings,
                  artifacts={name: sha(out / name) for name in ['context_geometry.jsonl.gz', 'context_geometry_counts.tsv']}, scientific_eligibility=False, scope=plan['scope'])
    with (out / 'receipt.json').open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', required=True, type=Path); run(parser.parse_args().plan)
