#!/usr/bin/env python3
"""Synthetic full-context linkage/null/parent/tie-policy software checks."""
import argparse
import copy
import csv
import gzip
import hashlib
import json
import tempfile
from pathlib import Path
import project_full_triad_context_geometry as producer
import readback_full_triad_context_geometry as reader
from full_triad_context_geometry_sources import GATES, DEFINITIONS
from run_ortholog_pair_guide_comparison import sha


def digest(models): return hashlib.sha256(json.dumps(models, separators=(',', ':')).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args(); assert not args.output.exists()
    screens = [dict(id='n50_c70', minimum_aligned_residues=50, minimum_original_coverage=.7)]
    with tempfile.TemporaryDirectory(prefix='full-triad-context-software-', dir='results') as folder:
        temp = Path(folder); original = temp / 'contexts.jsonl.gz'; measurements = temp / 'robustness.jsonl.gz'
        models = [[['a', 1], ['b', 1], ['r' + str(i), 1]] for i in range(3)]
        with gzip.open(measurements, 'wt') as handle:
            for i, roles in enumerate(models):
                for mask in ['full', 'plddt70']:
                    for definition in DEFINITIONS:
                        bits = 255 if i == 0 else 255 if i == 1 and mask == 'full' else 0 if i == 1 else 170 if mask == 'full' else 85
                        handle.write(json.dumps(dict(triad_id=digest(roles), models=roles, mask=mask, mapping_definition=definition, screens={'n50_c70': dict(order_pass_bits=bits)})) + '\n')
        records = []
        for guide in ['profile', 'mafft']:
            for case in range(10):
                parent = case != 1; all_pairs = case != 7
                choices = [('r0', 1)]
                if case == 2: choices = [('', ''), ('r0', 1)]
                if case == 3: choices = []
                if case == 4: choices = [('a', 2)]
                if case == 6: choices = [('r0', 1), ('r1', 1)]
                if case == 7: choices = [('r3', 1)]
                if case == 8: choices = [('r2', 1)]
                if case == 9: choices = [('r1', 1)]
                co = dict(profile='both', mafft='only_a') if case == 5 else dict(profile='both', mafft='both')
                designs, linked = {}, {}
                for design in ['availability', 'sequence_first']:
                    refs, links = [], []
                    for ix, (model, version) in enumerate(choices):
                        roles = [['a', 1], ['b', 1], [model, version]] if model else None
                        ref = dict(reference_gene='ref-' + str(ix), reference_model=model, reference_version=str(version), lexical_choice=ix == 0, native_coorthology=co)
                        refs.append(dict(reference=ref, side_work={side: dict(measurement_disposition='in_full_reference_measurement_design' if model and all_pairs else 'outside_full_reference_measurement_design') for side in ['a', 'b']}))
                        distinct = bool(roles and len({tuple(m) for m in roles}) == 3); ids = bool(roles and len({m[0] for m in roles}) == 3); in_design = bool(model and all_pairs)
                        ready = bool(parent and distinct and ids and in_design)
                        flags = [ready, ready and co[guide] == 'both', ready and co['profile'] == co['mafft'] == 'both']
                        links.append(dict(reference_gene=ref['reference_gene'], triad_id=digest(roles) if model else None, three_distinct_versioned_models=distinct,
                                          three_distinct_model_ids=ids, all_three_pairs_in_measurement_design=in_design, **dict(zip(GATES, flags)), work_exclusion_reasons=[] if ready else ['synthetic_source_exclusion']))
                    designs[design] = refs; linked[design] = links
                records.append(dict(source_design=dict(native_context=dict(source_guide=guide, source_row_number=case + 1, parent_context_eligible=parent,
                                                                           source=dict(family='synthetic-family', taxon_id='synthetic-taxon', gene_node='node-' + str(case))),
                                                        duplicate_models={side: dict(model_id=side, version=1, gene='gene-' + side) for side in ['a', 'b']},
                                                        duplicate_comparison_status='queued_distinct_models', measurement_designs=designs),
                                    duplicate_pair_design=dict(measurement_disposition='in_full_primary_measurement_design' if all_pairs else 'outside_full_primary_measurement_design'), triad_designs=linked))
        with gzip.open(original, 'wt') as handle:
            for value in records: handle.write(json.dumps(value) + '\n')
        ties = sum(len(r['triad_designs'][d]) for r in records for d in ['availability', 'sequence_first'])
        def sources(plan, path): return original, measurements, dict(guide_contexts=dict(profile=10, mafft=10)), {str(path): sha(path), str(original): sha(original), str(measurements): sha(measurements)}
        producer.load_sources = reader.load_sources = sources
        plan = dict(output=str(temp / 'baseline'), screens=screens, resources=dict(minimum_free_disk_gib=0), pins={},
                    expected=dict(target_contexts=20, reference_tie_records=ties, duplicate_reference_links=2 * ties, measured_triads=3, measured_robustness_groups=12),
                    scope='Synthetic full-context software fixture only; source/proof I/O stubbed; not a biological pilot or source qualification.')
        pp = temp / 'baseline-plan.json'; pp.write_text(json.dumps(plan)); producer.run(pp); baseline = reader.run(pp, temp / 'baseline-readback.json')
        with gzip.open(Path(plan['output']) / 'context_geometry.jsonl.gz', 'rt') as handle: groups = [json.loads(line) for line in handle]
        def policies(i, mask='both_masks'): return groups[i]['structural_designs']['availability']['context_policy_flags'][mask]['reference_common']['n50_c70']
        assert policies(0) == [True] * 5 and policies(1) == [False] * 5
        assert groups[1]['structural_designs']['availability']['references'][0]['physical_measurement_present'] is True
        assert policies(2) == [False, False, False, True, False] and policies(3) == [False] * 5
        assert policies(4) == [False] * 5 and policies(5) == [True, True, False, False, False] and policies(15) == [True, False, False, False, False]
        assert policies(6) == [True, True, True, True, False] and policies(6, 'full') == [True] * 5
        assert policies(7) == policies(8) == policies(9) == [False] * 5
        assert policies(9, 'full') == [True] * 5
        def reference(g, i): return g[i]['structural_designs']['availability']['references'][0]
        mutations = {
            'promoted_excluded_parent': lambda g, r: reference(g, 1)['source_eligibility_flags'].__setitem__(0, True),
            'unmeasured_bitmap_as_zero': lambda g, r: reference(g, 7)['order_pass_bits']['full']['reference_common'].__setitem__('n50_c70', 0),
            'favorable_reference_reselection': lambda g, r: g[2]['structural_designs']['availability']['references'].reverse(),
            'dropped_unmodeled_tie': lambda g, r: g[2]['structural_designs']['availability']['references'].pop(0),
            'vacuous_empty_all_ties': lambda g, r: g[3]['structural_designs']['availability']['context_policy_flags']['both_masks']['reference_common']['n50_c70'].__setitem__(4, True),
            'promoted_native_both_guides': lambda g, r: reference(g, 5)['source_eligibility_flags'].__setitem__(2, True),
            'invented_both_mask_bitmap': lambda g, r: reference(g, 8)['order_pass_bits']['both_masks']['reference_common'].__setitem__('n50_c70', 255),
            'changed_physical_result_key': lambda g, r: reference(g, 0)['physical_robustness_keys'][0].__setitem__(0, 'wrong'),
            'changed_original_gene_context': lambda g, r: g[0]['source_work_design']['source_design']['native_context']['source'].update(taxon_id='changed'),
            'duplicated_context': lambda g, r: g.insert(1, copy.deepcopy(g[0])),
            'missing_context': lambda g, r: g.pop(),
            'invented_scientific_eligibility': lambda g, r: r.update(scientific_eligibility=True),
        }
        base_receipt = json.loads((Path(plan['output']) / 'receipt.json').read_text()); rejected = []
        for name, mutate in mutations.items():
            candidate = copy.deepcopy(groups); receipt = copy.deepcopy(base_receipt); mutate(candidate, receipt)
            case = temp / name; case.mkdir(); cp = temp / (name + '-plan.json'); config = {**plan, 'output': str(case)}; cp.write_text(json.dumps(config))
            with gzip.open(case / 'context_geometry.jsonl.gz', 'wt') as handle:
                for record in candidate: handle.write(json.dumps(record) + '\n')
            (case / 'context_geometry_counts.tsv').write_bytes((Path(plan['output']) / 'context_geometry_counts.tsv').read_bytes())
            receipt.update(plan_sha256=sha(cp), source_hashes=sources(config, cp)[-1], artifacts={n: sha(case / n) for n in ['context_geometry.jsonl.gz', 'context_geometry_counts.tsv']})
            (case / 'receipt.json').write_text(json.dumps(receipt))
            try: reader.run(cp, case / 'readback.json')
            except AssertionError: rejected.append(name)
            else: raise AssertionError('False export accepted: ' + name)
        assert len(rejected) == len(mutations)
    result = dict(status='passed_synthetic_full_triad_context_geometry_software_checks', synthetic_contexts=20, synthetic_reference_ties=ties,
                  baseline_readback_status=baseline['status'], rejected_rehashed_false_exports=rejected,
                  script_hashes={str(path): sha(path) for path in [Path(__file__), Path(producer.__file__), Path(reader.__file__), Path('scripts/full_triad_context_geometry_sources.py')]},
                  scope='All10synthetic context cases in both guides/bothdesigns: measured excluded parents, missing lexical model followed by measured ties, empty tie sets, same modelID/differentversions, own vsboth native guides, full/p70 discordance, complementary orderbits, unscheduled links. Independent SQL baseline and12 rehashed false export rejections. Source/proof I/O stubbed only here; not a biological pilot or real source/geometry/scientific qualification.')
    with args.output.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__': main()
