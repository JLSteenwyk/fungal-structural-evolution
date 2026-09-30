#!/usr/bin/env python3
"""Test all context coverage policies, exclusions and corruption rejection."""
import csv
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True); path.write_text(json.dumps(value, indent=2) + '\n')


def table(path, rows):
    with path.open('w') as handle:
        w = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)


def key(ends):
    return hashlib.sha256(json.dumps(sorted(ends), separators=(',', ':')).encode()).hexdigest()


def main():
    with tempfile.TemporaryDirectory(prefix='full-reference-context-coverage-fixture-') as temp:
        root = Path(temp); source = root / 'contexts'; cov = root / 'coverage'; source.mkdir(); cov.mkdir()
        screens = json.loads(Path('metadata/whole_protein_common_fits_plan_20260927.json').read_text())['screens']
        pair_ends = [('R', 10, 'A', 6), ('B', 10, 'R', 10), ('B', 10, 'A', 6)]
        rows = []; pair_info = {}
        for ix, (a, av, b, bv) in enumerate(pair_ends):
            pair = key([(a, av), (b, bv)]); pair_info[pair] = [(a, av), (b, bv)]
            for mask in ['full', 'plddt70']:
                r = dict(pair_key=pair, mask=mask, model_a=a, version_a=str(av), model_b=b, version_b=str(bv))
                for order in [0, 1]:
                    status = 'aligned'; native = 'aligned'; why = ''
                    if ix == 0 and mask == 'plddt70' and order == 1: status = native = 'input_unavailable'; why = 'input_unavailable'
                    if ix == 2:
                        if mask == 'full' and order == 0: status = 'excluded_numerically'; why = 'rmsd_discrepancy'
                        if mask == 'plddt70': status = native = 'native_error' if order == 0 else 'parse_error'; why = native
                    r.update({f'order{order}_status': status, f'order{order}_native_status': native, f'order{order}_numerical_exclusion_reasons': why})
                why = [f'order{o}_not_numerically_usable' for o in [0, 1] if r[f'order{o}_status'] != 'aligned']
                for spec in screens: r[spec['id'] + '_pass'] = '0' if why else '1'; r[spec['id'] + '_exclusions'] = ';'.join(why)
                rows.append(r)
        table(cov / 'pair_mask_order_coverage.tsv', rows)
        native_plan = root / 'native-plan.json'; save(native_plan, dict(inventory='fixture-inventory', base_queue='fixture-queue'))
        union_plan = root / 'union-plan.json'; save(union_plan, dict(native_plan=str(native_plan)))
        cp = root / 'coverage-plan.json'; save(cp, dict(output=str(cov), screens=screens, union_plan=str(union_plan)))
        cov_summary = dict(full_pairs=3, directed_dispositions=12, pair_mask_rows=6, pair_screen_decisions=36,
             order_summary_counts={'synthetic': 6}, maximum_order_differences={}, pair_pass_counts={}, pair_exclusion_counts={}, pair_both_masks_pass_counts={}, source_dispositions={}, native_status_counts={}, numerical_counts={})
        save(cov / 'receipt.json', dict(status='complete_full_reference_order_and_original_coverage_pending_independent_readback', plan_sha256=sha(cp), screens=screens, **cov_summary, artifacts={'pair_mask_order_coverage.tsv': sha(cov / 'pair_mask_order_coverage.tsv')}))
        ca = root / 'coverage-readback.json'; save(ca, dict(status='passed_full_reference_order_and_original_coverage_readback', plan_sha256=sha(cp), producer_receipt_sha256=sha(cov / 'receipt.json'), **cov_summary))
        ar = root / 'coverage-archive.json'; save(ar, dict(services=[dict(synthetic_stub=True)] * 2, summary=cov_summary, source_hashes={str(cov / 'pair_mask_order_coverage.tsv'): sha(cov / 'pair_mask_order_coverage.tsv')}))
        cc = root / 'coverage-completed.json'; save(cc, dict(status='complete_verified_full_reference_order_and_original_coverage', exact_process_journals_checked=2, **cov_summary,
             producer_receipt=str(cov / 'receipt.json'), producer_receipt_sha256=sha(cov / 'receipt.json'), independent_readback=str(ca), independent_readback_sha256=sha(ca), full_hash_archive=str(ar), full_hash_archive_sha256=sha(ar)))
        contexts = []; counts = Counter(); ties = links = availability = 0
        for guide in ['profile', 'mafft']:
            for ix in range(6):
                parent = ix != 3
                native = dict(source_guide=guide, source_row_number=ix + 1, parent_context_eligible=parent, source=dict(family='OG0', gene_node='n' + str(ix), taxon_id='F1', source_tag='fixture'), designs={})
                designs = {}
                for design in ['availability', 'sequence_first']:
                    genes = [('r', 'R', '10')]
                    if ix == 1 and design == 'sequence_first': genes = [('0missing', '', ''), ('r', 'R', '10')]
                    if ix == 4: genes = []
                    if ix == 5: genes = [('a', 'A', '6'), ('c', 'C', '6')]
                    refs = []
                    for index, (gene, model, version) in enumerate(genes):
                        ref = dict(reference_gene='F2_' + gene, reference_model=model, reference_version=version, lexical_choice=index == 0,
                                   native_coorthology=dict(profile='both', mafft='only_a' if ix == 2 else 'both'), native_pair_keys=dict(a='fixture-a', b='fixture-b'))
                        native['designs'][design] = [x['reference'] for x in refs] + [ref]
                        work = {}
                        for side, focal in [('a', ('A', 6)), ('b', ('B', 10))]:
                            reference = (model, int(version)) if model else None; pair = key([focal, reference]) if reference else ''
                            if reference is None: status = 'reference_model_missing'; endpoint = ''
                            elif reference == focal: status = 'identical_model_not_independent'; endpoint = ''
                            elif pair not in pair_info: status = 'outside_full_reference_measurement_design'; endpoint = ''
                            else: status = 'in_full_reference_measurement_design'; endpoint = pair_info[pair].index(focal)
                            work[side] = dict(duplicate_gene=f'F1_{side}{ix}', duplicate_model=focal[0], duplicate_version=focal[1], reference_model=model,
                                              reference_version=int(version) if version else '', pair_key=pair, measurement_disposition=status,
                                              current_pair_focal_endpoint=endpoint, pair_work_disposition='fixture' if pair in pair_info else '', availability_ledger_link=design == 'availability' and parent)
                            counts[f'{guide}|{design}|parent={int(parent)}|{status}'] += 1; links += 1
                            availability += design == 'availability' and parent
                        refs.append(dict(reference=ref, side_work=work)); ties += 1
                    if not genes: counts[f'{guide}|{design}|parent={int(parent)}|no_reference_gene'] += 1
                    designs[design] = refs; native['designs'][design] = [r['reference'] for r in refs]
                contexts.append(dict(native_context=native, measurement_designs=designs, duplicate_models=dict(a=dict(model_id='A', version=6), b=dict(model_id='B', version=10)), duplicate_pair_key=key([('A', 6), ('B', 10)]), duplicate_comparison_status='queued_distinct_models'))
        source_file = source / 'context_measurement_design.jsonl'; source_file.write_text(''.join(json.dumps(r) + '\n' for r in contexts))
        sp = root / 'context-plan.json'; save(sp, dict(output=str(source), inventory='fixture-inventory', queue='fixture-queue'))
        summary = dict(target_contexts=12, context_design_records=24, reference_tie_records=ties, duplicate_reference_links=links,
                       availability_side_links_checked=availability, guide_contexts=dict(profile=6, mafft=6), measurement_disposition_counts=dict(counts))
        save(source / 'receipt.json', dict(status='complete_full_reference_context_measurement_design_pending_independent_readback', plan_sha256=sha(sp), **summary, artifacts={source_file.name: sha(source_file)}))
        sa = root / 'context-proof.json'; save(sa, dict(status='passed_full_reference_context_measurement_design_sql_readback', plan_sha256=sha(sp), producer_receipt_sha256=sha(source / 'receipt.json'), **summary))
        sc = root / 'context-complete.json'; save(sc, dict(status='complete_verified_full_reference_context_measurement_design', services=[dict(synthetic_stub=True)] * 2, summary=summary,
             source_hashes={str(source / 'receipt.json'): sha(source / 'receipt.json'), str(sa): sha(sa), str(source_file): sha(source_file)}))
        expected = dict(target_contexts=12, reference_tie_records=ties, duplicate_reference_links=links, availability_side_links=availability, model_pairs=3)
        plan = root / 'plan.json'; out = root / 'output'; save(plan, dict(context_plan=str(sp), context_readback=str(sa), context_completion=str(sc), coverage_plan=str(cp), coverage_completion=str(cc), screens=screens, expected=expected, output=str(out), pins={}))
        reader = [sys.executable, 'scripts/readback_full_reference_context_coverage.py', '--plan', str(plan), '--output']
        for command in [[sys.executable, 'scripts/project_full_reference_context_coverage.py', '--plan', str(plan)], reader + [str(root / 'proof.json')]]:
            r = subprocess.run(command, capture_output=True, text=True)
            if r.returncode: raise RuntimeError(r.stdout + r.stderr)
        with gzip.open(out / 'context_coverage.jsonl.gz', 'rt') as handle: original = handle.read()
        records = [json.loads(line) for line in original.splitlines()]; receipt = json.loads((out / 'receipt.json').read_text())
        assert records[0]['coverage_designs']['availability']['context_policy_flags']['full']['n50_c70'] == [True] * 5
        assert records[0]['coverage_designs']['availability']['context_policy_flags']['plddt70']['n50_c70'] == [False] * 5
        assert records[1]['coverage_designs']['sequence_first']['context_policy_flags']['full']['n50_c70'] == [False, False, False, True, False]
        assert records[2]['coverage_designs']['availability']['context_policy_flags']['full']['n50_c70'] == [True, True, False, False, False]
        assert records[3]['coverage_designs']['availability']['references'][0]['side_masks']['a']['full']['screen_pass']['n50_c70']
        assert records[3]['coverage_designs']['availability']['context_policy_flags']['full']['n50_c70'] == [False] * 5
        assert records[4]['coverage_designs']['sequence_first']['context_policy_flags']['full']['n50_c70'] == [False] * 5
        rejected = []
        for label in ['promoted_excluded_parent', 'favorable_order', 'modeled_tie_reselection', 'removed_missing_reference', 'empty_all_ties_pass', 'cleared_numerical_flag', 'swapped_side_pair', 'changed_denominator']:
            records = [json.loads(line) for line in original.splitlines()]
            if label == 'promoted_excluded_parent': records[3]['coverage_designs']['availability']['context_policy_flags']['full']['n50_c70'] = [True] * 5
            elif label == 'favorable_order': r = records[0]['coverage_designs']['availability']['references'][0]['side_masks']['a']['plddt70']; r['screen_pass']['n50_c70'] = True; r['screen_exclusions']['n50_c70'] = []
            elif label == 'modeled_tie_reselection': records[1]['coverage_designs']['sequence_first']['context_policy_flags']['full']['n50_c70'] = [True] * 5
            elif label == 'removed_missing_reference': records[1]['coverage_designs']['sequence_first']['references'].pop(0)
            elif label == 'empty_all_ties_pass': records[4]['coverage_designs']['availability']['context_policy_flags']['full']['n50_c70'][-1] = True
            elif label == 'cleared_numerical_flag': records[5]['coverage_designs']['availability']['references'][0]['side_masks']['b']['full']['order_numerical_exclusions'][0] = []
            elif label == 'swapped_side_pair':
                sides = records[0]['coverage_designs']['availability']['references'][0]['side_masks']; sides['a'], sides['b'] = sides['b'], sides['a']
            else:
                summary_rows = list(csv.DictReader((out / 'context_screen_counts.tsv').open(), delimiter='\t')); summary_rows[0]['source_contexts'] = '1'; table(out / 'context_screen_counts.tsv', summary_rows)
            with gzip.open(out / 'context_coverage.jsonl.gz', 'wt') as handle: handle.write(''.join(json.dumps(r) + '\n' for r in records))
            changed = dict(receipt, artifacts={n: sha(out / n) for n in receipt['artifacts']}); save(out / 'receipt.json', changed)
            r = subprocess.run(reader + [str(root / (label + '.json'))], capture_output=True, text=True)
            assert r.returncode != 0 and not (root / (label + '.json')).exists(), label; rejected.append(label)
        print(json.dumps(dict(status='passed_full_synthetic_reference_context_coverage_checks', **expected, context_screen_rows=288, context_policy_decisions=1440,
             rejected_rehashed_exports=rejected, scope='Synthetic closed context/coverage/proof/journal records only, not production qualification or a pilot. All six screens/two masks/both guides/designs, hard excluded-parent gate despite measured pair, missing lexical model with modeled alternate, any/all tie policies, empty all-tie exclusion, native guide disagreement, measured versus null unmeasured orders, raw numerical/error flags and model-version endpoint roles checked. Producer policies independently reconstructed using SQL aggregates.'), indent=2))


if __name__ == '__main__': main()
