#!/usr/bin/env python3
"""Exercise full triad work design invariants and rehashed false-export rejection."""
import csv
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def table(path, rows):
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n'); writer.writeheader(); writer.writerows(rows)


def digest(ends):
    return hashlib.sha256(json.dumps(sorted(ends), separators=(',', ':')).encode()).hexdigest()


def main():
    with tempfile.TemporaryDirectory(prefix='full-reference-triad-fixture-') as temp:
        root = Path(temp); queue = root / 'queue'; inventory = root / 'inventory'; context_root = root / 'context'
        for path in [queue, inventory, context_root]: path.mkdir()
        a, b, r, av, c = ('A', 6), ('B', 10), ('R', 10), ('A', 10), ('C', 6)
        def pair(left, right, work=None):
            row = dict(pair_key=digest([left, right]), model_a=left[0], version_a=left[1], model_b=right[0], version_b=right[1])
            if work: row['work_disposition'] = work
            return row
        primary = [pair(b, a), pair(av, a)]; refs = [pair(a, r, 'additional_pair'), pair(r, b, 'additional_pair'), pair(av, r, 'additional_pair'), pair(b, a, 'existing_duplicate_pair')]
        table(queue / 'model_pairs.tsv', primary); table(inventory / 'model_pairs.tsv', refs)
        ref_catalog = {x['pair_key']: x for x in refs}; model_by_gene = dict(ref=r, ref_alias=r, ref_a=a, ref_c=c, missing=None)
        specs = [(a, b, ['ref', 'ref_alias'], ['ref', 'ref_alias'], True),
                 (a, a, ['ref'], ['ref'], True),
                 (a, b, [], ['missing', 'ref'], True),
                 (a, b, ['ref'], ['ref'], False),
                 (a, b, [], [], True),
                 (a, b, ['ref_a'], ['ref_a'], True),
                 (a, b, [], ['ref_c'], False),
                 (a, av, ['ref'], ['ref'], True)]
        rows = []; availability = ties = 0
        for guide in ['profile', 'mafft']:
            for ix, (am, bm, available, sequence, parent) in enumerate(specs):
                source = dict(family='OG0', gene_node=f'n{ix}', gene_a=f'g{ix}a', gene_b=f'g{ix}b', model_a=am[0], model_b=bm[0],
                              status='provisional_reference_available' if parent else 'parent_reported_duplication', support='0.99', taxon_id='F1')
                models = {s: dict(gene=source['gene_' + s], model_id=m[0], version=m[1]) for s, m in [('a', am), ('b', bm)]}
                designs, native = {}, {}
                for design, genes in [('availability', available), ('sequence_first', sequence)]:
                    records = []; native_records = []
                    for index, gene in enumerate(genes):
                        rm = model_by_gene[gene]; membership = {'profile': 'both', 'mafft': 'only_a' if ix == 0 and gene == 'ref_alias' else 'both'}
                        reference = dict(reference_gene=gene, reference_model=rm[0] if rm else '', reference_version=str(rm[1]) if rm else '',
                                         lexical_choice=index == 0, native_coorthology=membership, native_pair_keys=dict(a='aa', b='bb'))
                        sides = {}
                        for side, focal in [('a', am), ('b', bm)]:
                            key = digest([focal, rm]) if rm else ''; p = ref_catalog.get(key)
                            if rm is None: state = 'reference_model_missing'; endpoint = ''; work = ''
                            elif rm == focal: state = 'identical_model_not_independent'; endpoint = ''; work = ''
                            elif p:
                                state = 'in_full_reference_measurement_design'; endpoint = int(focal != (p['model_a'], p['version_a'])); work = p['work_disposition']
                            else: state = 'outside_full_reference_measurement_design'; endpoint = ''; work = ''
                            linked = design == 'availability'; availability += int(linked)
                            sides[side] = dict(duplicate_gene=models[side]['gene'], duplicate_model=focal[0], duplicate_version=focal[1],
                                               reference_model=rm[0] if rm else '', reference_version=rm[1] if rm else '', pair_key=key, measurement_disposition=state,
                                               current_pair_focal_endpoint=endpoint, pair_work_disposition=work, availability_ledger_link=linked)
                        records.append(dict(reference=reference, side_work=sides)); native_records.append(reference); ties += 1
                    designs[design] = records; native[design] = native_records
                context = dict(source=source, source_guide=guide, source_row_number=ix + 1, parent_context_eligible=parent,
                               source_row_sha256=hashlib.sha256(json.dumps(source, sort_keys=True, separators=(',', ':')).encode()).hexdigest(), designs=native)
                rows.append(dict(native_context=context, duplicate_models=models, duplicate_pair_key=digest([am, bm]),
                                 duplicate_comparison_status='queued_distinct_models' if am != bm else 'identical_model_no_alignment', measurement_designs=designs))
        (queue / 'event_model_pair_links.tsv').write_text('synthetic fixture: full upstream joins are separately tested\n')
        (queue / 'models.jsonl').write_text('synthetic fixture: no native structural coordinates\n')
        save(queue / 'receipt.json', dict(status='complete_reviewed_duplication_model_pair_queue', artifacts={p.name: sha(p) for p in queue.iterdir()}))
        qp = root / 'queue-proof.json'; save(qp, dict(status='passed_full_duplication_model_pair_queue_export_readback', source_sha256={str(p): sha(p) for p in queue.iterdir()}))
        save(inventory / 'receipt.json', dict(unique_distinct_model_pairs=len(refs), artifacts={'model_pairs.tsv': sha(inventory / 'model_pairs.tsv')}))
        source_path = context_root / 'context_measurement_design.jsonl'; source_path.write_text(''.join(json.dumps(row) + '\n' for row in rows))
        cp = root / 'context-plan.json'; save(cp, dict(queue=str(queue), inventory=str(inventory), queue_readback=str(qp), output=str(context_root)))
        summary = dict(target_contexts=len(rows), context_design_records=2 * len(rows), reference_tie_records=ties, duplicate_reference_links=2 * ties,
                       availability_side_links_checked=availability, guide_contexts=dict(profile=8, mafft=8))
        rp = context_root / 'receipt.json'; save(rp, dict(status='complete_full_reference_context_measurement_design_pending_independent_readback', plan_sha256=sha(cp), **summary, artifacts={source_path.name: sha(source_path)}))
        proof = root / 'context-proof.json'; save(proof, dict(status='passed_full_reference_context_measurement_design_sql_readback', plan_sha256=sha(cp), producer_receipt_sha256=sha(rp), **summary))
        closed = root / 'context-complete.json'; save(closed, dict(status='complete_verified_full_reference_context_measurement_design', summary=summary, services=[dict(synthetic_stub=True)] * 2,
                        source_hashes={str(p): sha(p) for p in [rp, proof, source_path, inventory / 'model_pairs.tsv', *queue.iterdir()]}))
        plan = root / 'plan.json'; out = root / 'out'; expected = dict(target_contexts=len(rows), reference_tie_records=ties, duplicate_reference_links=2 * ties,
                          availability_side_links=availability, primary_pairs=len(primary), reference_pairs=len(refs))
        save(plan, dict(context_plan=str(cp), context_readback=str(proof), context_completion=str(closed), expected=expected, pins={}, output=str(out), scope='Synthetic fixture only.'))
        reader = [sys.executable, 'scripts/readback_full_reference_triad_design.py', '--plan', str(plan), '--output']
        for command in [[sys.executable, 'scripts/prepare_full_reference_triad_design.py', '--plan', str(plan)], reader + [str(root / 'passed.json')]]:
            result = subprocess.run(command, capture_output=True, text=True)
            if result.returncode: raise RuntimeError(result.stdout + result.stderr)
        path = out / 'context_triad_design.jsonl.gz'; physical_path = out / 'ordered_model_triads.jsonl'
        with gzip.open(path, 'rt') as handle: original = [json.loads(line) for line in handle]
        physical_original = [json.loads(line) for line in physical_path.read_text().splitlines()]; receipt = json.loads((out / 'receipt.json').read_text())
        assert original[0]['duplicate_pair_design']['desired_direction_order'] == 1
        assert original[0]['triad_designs']['availability'][0]['source_design_ready_for_correspondence']
        assert not original[1]['triad_designs']['availability'][0]['source_design_ready_for_correspondence']
        assert original[2]['triad_designs']['sequence_first'][0]['triad_id'] is None and original[2]['triad_designs']['sequence_first'][1]['source_design_ready_for_correspondence']
        assert not original[3]['triad_designs']['availability'][0]['source_design_ready_for_correspondence']
        assert original[4]['triad_designs'] == dict(availability=[], sequence_first=[])
        assert not original[5]['triad_designs']['availability'][0]['three_distinct_versioned_models']
        assert not original[6]['triad_designs']['sequence_first'][0]['all_three_pairs_in_measurement_design']
        version_case = original[7]['triad_designs']['availability'][0]
        assert version_case['three_distinct_versioned_models'] and not version_case['three_distinct_model_ids'] and not version_case['source_design_ready_for_correspondence']
        rejected = []
        for label in ['promoted_parent', 'promoted_duplicate_identity', 'dropped_missing_reference', 'reselected_modeled_tie', 'changed_ab_direction', 'swapped_ar_br_roles', 'promoted_model_versions', 'changed_logical_link_count', 'removed_empty_context', 'changed_pair_hash']:
            altered = json.loads(json.dumps(original)); triads = json.loads(json.dumps(physical_original))
            if label == 'promoted_parent': altered[3]['triad_designs']['availability'][0]['source_design_ready_for_correspondence'] = True
            elif label == 'promoted_duplicate_identity': altered[1]['triad_designs']['availability'][0]['source_design_ready_for_correspondence'] = True
            elif label == 'dropped_missing_reference': altered[2]['triad_designs']['sequence_first'].pop(0)
            elif label == 'reselected_modeled_tie': altered[2]['source_design']['measurement_designs']['sequence_first'][0]['reference']['lexical_choice'] = False; altered[2]['source_design']['measurement_designs']['sequence_first'][1]['reference']['lexical_choice'] = True
            elif label == 'changed_ab_direction': altered[0]['duplicate_pair_design']['desired_direction_order'] = 0
            elif label == 'swapped_ar_br_roles': triads[0]['edges']['ar'], triads[0]['edges']['br'] = triads[0]['edges']['br'], triads[0]['edges']['ar']
            elif label == 'promoted_model_versions': altered[7]['triad_designs']['availability'][0]['source_design_ready_for_correspondence'] = True
            elif label == 'changed_logical_link_count': triads[0]['logical_reference_links'] += 1
            elif label == 'removed_empty_context': altered.pop(4)
            else: triads[0]['edges']['ab']['pair_key'] = '0' * 64
            with gzip.open(path, 'wt', compresslevel=1) as handle: handle.write(''.join(json.dumps(row) + '\n' for row in altered))
            physical_path.write_text(''.join(json.dumps(row) + '\n' for row in triads))
            save(out / 'receipt.json', dict(receipt, artifacts={p.name: sha(p) for p in [path, physical_path]}))
            result = subprocess.run(reader + [str(root / (label + '.json'))], capture_output=True, text=True)
            assert result.returncode and not (root / (label + '.json')).exists(), label; rejected.append(label)
        print(json.dumps(dict(status='passed_full_synthetic_reference_triad_work_design_checks', **expected,
             unique_ordered_model_triads=receipt['unique_ordered_model_triads'], correspondence_work_triads=receipt['correspondence_work_triads'],
             rejected_rehashed_exports=rejected, scope='Synthetic source/proof/journal stubs only; not a production qualification or a sampling pilot. Both guides/designs, original lexical ties and missing reference models, excluded parents with measured physical pairs, empty contexts, identical duplicate/reference models, distinct versions sharing one model ID, unequal model versions, reversed primary/reference endpoints, native guide disagreement and reference aliases exercised. SQL reader shares no producer projection functions.'), indent=2))


if __name__ == '__main__': main()
