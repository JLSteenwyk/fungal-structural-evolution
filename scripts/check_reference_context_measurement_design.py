#!/usr/bin/env python3
"""Test full context work states, reversed gene roles and corrupt projection rejection."""
import csv
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
    with tempfile.TemporaryDirectory(prefix='reference-context-measurement-fixture-') as temp:
        root = Path(temp); native = root / 'native'; queue = root / 'queue'; inventory = root / 'inventory'
        for p in [native, queue, inventory]: p.mkdir()
        models = dict(A=6, B=10, R=10, C=6)
        pairs = []
        for a, b, disposition in [('R', 'A', 'additional_pair'), ('B', 'R', 'additional_pair'), ('B', 'A', 'existing_duplicate_pair')]:
            pairs.append(dict(pair_key=key([(a, models[a]), (b, models[b])]), model_a=a, version_a=models[a], model_b=b, version_b=models[b], work_disposition=disposition))
        table(inventory / 'model_pairs.tsv', pairs)
        table(queue / 'model_pairs.tsv', [pairs[-1]])
        (queue / 'models.jsonl').write_text(''.join(json.dumps(dict(model_id=m, version=v)) + '\n' for m, v in models.items()))
        qrows, contexts, ledger = [], [], []
        specs = [(['F2_r'], ['F2_r'], 'provisional_reference_available'),
                 (['F2_a'], ['F2_a'], 'provisional_reference_available'),
                 ([], ['F2_0missing', 'F2_r'], 'no_modeled_nonfocal_sister'),
                 ([], ['F2_c'], 'parent_reported_duplication'),
                 ([], [], 'no_modeled_nonfocal_sister'),
                 (['F2_r', 'F2_ralias'], ['F2_r', 'F2_ralias'], 'provisional_reference_available')]
        refmodels = {'F2_r': ('R', '10'), 'F2_ralias': ('R', '10'), 'F2_a': ('A', '6'), 'F2_c': ('C', '6'), 'F2_0missing': ('', '')}
        for guide in ['profile', 'mafft']:
            for ix, (available, sequence, status) in enumerate(specs):
                ma, mb = 'A', 'A' if ix == 5 else 'B'
                source = dict(family='OG0000000', taxon_id='F1', gene_node='n' + str(ix), gene_a=f'F1_a{ix}', gene_b=f'F1_b{ix}',
                              model_a=ma, model_b=mb, status=status, parent_reported_duplication=str(int(ix == 3)), support='0.9')
                q = dict(source, guide=guide, version_a=str(models[ma]), version_b=str(models[mb]), pair_key=key([(ma, models[ma]), (mb, models[mb])]),
                         same_model=str(int(ma == mb)), comparison_status='identical_model_no_alignment' if ma == mb else 'queued_distinct_models')
                # Primary queue positions differ from the reference-source positions.
                if ix == 0:
                    for field in ['gene', 'model', 'version']: q[field + '_a'], q[field + '_b'] = q[field + '_b'], q[field + '_a']
                qrows.append(q)
                designs = {}
                for design, genes in [('availability', available), ('sequence_first', sequence)]:
                    refs = []
                    for i, gene in enumerate(genes):
                        rm, rv = refmodels[gene]
                        refs.append(dict(reference_gene=gene, reference_model=rm, reference_version=rv, lexical_choice=i == 0,
                                         native_pair_keys=dict(a='000001000002', b='000002000003'),
                                         native_coorthology=dict(profile='only_a' if ix == 2 else 'both', mafft='both')))
                        if design == 'availability':
                            for side, focal in [('a', ma), ('b', mb)]:
                                k = key([(focal, models[focal]), (rm, int(rv))])
                                work = 'identical_model' if focal == rm else next(p['work_disposition'] for p in pairs if p['pair_key'] == k)
                                ledger.append(dict(guide=guide, family=source['family'], gene_node=source['gene_node'], gene_a=source['gene_a'], gene_b=source['gene_b'],
                                              reference_gene=gene, lexical_representative=str(int(i == 0)), focal_side=side, focal_gene=source['gene_' + side],
                                              focal_model=focal, focal_version=str(models[focal]), reference_model=rm, reference_version=rv, pair_key=k, work_disposition=work))
                    designs[design] = refs
                contexts.append(dict(source=source, source_guide=guide, source_row_number=ix + 1,
                                     source_row_sha256=hashlib.sha256(json.dumps(source, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
                                     parent_context_eligible=status != 'parent_reported_duplication', designs=designs))
        table(queue / 'event_model_pair_links.tsv', qrows); table(inventory / 'event_reference_comparisons.tsv', ledger)
        save(queue / 'receipt.json', dict(status='complete_reviewed_duplication_model_pair_queue', counts=dict(event_links=12), artifacts={p.name: sha(p) for p in queue.iterdir()}))
        qp = root / 'queue-proof.json'; save(qp, dict(status='passed_full_duplication_model_pair_queue_export_readback', source_sha256={str(p): sha(p) for p in queue.iterdir()}))
        save(inventory / 'receipt.json', dict(status='complete_provisional_reference_comparison_inventory', unique_distinct_model_pairs=3,
             source_pins={str(queue / 'receipt.json'): sha(queue / 'receipt.json')}, artifacts={p.name: sha(p) for p in inventory.iterdir()}))
        ip = root / 'inventory-proof.json'; save(ip, dict(status='passed_full_reference_comparison_ledger_and_native_model_readback', producer_receipt_sha256=sha(inventory / 'receipt.json'),
             event_reference_comparisons=len(ledger), distinct_model_pairs=3, primary_queue_readback_sha256=sha(qp)))
        ic = root / 'inventory-complete.json'; save(ic, dict(status='complete_verified_expanded_duplication_reference_pipeline', services=[dict(synthetic_stub=True)] * 6,
             summary=dict(event_reference_comparisons=len(ledger)), source_hashes={str(inventory / 'receipt.json'): sha(inventory / 'receipt.json'), str(ip): sha(ip)}))
        context_file = native / 'contexts.jsonl'; context_file.write_text(''.join(json.dumps(r) + '\n' for r in contexts))
        ties = sum(len(refs) for c in contexts for refs in c['designs'].values())
        expected = dict(target_contexts=12, reference_tie_records=ties, duplicate_reference_links=2 * ties, availability_side_links=len(ledger), model_pairs=3)
        summary = {k: expected[k] for k in ['target_contexts', 'reference_tie_records', 'duplicate_reference_links']}
        save(native / 'receipt.json', dict(status='complete_full_reference_native_orthology_pending_readback', **summary, guide_contexts=dict(profile=6, mafft=6), artifacts={context_file.name: sha(context_file)}))
        np = root / 'native-proof.json'; save(np, dict(status='passed_full_reference_native_orthology_readback', producer_receipt_sha256=sha(native / 'receipt.json'), **summary))
        nc = root / 'native-complete.json'; save(nc, dict(status='complete_verified_full_reference_native_orthology', services=[dict(synthetic_stub=True)] * 2,
             summary=dict(**summary, guide_contexts=dict(profile=6, mafft=6)), source_hashes={str(native / 'receipt.json'): sha(native / 'receipt.json'), str(np): sha(np), str(context_file): sha(context_file)}))
        plan = root / 'plan.json'; out = root / 'output'; save(plan, dict(native_contexts=str(native), native_readback=str(np), native_completion=str(nc),
             queue=str(queue), queue_readback=str(qp), inventory=str(inventory), inventory_readback=str(ip), inventory_completion=str(ic), expected=expected, output=str(out), pins={}))
        reader = [sys.executable, 'scripts/readback_reference_context_measurement_design.py', '--plan', str(plan), '--output']
        for command in [[sys.executable, 'scripts/prepare_reference_context_measurement_design.py', '--plan', str(plan)], reader + [str(root / 'proof.json')]]:
            r = subprocess.run(command, capture_output=True, text=True)
            if r.returncode: raise RuntimeError(r.stdout + r.stderr)
        output = out / 'context_measurement_design.jsonl'; original = output.read_text(); receipt = json.loads((out / 'receipt.json').read_text())
        rows = [json.loads(line) for line in original.splitlines()]
        assert rows[0]['duplicate_models']['a']['model_id'] == 'A' and rows[0]['duplicate_models']['b']['version'] == 10
        assert rows[0]['measurement_designs']['availability'][0]['side_work']['a']['current_pair_focal_endpoint'] == 1
        assert rows[1]['measurement_designs']['availability'][0]['side_work']['a']['measurement_disposition'] == 'identical_model_not_independent'
        assert rows[2]['measurement_designs']['sequence_first'][0]['reference']['lexical_choice'] and rows[2]['measurement_designs']['sequence_first'][0]['side_work']['a']['measurement_disposition'] == 'reference_model_missing'
        assert rows[3]['measurement_designs']['sequence_first'][0]['side_work']['a']['measurement_disposition'] == 'outside_full_reference_measurement_design'
        assert rows[4]['measurement_designs'] == dict(availability=[], sequence_first=[])
        rejected = []
        for label in ['swapped_duplicate_models', 'changed_version', 'changed_current_endpoint', 'cleared_parent_exclusion', 'reselected_modeled_tie', 'dropped_empty_context', 'deleted_alias_reference', 'changed_pair_hash']:
            rows = [json.loads(line) for line in original.splitlines()]
            if label == 'swapped_duplicate_models': rows[0]['duplicate_models']['a'], rows[0]['duplicate_models']['b'] = rows[0]['duplicate_models']['b'], rows[0]['duplicate_models']['a']
            elif label == 'changed_version': rows[0]['duplicate_models']['b']['version'] = 6
            elif label == 'changed_current_endpoint': rows[0]['measurement_designs']['availability'][0]['side_work']['a']['current_pair_focal_endpoint'] = 0
            elif label == 'cleared_parent_exclusion': rows[3]['native_context']['parent_context_eligible'] = True
            elif label == 'reselected_modeled_tie': rows[2]['measurement_designs']['sequence_first'][0]['reference']['lexical_choice'] = False; rows[2]['measurement_designs']['sequence_first'][1]['reference']['lexical_choice'] = True
            elif label == 'dropped_empty_context': rows.pop(4)
            elif label == 'deleted_alias_reference': rows[5]['measurement_designs']['availability'].pop()
            else: rows[0]['measurement_designs']['availability'][0]['side_work']['a']['pair_key'] = '0' * 64
            output.write_text(''.join(json.dumps(r) + '\n' for r in rows)); changed = dict(receipt, artifacts={output.name: sha(output)}); save(out / 'receipt.json', changed)
            r = subprocess.run(reader + [str(root / (label + '.json'))], capture_output=True, text=True)
            assert r.returncode != 0 and not (root / (label + '.json')).exists(), label; rejected.append(label)
        print(json.dumps(dict(status='passed_full_synthetic_reference_context_measurement_design_checks', **expected, rejected_rehashed_exports=rejected,
             scope='Synthetic context/queue/native-proof/journal records only, not production qualification or a pilot. Both guides/designs, reversed queue gene roles, unequal versions, nonlexical model endpoints, unmodeled lexical tie with modeled alternate, excluded parent, no reference, identical duplicate/reference models, shared-model reference aliases, native guide disagreement and all four measurement work states exercised. Independent SQL reconstruction shares only source/proof I/O.'), indent=2))


if __name__ == '__main__': main()
