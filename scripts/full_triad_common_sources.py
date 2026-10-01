"""Closed full three-way design, written input and measured-edge source I/O."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from reference_measurement_union_sources import bind, verify, artifacts
from reference_triad_design_sources import SUMMARY_FIELDS
from run_ortholog_pair_guide_comparison import sha


def load_design_inputs(plan):
    bindings = dict(plan['pins']); config = json.loads(Path(plan['triad_plan']).read_text())
    root = Path(config['output']); rp = root / 'receipt.json'
    r = json.loads(rp.read_text()); a = json.loads(Path(plan['triad_readback']).read_text()); c = json.loads(Path(plan['triad_completion']).read_text())
    assert r['status'] == 'complete_full_reference_triad_work_design_pending_independent_readback'
    assert a['status'] == 'passed_full_reference_triad_work_design_sql_readback'
    assert c['status'] == 'complete_verified_full_reference_triad_work_design' and len(c['services']) == 2
    assert r['plan_sha256'] == a['plan_sha256'] == sha(plan['triad_plan'])
    assert a['producer_receipt_sha256'] == sha(rp) == c['source_hashes'][str(rp)]
    assert c['source_hashes'][plan['triad_readback']] == sha(plan['triad_readback'])
    assert all(r[k] == a[k] == c['summary'][k] for k in SUMMARY_FIELDS)
    for k in ['target_contexts', 'reference_tie_records', 'duplicate_reference_links', 'unique_ordered_model_triads', 'correspondence_work_triads', 'potential_correspondence_states']:
        assert r[k] == plan['expected'][k]
    for p, h in c['source_hashes'].items(): bind(bindings, p, h)
    artifacts(bindings, root, r)
    context_config = json.loads(Path(config['context_plan']).read_text())
    native_config = json.loads(Path(plan['reference_native_plan']).read_text())
    primary_config = json.loads(Path(plan['primary_native_plan']).read_text())
    assert native_config['inventory'] == context_config['inventory'] and native_config['base_queue'] == primary_config['queue'] == context_config['queue']
    triads = []
    with (root / 'ordered_model_triads.jsonl').open() as handle:
        for line in handle:
            row = json.loads(line)
            if row['source_design_ready_links']:
                assert row['distinct_versioned_models'] == row['distinct_model_ids'] == 3 and row['all_three_pairs_in_measurement_design']
                assert not row['model_work_exclusion_reasons']; triads.append(row)
    assert len(triads) == plan['expected']['correspondence_work_triads']
    needed = {tuple(model) for t in triads for model in t['models']}
    catalog = {}
    inventory = Path(context_config['inventory'])
    with (inventory / 'models.jsonl').open() as handle:
        for line in handle:
            row = json.loads(line); key = row['model_id'], row['version']; assert key not in catalog
            catalog[key] = row
    assert len(catalog) == plan['expected']['reference_models'] and needed <= set(catalog)
    inputs, seen = {}, set(); source_counts = {}; row_count = 0
    assert primary_config['inputs'] == native_config['input_sources']['primary']['inputs']
    for label, spec in native_config['input_sources'].items():
        folder = Path(spec['inputs']); rp = folder / 'receipt.json'; ir = json.loads(rp.read_text()); ip = json.loads(Path(spec['readback']).read_text())
        expected_status = 'complete_duplication_alignment_input_materialization' if label == 'primary' else 'complete_additional_reference_alignment_input_materialization'
        expected_proof = 'passed_full_duplication_alignment_input_readback' if label == 'primary' else 'passed_full_reference_alignment_input_readback'
        proof_field = 'source_receipt_sha256' if label == 'primary' else 'producer_receipt_sha256'
        assert ir['status'] == expected_status and ip['status'] == expected_proof and ip[proof_field] == sha(rp)
        assert ir['plan_sha256'] == sha(spec['input_plan']); artifacts(bindings, folder, ir)
        bind(bindings, spec['readback']); bind(bindings, spec['input_plan']); counts = Counter(); n = 0
        with (folder / 'inputs.jsonl').open() as handle:
            for line in handle:
                row = json.loads(line); key = row['model_id'], row['version'], row['mask']
                assert isinstance(key[1], int) and key not in seen and key[2] in ['full', 'plddt70']; seen.add(key)
                counts[key[2] + ':' + row['status']] += 1; n += 1; row_count += 1
                if key[:2] not in needed: continue
                m = catalog[key[:2]]; positions = row['original_positions']
                assert row['original_length'] == m['length'] and row['source_sha256'] == m['sha256']
                assert positions == sorted(set(positions)) and len(positions) == len(row['sequence']) == row['retained_residues']
                assert all(isinstance(v, int) and 1 <= v <= row['original_length'] for v in positions)
                if row['status'] == 'ready': assert len(positions) >= 3 and row['sha256']
                elif row['status'] == 'too_few_retained_residues': assert len(positions) < 3
                else: assert row['status'] == 'source_rejected'
                if key[2] == 'full' and row['status'] != 'source_rejected':
                    assert positions == list(range(1, m['length'] + 1))
                    assert hashlib.sha256(row['sequence'].encode()).hexdigest() == m['sequence_sha256']
                inputs[key] = row
        assert n == ir['input_dispositions'] == 2 * ir['models'] and dict(counts) == ir['counts']
        source_counts[label] = dict(input_dispositions=n, counts=dict(counts))
    assert row_count == len(seen) == plan['expected']['all_input_dispositions']
    assert set(inputs) == {(*model, mask) for model in needed for mask in ['full', 'plddt70']}
    for model in needed:
        full, masked = inputs[(*model, 'full')], inputs[(*model, 'plddt70')]
        if full['status'] == 'source_rejected': assert masked['status'] == 'source_rejected'
        else: assert masked['sequence'] == ''.join(full['sequence'][p - 1] for p in masked['original_positions'])
    for p in [plan['triad_plan'], plan['triad_readback'], plan['triad_completion'], config['context_plan'], plan['reference_native_plan'], plan['primary_native_plan']]: bind(bindings, p)
    verify(bindings)
    return triads, inputs, r, source_counts, bindings


def load_measured_edges(plan, triads, bindings):
    pc = json.loads(Path(plan['primary_coverage_plan']).read_text()); pr = Path(pc['output'])
    cr = json.loads((pr / 'receipt.json').read_text()); ca = json.loads(Path(plan['primary_coverage_readback']).read_text()); cc = json.loads(Path(plan['primary_coverage_completion']).read_text())
    assert cr['status'] == 'complete_full_expanded_duplication_coverage_pending_independent_readback'
    assert ca['status'] == 'passed_full_expanded_pair_event_and_taxon_coverage_readback'
    assert cc['status'] == 'complete_verified_full_expanded_pair_event_and_taxon_coverage_screens' and len(cc['services']) == 2
    assert ca['producer_receipt_sha256'] == sha(pr / 'receipt.json') == cc['source_hashes'][str(pr / 'receipt.json')]
    assert cr['plan_sha256'] == ca['plan_sha256'] == sha(plan['primary_coverage_plan'])
    assert cc['source_hashes'][plan['primary_coverage_readback']] == sha(plan['primary_coverage_readback'])
    assert cr['pairs'] == cc['summary']['pairs'] == plan['expected']['primary_pairs']
    assert ca['pair_mask_rows'] == 2 * cr['pairs']
    artifacts(bindings, pr, cr)
    # Check a previously closed complete primary lineage without copying its
    # raw checkpoint hashes. Every checkpoint actually used is bound below.
    for p, h in cc['source_hashes'].items(): bind(bindings, p, h)
    rc = json.loads(Path(plan['reference_coverage_plan']).read_text()); rr = Path(rc['output'])
    r = json.loads((rr / 'receipt.json').read_text()); c = json.loads(Path(plan['reference_coverage_completion']).read_text())
    assert r['status'] == 'complete_full_reference_order_and_original_coverage_pending_independent_readback'
    assert c['status'] == 'complete_verified_full_reference_order_and_original_coverage' and c['exact_process_journals_checked'] == 2
    a = json.loads(Path(c['independent_readback']).read_text()); archive_path = c['full_hash_archive']
    bind(bindings, archive_path, c['full_hash_archive_sha256']); archive = json.loads(Path(archive_path).read_text())
    assert a['status'] == 'passed_full_reference_order_and_original_coverage_readback' and len(archive['services']) == 2
    assert r['plan_sha256'] == a['plan_sha256'] == sha(plan['reference_coverage_plan'])
    assert c['producer_receipt_sha256'] == a['producer_receipt_sha256'] == sha(rr / 'receipt.json')
    assert r['full_pairs'] == a['full_pairs'] == c['full_pairs'] == plan['expected']['reference_pairs']
    artifacts(bindings, rr, r); bind(bindings, c['independent_readback'], c['independent_readback_sha256'])
    assert archive['source_hashes'][str(rr / 'pair_mask_order_coverage.tsv')] == r['artifacts']['pair_mask_order_coverage.tsv']
    assert pc['screens'] == rc['screens'] == plan['screens'] and len(plan['screens']) == 6
    primary_config = json.loads(Path(plan['primary_native_plan']).read_text()); reference_config = json.loads(Path(plan['reference_native_plan']).read_text())
    source_config = json.loads(Path(json.loads(Path(plan['triad_plan']).read_text())['context_plan']).read_text())
    assert pc['queue'] == primary_config['queue'] == source_config['queue']
    union_config = json.loads(Path(rc['union_plan']).read_text()); assert union_config['native_plan'] == plan['reference_native_plan']
    assert reference_config['inventory'] == source_config['inventory']
    union_root = Path(union_config['output']); ur = json.loads((union_root / 'receipt.json').read_text())
    assert ur['status'] == 'complete_full_reference_measurement_union_pending_independent_readback'
    artifacts(bindings, union_root, ur)
    union_closure = json.loads(Path(rc['union_completion']).read_text())
    assert union_closure['status'] == 'complete_verified_full_reference_measurement_union'
    assert union_closure['producer_receipt_sha256'] == sha(union_root / 'receipt.json')
    assert archive['source_hashes'][rc['union_completion']] == sha(rc['union_completion'])
    assert archive['source_hashes'][str(union_root / 'receipt.json')] == sha(union_root / 'receipt.json')
    assert archive['source_hashes'][str(union_root / 'reference_measurement_dispositions.jsonl')] == ur['artifacts']['reference_measurement_dispositions.jsonl']
    bind(bindings, union_closure['full_hash_archive'], union_closure['full_hash_archive_sha256'])
    union_archive = json.loads(Path(union_closure['full_hash_archive']).read_text())
    assert union_archive['source_hashes'][str(union_root / 'reference_measurement_dispositions.jsonl')] == ur['artifacts']['reference_measurement_dispositions.jsonl']
    required = {kind: {t['edges'][label]['pair_key'] for t in triads for label in labels} for kind, labels in [('primary', ['ab']), ('reference', ['ar', 'br'])]}
    coverage = {}
    for kind, path, n in [('primary', pr / 'pair_mask_coverage.tsv', plan['expected']['primary_pairs']), ('reference', rr / 'pair_mask_order_coverage.tsv', plan['expected']['reference_pairs'])]:
        seen = set()
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                key = row['pair_key'], row['mask']; assert key not in seen; seen.add(key)
                if key[0] in required[kind]: coverage[(kind, *key)] = row
        assert len(seen) == 2 * n and seen == {(k, mask) for k, _ in seen for mask in ['full', 'plddt70']}
    native_root = Path(primary_config['output']); nr = json.loads((native_root / 'receipt.json').read_text())
    assert nr['status'] == 'complete_duplication_alignment_dispositions_pending_readback'
    assert nr['plan_sha256'] == sha(plan['primary_native_plan']); artifacts(bindings, native_root, nr)
    source_summary = json.loads((Path(pc['summary']) / 'receipt.json').read_text()); assert source_summary['native_receipt_sha256'] == sha(native_root / 'receipt.json')
    checkpoints = {}; states = set(); primary_manifest_sha = sha(Path(primary_config['inputs']) / 'inputs.jsonl')
    with (native_root / 'checkpoint_manifest.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            pair, mask, order = Path(row['path']).stem.split('-'); key = pair, mask, int(order)
            assert key not in states; states.add(key)
            if pair in required['primary']:
                cov = coverage[('primary', pair, mask)]; assert row['status'] == cov[f'order{order}_native_status']
                checkpoints[('primary', *key)] = dict(source_checkpoint=str(native_root / row['path']), source_checkpoint_sha256=row['sha256'], source_order=int(order), source_native_status=row['status'],
                                                       numerical_usable=cov[f'order{order}_status'] == 'aligned', numerical_exclusion_reasons=cov[f'order{order}_numerical_exclusion_reasons'].split(';') if cov[f'order{order}_numerical_exclusion_reasons'] else [],
                                                       selected_source='primary_native', source_plan_sha256=nr['plan_sha256'], source_input_manifest_sha256=primary_manifest_sha)
    assert len(states) == nr['directed_dispositions'] == 4 * plan['expected']['primary_pairs']
    states = set()
    with (union_root / 'reference_measurement_dispositions.jsonl').open() as handle:
        for line in handle:
            row = json.loads(line); key = row['pair_key'], row['mask'], row['order']; assert key not in states; states.add(key)
            if key[0] in required['reference']: checkpoints[('reference', *key)] = row
    assert len(states) == ur['directed_dispositions'] == 4 * plan['expected']['reference_pairs']
    assert set(coverage) == {(kind, pair, mask) for kind, pairs in required.items() for pair in pairs for mask in ['full', 'plddt70']}
    assert set(checkpoints) == {(kind, pair, mask, order) for kind, pairs in required.items() for pair in pairs for mask in ['full', 'plddt70'] for order in [0, 1]}
    for p in [plan['primary_coverage_plan'], plan['primary_coverage_readback'], plan['primary_coverage_completion'], plan['reference_coverage_plan'], plan['reference_coverage_completion'], rc['union_plan'], rc['union_completion']]: bind(bindings, p)
    verify(bindings)
    return checkpoints, coverage
