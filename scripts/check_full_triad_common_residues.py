#!/usr/bin/env python3
"""Check mapping algorithms, masked positions, failed orders and false-export rejection."""
import contextlib
import gzip
import hashlib
import io
import json
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch
import map_full_triad_common_residues as producer
import readback_full_triad_common_residues as reader
from run_ortholog_pair_guide_comparison import sha


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def digest(value):
    return hashlib.sha256(json.dumps(value, separators=(',', ':')).encode()).hexdigest()


def main():
    with tempfile.TemporaryDirectory(prefix='full-triad-common-native-fixture-') as temp:
        root = Path(temp); am, bm, rm, cm = ('A', 6), ('B', 10), ('R', 10), ('C', 6)
        inputs = {}
        for model in [am, bm, rm, cm]:
            for mask in ['full', 'plddt70']:
                selected = list(range(1, 9)) if mask == 'full' else {'A': [2, 4, 7], 'B': [1, 4, 8], 'R': [3, 5, 6], 'C': [2, 4]}[model[0]]
                sequence = ''.join('ACDEFGHI'[p - 1] for p in selected)
                inputs[(*model, mask)] = dict(model_id=model[0], version=model[1], mask=mask, original_length=8, original_positions=selected,
                                             sequence=sequence, retained_residues=len(selected), status='ready' if len(selected) >= 3 else 'too_few_retained_residues', sha256=digest([model, mask]))
        triads = []
        for ref in [rm, cm]:
            models = [am, bm, ref]; triads.append(dict(triad_id=digest(models), models=[list(m) for m in models],
                   source_design_ready_links=1, edges={label: dict(pair_key=digest(sorted(ends))) for label, ends in [('ab', [am, bm]), ('ar', [am, ref]), ('br', [bm, ref])]}))
        screens = [dict(id=sid, minimum_aligned_residues=n, minimum_original_coverage=c) for n in [30, 50] for c, suffix in [(.5, '50'), (.7, '70'), (.9, '90')] for sid in [f'n{n}_c{suffix}']]
        coverage, checkpoints = {}, {}
        for triad in triads:
            ref = tuple(triad['models'][2])
            for label, kind, ends in [('ab', 'primary', [bm, am]), ('ar', 'reference', [am, ref]), ('br', 'reference', [ref, bm])]:
                pair = triad['edges'][label]['pair_key']
                for mask in ['full', 'plddt70']:
                    covkey = kind, pair, mask
                    if covkey in coverage: continue
                    cov = dict(model_a=ends[0][0], version_a=str(ends[0][1]), model_b=ends[1][0], version_b=str(ends[1][1]))
                    for order in [0, 1]:
                        directed = ends if order == 0 else ends[::-1]; rows = [inputs[(*model, mask)] for model in directed]
                        status = 'aligned'; reasons = []
                        if mask == 'plddt70':
                            if any(r['status'] != 'ready' for r in rows): status = 'input_unavailable'
                            elif label == 'ab' and order == 0: reasons = ['rmsd_discrepancy', 'nonunique_rotation']
                            elif label == 'ab': status = 'native_error'
                            elif label == 'ar': status = 'parse_error' if order == 0 else 'timeout'
                        strings = [r['sequence'] for r in rows]
                        if mask == 'full' and label == 'ab' and order == 1: strings = [strings[0] + '-', '-' + strings[1]]
                        metrics = dict(alignment_left=strings[0], alignment_right=strings[1], aligned_length=sum(a != '-' and b != '-' for a, b in zip(*strings)))
                        original_order = 1 - order if label == 'ar' else order
                        path = root / f'{kind}-{pair}-{mask}-{order}.json'
                        raw = dict(pair_key=pair, mask=mask, order=original_order, status=status, plan_sha256='original-plan', input_manifest_sha256='original-manifest',
                                   inputs=[{k: r[k] for k in ['model_id', 'version', 'status', 'sha256']} for r in rows])
                        if status == 'aligned': raw['metrics'] = metrics
                        save(path, raw)
                        checkpoints[(kind, pair, mask, order)] = dict(source_checkpoint=str(path), source_checkpoint_sha256=sha(path), source_order=original_order, source_native_status=status,
                                                numerical_usable=status == 'aligned' and not reasons, numerical_exclusion_reasons=reasons, selected_source='reversed_old_reference' if label == 'ar' else 'primary_native',
                                                source_plan_sha256='original-plan', source_input_manifest_sha256='original-manifest')
                        cov[f'order{order}_native_status'] = status; cov[f'order{order}_status'] = 'excluded_numerically' if reasons else status
                        cov[f'order{order}_numerical_exclusion_reasons'] = ';'.join(reasons)
                    for s in screens:
                        sid = s['id']; cov[sid + '_pass'] = '1' if mask == 'full' else '0'; cov[sid + '_exclusions'] = '' if mask == 'full' else 'fixture_both_order_exclusion'
                    coverage[covkey] = cov
        completion = root / 'design-complete.json'; save(completion, dict(synthetic_source_stub=True))
        out = root / 'out'; plan = root / 'plan.json'
        save(plan, dict(output=str(out), triad_completion=str(completion), screens=screens, expected=dict(potential_correspondence_states=32, all_input_dispositions=16),
                        resources=dict(minimum_free_disk_gib=0), scope='Synthetic mapping algorithms/native states only. Source proof loaders are injected; no actual source qualification.'))
        design = dict(target_contexts=5, reference_tie_records=7, duplicate_reference_links=14, unique_ordered_model_triads=3)
        counts = {'synthetic': dict(input_dispositions=16)}
        def load_design(_): return triads, inputs, design, counts, {str(completion): sha(completion)}
        def load_measured(_, __, ___): return checkpoints, coverage
        def run(module, output=None):
            argv = ['fixture', '--plan', str(plan)] + (['--output', str(output)] if output else [])
            with patch.object(module, 'load_design_inputs', load_design), patch.object(module, 'load_measured_edges', load_measured), patch.object(sys, 'argv', argv), contextlib.redirect_stdout(io.StringIO()):
                module.main()
        run(producer); run(reader, root / 'passed.json')
        path = out / 'common_residue_maps.jsonl.gz'
        with gzip.open(path, 'rt') as f: original = [json.loads(line) for line in f]
        receipt = json.loads((out / 'receipt.json').read_text()); assert len(original) == 32
        assert original[0]['reference_common_triples'] == [[i, i, i] for i in range(1, 9)]
        assert original[0]['cycle_consistent_triples'] == original[0]['reference_common_triples']
        assert original[4]['reference_common_triples'] == original[0]['reference_common_triples'] and original[4]['cycle_consistent_triples'] == []
        assert original[0]['edge_provenance'][1]['source_order'] != original[0]['orders'][1]
        p70_ab = checkpoints[('primary', triads[0]['edges']['ab']['pair_key'], 'plddt70', 0)]
        _, pairs, trace = producer.recover_native(p70_ab, coverage[('primary', triads[0]['edges']['ab']['pair_key'], 'plddt70')], inputs, triads[0]['edges']['ab']['pair_key'], 'plddt70', 0, {})
        assert pairs == [(1, 2), (4, 4), (8, 7)] and trace['mapping_exclusions'] == ['rmsd_discrepancy', 'nonunique_rotation']
        all_statuses = {t['source_native_status'] for row in original for t in row['edge_provenance']}
        assert all_statuses == {'aligned', 'input_unavailable', 'native_error', 'parse_error', 'timeout'}
        rejected = []
        for label in ['changed_original_position', 'promoted_cycle', 'cleared_numerical_flag', 'cleared_failed_order', 'changed_source_order', 'changed_directed_endpoint', 'favorable_coverage', 'changed_original_length', 'removed_state', 'duplicated_state']:
            altered = json.loads(json.dumps(original))
            if label == 'changed_original_position': altered[0]['reference_common_triples'][0][0] = 2
            elif label == 'promoted_cycle': altered[4]['cycle_consistent_triples'] = altered[4]['reference_common_triples']; altered[4]['cycle_consistent_count'] = altered[4]['common_reference_count']
            elif label == 'cleared_numerical_flag': altered[8]['edge_provenance'][0]['numerical_exclusion_reasons'] = []
            elif label == 'cleared_failed_order': altered[8]['edge_provenance'][1]['mapping_exclusions'] = []
            elif label == 'changed_source_order': altered[0]['edge_provenance'][1]['source_order'] = altered[0]['orders'][1]
            elif label == 'changed_directed_endpoint': altered[0]['edge_provenance'][1]['directed_endpoints'].reverse()
            elif label == 'favorable_coverage': altered[8]['all_three_pair_screen_pass'][screens[0]['id']] = True
            elif label == 'changed_original_length': altered[0]['original_lengths'][0] = 3
            elif label == 'removed_state': altered.pop()
            else: altered.append(altered[-1])
            with gzip.open(path, 'wt', compresslevel=1) as f: f.write(''.join(json.dumps(row) + '\n' for row in altered))
            save(out / 'receipt.json', dict(receipt, artifacts={path.name: sha(path)}))
            target = root / (label + '.json')
            try: run(reader, target)
            except (AssertionError, ValueError): pass
            else: raise AssertionError('False export accepted: ' + label)
            assert not target.exists(); rejected.append(label)
        print(json.dumps(dict(status='passed_full_synthetic_original_residue_mapping_checks', triads=2, mask_order_states=32, numerical_raw_map_pairs=3,
                             native_statuses=sorted(all_statuses), rejected_rehashed_exports=rejected,
                             scope='Synthetic raw checkpoint/masked input/projection fixtures with injected source-I/O loaders; no production source proofs or sampling pilot. Indexed producer and independent iterator reader both exercised on full grid, alignment gaps, unequal model versions, reversed original source order/current endpoints, nonconsecutive masked protein positions, raw numerically excluded maps and all four native failure statuses. Full source/proof loaders require separate actual input preflight and closed native production validation.'), indent=2))


if __name__ == '__main__': main()
