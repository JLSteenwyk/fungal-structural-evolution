#!/usr/bin/env python3
"""Independently rebuild every full-reference union field, source, direction and flag."""
import argparse
import json
from collections import Counter
from pathlib import Path
from reference_measurement_union_sources import load_union_sources, bind, verify
from run_ortholog_pair_guide_comparison import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    plan = json.loads(args.plan.read_text())
    original, native_source, bindings = load_union_sources(plan)
    root = Path(plan['output'])
    receipt = json.loads((root / 'receipt.json').read_text())
    assert receipt['status'] == 'complete_full_reference_measurement_union_pending_independent_readback'
    assert receipt['plan_sha256'] == sha(args.plan)
    bind(bindings, args.plan)
    bind(bindings, root / 'receipt.json')
    for name, digest in receipt['artifacts'].items(): bind(bindings, root / name, digest)
    verify(bindings)
    seen, source_counts, statuses, numeric_counts = set(), Counter(), Counter(), Counter()
    with original.open() as originals, (root / 'reference_measurement_dispositions.jsonl').open() as exported:
        for line in originals:
            source = json.loads(line)
            row = json.loads(exported.readline())
            pair, mask, order = source['pair_key'], source['mask'], source['order']
            key = pair, mask, order
            assert key not in seen and pair in native_source['full'] and mask in ['full', 'plddt70'] and order in [0, 1]
            models = native_source['full'][pair]
            assert models == [(source['model_a'], source['version_a']), (source['model_b'], source['version_b'])]
            if pair in native_source['new']:
                assert source['reuse_status'] == 'new_native_measurement_pending' and source['selected_source'] == ''
                proof = native_source['checkpoints'][key]
                label, path, digest, native_order = 'reference_new_native', proof['path'], proof['sha256'], order
                numeric = native_source['numeric'].get(key)
                geometry = native_source['geometry'].get(key)
            else:
                assert source['reuse_status'] == 'verified_identical_input_checkpoint_and_retained_disposition'
                label, path, digest, native_order = [source[k] for k in ['selected_source', 'source_checkpoint', 'source_checkpoint_sha256', 'source_order']]
                numeric, geometry = source['source_numeric'], source['source_geometry']
            assert native_order in [0, 1] and sha(path) == digest
            bind(bindings, path, digest)
            raw = json.loads(Path(path).read_text())
            assert raw['pair_key'] == pair and raw['mask'] == mask and raw['order'] == native_order
            target = models if order == 0 else list(reversed(models))
            assert [(v['model_id'], v['version']) for v in raw['inputs']] == target
            status = raw['status']
            if label == 'reference_new_native':
                assert status == proof['status'] and raw['plan_sha256'] == sha(plan['native_plan'])
                assert raw['input_manifest_sha256'] == native_source['native_receipt']['input_bundle_sha256']
            else:
                assert label in ['primary_expanded_completed', 'reference_old'] and status == source['source_native_status']
            if status == 'aligned':
                assert raw['returncode'] == 0 and numeric is not None and geometry is not None
                assert numeric['pair_key'] == geometry['pair_key'] == pair and numeric['mask'] == geometry['mask'] == mask
                assert int(numeric['order']) == int(geometry['order']) == native_order
                assert numeric['aligned_length'] == geometry['aligned_length'] and numeric['rmsd_status'] == geometry['rmsd_status']
                conditions = [('rmsd_discrepancy', numeric['rmsd_status'] != 'within_printed_rounding'),
                              ('fewer_than_three_pairs', int(numeric['aligned_length']) < 3),
                              ('nonunique_rotation', geometry['geometry_status'] != 'unique_at_numeric_tolerance')]
                why = [name for name, excluded in conditions if excluded]
                metrics = raw['metrics']
                assert metrics['aligned_length'] == int(numeric['aligned_length'])
            else:
                assert status in ['input_unavailable', 'native_error', 'parse_error', 'timeout']
                assert numeric is None and geometry is None and 'metrics' not in raw
                why, metrics = [status], None
            if label != 'reference_new_native':
                assert why == source['numerical_exclusion_reasons'] and (len(why) == 0) == source['numerical_usable']
            expected = dict(pair_key=pair, model_a=models[0][0], version_a=models[0][1], model_b=models[1][0], version_b=models[1][1],
                            mask=mask, order=order, directed_endpoints=[dict(model_id=v['model_id'], version=v['version']) for v in raw['inputs']],
                            selected_source=label, source_checkpoint=path, source_checkpoint_sha256=digest, source_order=native_order,
                            source_plan_sha256=raw['plan_sha256'], source_input_manifest_sha256=raw['input_manifest_sha256'],
                            source_native_status=status, source_numeric=numeric, source_geometry=geometry, native_metrics=metrics,
                            numerical_exclusion_reasons=why, numerical_usable=len(why) == 0)
            assert row == expected
            source_counts[label] += 1; statuses[mask + ':' + status] += 1
            numeric_counts['usable' if not why else ';'.join(why)] += 1
            seen.add(key)
            if len(seen) % 10000 == 0: print('Union readback', len(seen), '/', 4 * plan['full_pairs'], flush=True)
        assert exported.readline() == ''
    required = {(p, m, o) for p in native_source['full'] for m in ['full', 'plddt70'] for o in range(2)}
    assert seen == required and receipt['directed_dispositions'] == len(seen) == 4 * plan['full_pairs']
    assert receipt['full_pairs'] == plan['full_pairs'] and receipt['source_dispositions'] == dict(source_counts)
    assert receipt['native_status_counts'] == dict(statuses) and receipt['numerical_counts'] == dict(numeric_counts)
    # Union bindings retain every old/new checkpoint, not only the currently usable subset.
    for path, digest in bindings.items():
        if path not in [str(root / 'receipt.json'), *(str(root / k) for k in receipt['artifacts'])]:
            assert receipt['source_hashes'].get(path) == digest
    verify(bindings)
    result = dict(status='passed_full_reference_measurement_union_readback', plan_sha256=sha(args.plan),
                  producer_receipt_sha256=sha(root / 'receipt.json'), full_pairs=plan['full_pairs'], directed_dispositions=len(seen),
                  source_dispositions=dict(source_counts), native_status_counts=dict(statuses), numerical_counts=dict(numeric_counts),
                  source_hashes=bindings, checker_sha256=sha(__file__), scientific_eligibility=False,
                  scope='Every full reference pair/mask/order independently rebuilt from the verified complete reuse ledger '
                        'and newly completed native source, including actual checkpoint byte hashes, ordered model/version '
                        'endpoints, original native metrics, every numerical/geometry field and all exclusions. '
                        'Source I/O/proof loader shared, producer projection/exclusion logic not shared. No native '
                        'optimization repeated, scientific filtering or biological effect inferred.')
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True)


if __name__ == '__main__': main()
