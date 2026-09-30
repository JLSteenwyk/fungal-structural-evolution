#!/usr/bin/env python3
"""Union every reference pair/mask/order with intact native provenance and exclusions."""
import argparse
import json
from collections import Counter
from pathlib import Path

from reference_measurement_union_sources import load_union_sources, bind, verify
from run_ortholog_pair_guide_comparison import sha


def exclusions(status, numeric, geometry):
    if status != 'aligned': return [status]
    why = []
    if numeric['rmsd_status'] != 'within_printed_rounding': why.append('rmsd_discrepancy')
    if int(numeric['aligned_length']) < 3: why.append('fewer_than_three_pairs')
    if geometry['geometry_status'] != 'unique_at_numeric_tolerance': why.append('nonunique_rotation')
    return why


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    reuse_path, new, bindings = load_union_sources(plan)
    bind(bindings, args.plan)
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    seen, sources, statuses, reasons = set(), Counter(), Counter(), Counter()
    with reuse_path.open() as source, (out / 'reference_measurement_dispositions.jsonl').open('x') as target:
        for line in source:
            row = json.loads(line)
            key = row['pair_key'], row['mask'], row['order']
            pair, mask, order = key
            assert key not in seen and pair in new['full'] and mask in ['full', 'plddt70'] and order in [0, 1]
            ends = [(row['model_a'], row['version_a']), (row['model_b'], row['version_b'])]
            assert ends == new['full'][pair]
            directed = ends if order == 0 else ends[::-1]
            if pair in new['new']:
                assert row['reuse_status'] == 'new_native_measurement_pending' and not row['selected_source']
                checkpoint = new['checkpoints'][key]
                path, digest, source_order = checkpoint['path'], checkpoint['sha256'], order
                label = 'reference_new_native'
                numeric, geometry = new['numeric'].get(key), new['geometry'].get(key)
            else:
                assert row['reuse_status'] == 'verified_identical_input_checkpoint_and_retained_disposition'
                assert row['selected_source'] in ['primary_expanded_completed', 'reference_old']
                path, digest, source_order = row['source_checkpoint'], row['source_checkpoint_sha256'], row['source_order']
                label = row['selected_source']
                numeric, geometry = row['source_numeric'], row['source_geometry']
            assert source_order in [0, 1] and sha(path) == digest
            bind(bindings, path, digest)
            native = json.loads(Path(path).read_text())
            assert (native['pair_key'], native['mask'], native['order']) == (pair, mask, source_order)
            assert [(r['model_id'], r['version']) for r in native['inputs']] == directed
            status = native['status']
            assert status in ['aligned', 'input_unavailable', 'native_error', 'parse_error', 'timeout']
            if pair in new['new']:
                assert status == checkpoint['status']
                assert native['plan_sha256'] == sha(plan['native_plan'])
                assert native['input_manifest_sha256'] == new['native_receipt']['input_bundle_sha256']
            else: assert status == row['source_native_status']
            if status == 'aligned':
                assert native['returncode'] == 0 and numeric is not None and geometry is not None
                assert (numeric['pair_key'], numeric['mask'], int(numeric['order'])) == (pair, mask, source_order)
                assert (geometry['pair_key'], geometry['mask'], int(geometry['order'])) == (pair, mask, source_order)
                assert numeric['aligned_length'] == geometry['aligned_length'] and numeric['rmsd_status'] == geometry['rmsd_status']
                assert int(numeric['aligned_length']) == native['metrics']['aligned_length']
            else: assert numeric is None and geometry is None and 'metrics' not in native
            why = exclusions(status, numeric, geometry)
            if pair not in new['new']:
                assert why == row['numerical_exclusion_reasons'] and (not why) == row['numerical_usable']
            exported = dict(pair_key=pair, model_a=ends[0][0], version_a=ends[0][1], model_b=ends[1][0], version_b=ends[1][1],
                            mask=mask, order=order, directed_endpoints=[dict(model_id=m, version=v) for m, v in directed],
                            selected_source=label, source_checkpoint=path, source_checkpoint_sha256=digest, source_order=source_order,
                            source_plan_sha256=native['plan_sha256'], source_input_manifest_sha256=native['input_manifest_sha256'],
                            source_native_status=status, source_numeric=numeric, source_geometry=geometry,
                            native_metrics=native.get('metrics'), numerical_exclusion_reasons=why, numerical_usable=not why)
            target.write(json.dumps(exported, separators=(',', ':')) + '\n')
            seen.add(key); sources[label] += 1; statuses[mask + ':' + status] += 1
            reasons['usable' if not why else ';'.join(why)] += 1
            if len(seen) % 10000 == 0: print('Union dispositions', len(seen), '/', 4 * plan['full_pairs'], flush=True)
    assert seen == {(p, m, o) for p in new['full'] for m in ['full', 'plddt70'] for o in [0, 1]}
    assert sources['reference_new_native'] == 4 * plan['new_pairs']
    verify(bindings)
    result = dict(status='complete_full_reference_measurement_union_pending_independent_readback',
                  plan_sha256=sha(args.plan), full_pairs=plan['full_pairs'], directed_dispositions=len(seen),
                  source_dispositions=dict(sources), native_status_counts=dict(statuses), numerical_counts=dict(reasons),
                  artifacts={'reference_measurement_dispositions.jsonl': sha(out / 'reference_measurement_dispositions.jsonl')},
                  source_hashes=bindings, scientific_eligibility=False,
                  scope='All reference pairs/two masks/two orders joined from independently verified identical-input old '
                        'sources and completed new native measurements. Original checkpoints/metrics/numeric/geometry '
                        'records and source order remain traceable; target order maps actual ordered model/version endpoints. '
                        'No outcome-based source switch, cleared numerical flag, retried exclusion or zero for missing/error '
                        'results. This numerical/provenance union is not coverage/confidence, biological orthology, '
                        'common-residue/sequence-locked/domain/PAE or calibrated asymmetry qualification.')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'artifacts']}, indent=2), flush=True)


if __name__ == '__main__': main()
