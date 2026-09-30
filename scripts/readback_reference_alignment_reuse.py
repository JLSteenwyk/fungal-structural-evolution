#!/usr/bin/env python3
"""Independently reconstruct every full reference reuse state from original inputs/results."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

from qualify_reference_alignment_reuse import load_sources
from screen_duplication_alignment_reuse import sha
from run_cross_clan_alignments import parse_output


def signature(row):
    names = ('model_id', 'version', 'mask', 'status', 'source_sha256', 'sequence',
             'original_positions', 'retained_residues', 'original_length', 'reason', 'sha256')
    return {k: row[k] for k in names if k in row}


def identity_hash(row):
    return hashlib.sha256(json.dumps(signature(row), sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def native_fields(source, pair, mask, order, desired):
    key = pair, mask, order
    proof = source['proofs'][key]
    path = source['root'] / proof['path']
    if sha(path) != proof['sha256'] or proof['path'] != f'pairs/{pair[:2]}/{pair}-{mask}-{order}.json':
        raise ValueError('Invalid original checkpoint path/hash')
    native = json.loads(path.read_text())
    rows = [source['inputs'][(*end, mask)] for end in desired]
    identities = [dict(model_id=end[0], version=end[1], status=row['status'], path=row.get('path'),
                       sha256=row.get('sha256')) for end, row in zip(desired, rows)]
    ready = all(r['status'] == 'ready' for r in rows)
    command = [source['plan']['usalign'], *[r['path'] for r in rows], *source['plan']['options']] if ready else None
    expected = dict(pair_key=pair, mask=mask, order=order, plan_sha256=source['plan_sha256'],
                    input_manifest_sha256=source['bundle'], inputs=identities, command=command,
                    status=proof['status'])
    if any(native[k] != value for k, value in expected.items()):
        raise ValueError('Original checkpoint binding differs')
    numeric, geometry = source['numeric'].get(key), source['geometry'].get(key)
    exclusions = []
    state = native['status']
    if state == 'aligned':
        if (not ready or native['returncode'] != 0 or numeric is None or geometry is None
                or parse_output(native['stdout'], [r['sequence'] for r in rows]) != native['metrics']
                or numeric['aligned_length'] != geometry['aligned_length']
                or numeric['rmsd_status'] != geometry['rmsd_status']):
            raise ValueError('Invalid original successful disposition')
        if numeric['rmsd_status'] == 'outside_printed_rounding':
            exclusions.append('rmsd_discrepancy')
        elif numeric['rmsd_status'] != 'within_printed_rounding':
            raise ValueError('Unknown original RMSD flag')
        if int(numeric['aligned_length']) < 3:
            exclusions.append('fewer_than_three_pairs')
        if geometry['geometry_status'] == 'degenerate_at_numeric_tolerance':
            exclusions.append('nonunique_rotation')
        elif geometry['geometry_status'] != 'unique_at_numeric_tolerance':
            raise ValueError('Unknown original rotation flag')
    else:
        if numeric is not None or geometry is not None:
            raise ValueError('Nonalignment has numeric rows')
        if state == 'input_unavailable':
            if ready or native['input_statuses'] != [r['status'] for r in rows]:
                raise ValueError('Invalid original short/rejected inputs')
        elif state == 'native_error':
            if not ready or native['returncode'] == 0:
                raise ValueError('Invalid original native failure')
        elif state == 'parse_error':
            if not ready or native['returncode'] != 0 or not native.get('error'):
                raise ValueError('Invalid original parse failure')
        elif state == 'timeout':
            if not ready or native['timeout_seconds'] != source['plan']['per_pair_timeout_seconds']:
                raise ValueError('Invalid original timeout')
        else:
            raise ValueError('Unknown original disposition')
        exclusions.append(state)
    return dict(source_checkpoint=str(path), source_checkpoint_sha256=proof['sha256'], source_order=order,
                source_native_status=state, source_numeric=numeric, source_geometry=geometry,
                numerical_exclusion_reasons=exclusions, numerical_usable=len(exclusions) == 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    root = Path(plan['output'])
    rp = root / 'receipt.json'
    record = json.loads(rp.read_text())
    if (record['status'] != 'complete_full_reference_reuse_qualification_pending_independent_readback'
            or record['plan_sha256'] != sha(args.plan)):
        raise ValueError('Incomplete or unbound producer')
    bindings = {str(args.plan): sha(args.plan), **plan['pins'], **record['source_hashes'], str(rp): sha(rp)}
    bindings.update({str(root / name): value for name, value in record['artifacts'].items()})
    def verify():
        for path, value in bindings.items():
            if sha(path) != value:
                raise ValueError('Changed original source/output: ' + path)
    verify()
    full, loader_choices, current, curmodels, sources, loaded_bindings = load_sources(plan)
    for path, value in loaded_bindings.items():
        if bindings.get(path) != value:
            raise ValueError('Missing/inconsistent loader source binding')
    # Independent source selection from the original outer-union catalog table.
    chosen = {}
    for row in csv.DictReader((Path(plan['reuse_candidates']) / 'pair_reuse_candidates.tsv').open(), delimiter='\t'):
        if 'reference_expanded' not in json.loads(row['new_sources']):
            continue
        previous = set(json.loads(row['matching_old_sources']))
        choice = ('primary_expanded_completed' if 'primary_expanded_completed' in previous else
                  'reference_old' if 'reference_old' in previous else '')
        if previous and not choice:
            raise ValueError('No audited native source')
        if row['pair_key'] in chosen:
            raise ValueError('Repeated catalog pair')
        chosen[row['pair_key']] = choice
    if chosen != loader_choices or set(chosen) != set(full):
        raise ValueError('Selected source/complete pair universe differs')
    target = json.loads(Path(plan['target_alignment_plan']).read_text())
    for source in sources.values():
        p = source['plan']
        h = sha(p['usalign'])
        if h != sha(target['usalign']) or p['pins'][p['usalign']] != h or p['options'] != target['options']:
            raise ValueError('Changed executable/options')
    model_checks = {}
    for pair, row in full.items():
        owner = chosen[pair]
        if not owner:
            continue
        for end in [(row['model_a'], int(row['version_a'])), (row['model_b'], int(row['version_b']))]:
            for mask in ['full', 'plddt70']:
                key = owner, *end, mask
                if key in model_checks:
                    continue
                a, b = current[(*end, mask)], sources[owner]['inputs'][(*end, mask)]
                m, old = curmodels[end], sources[owner]['models'][end]
                same_model = (m['sha256'], m['sequence_sha256'], m['length']) == (old['sha256'], old['sequence_sha256'], old['length'])
                equal = same_model and signature(a) == signature(b)
                for value in [m, old]:
                    if sha(value['path']) != value['sha256'] or bindings.get(value['path']) != value['sha256']:
                        raise ValueError('Raw coordinate bytes lack exact verified binding')
                for value in [a, b]:
                    if value['source_sha256'] != m['sha256']:
                        raise ValueError('Wrong materialized source')
                    if value['status'] == 'ready' and (sha(value['path']) != value['sha256'] or bindings.get(value['path']) != value['sha256']):
                        raise ValueError('PDB bytes lack exact verified binding')
                model_checks[key] = dict(source=owner, model_id=end[0], version=end[1], mask=mask,
                                         compatible=equal, current_projection_sha256=identity_hash(a),
                                         old_projection_sha256=identity_hash(b), current_raw_path=m['path'],
                                         old_raw_path=old['path'], raw_sha256=m['sha256'],
                                         current_pdb_path=a.get('path'), old_pdb_path=b.get('path'))
    with (root / 'model_input_identity_checks.jsonl').open() as handle:
        for key, expected in sorted(model_checks.items()):
            if json.loads(next(handle)) != expected:
                raise ValueError('Actual model-mask input identity export differs')
        if handle.read():
            raise ValueError('Unexpected model-mask input identity rows')
    counts, numerical, owners = Counter(), Counter(), Counter()
    checked = 0
    with (root / 'reference_reuse_dispositions.jsonl').open() as handle:
        for pair, original in sorted(full.items()):
            ends = [(original['model_a'], int(original['version_a'])), (original['model_b'], int(original['version_b']))]
            owner = chosen[pair]
            for mask in ['full', 'plddt70']:
                compatible = bool(owner) and all(model_checks[(owner, *e, mask)]['compatible'] for e in ends)
                for order in [0, 1]:
                    row = json.loads(next(handle))
                    expected = dict(pair_key=pair, model_a=ends[0][0], version_a=ends[0][1],
                                    model_b=ends[1][0], version_b=ends[1][1], mask=mask, order=order,
                                    selected_source=owner, numerical_usable=False, source_checkpoint=None,
                                    source_checkpoint_sha256=None, source_order=None, source_native_status=None,
                                    source_numeric=None, source_geometry=None, numerical_exclusion_reasons=[])
                    if owner == '':
                        expected['reuse_status'] = 'new_native_measurement_pending'
                    elif compatible is False:
                        expected['reuse_status'] = 'incompatible_inputs_require_new_measurement'
                    else:
                        source = sources[owner]
                        wanted_ends = ends if order == 0 else list(reversed(ends))
                        mappings = [source['pairs'][pair], list(reversed(source['pairs'][pair]))]
                        if mappings.count(wanted_ends) != 1:
                            raise ValueError('Ambiguous directed endpoint correspondence')
                        old_order = mappings.index(wanted_ends)
                        expected.update(native_fields(source, pair, mask, old_order, wanted_ends))
                        expected['reuse_status'] = 'verified_identical_input_checkpoint_and_retained_disposition'
                        if bindings.get(expected['source_checkpoint']) != expected['source_checkpoint_sha256']:
                            raise ValueError('Unbound reused native checkpoint')
                        numerical['usable' if expected['numerical_usable'] else ';'.join(expected['numerical_exclusion_reasons'])] += 1
                    if row != expected:
                        raise ValueError('Actual full directed reuse export differs: ' + pair)
                    counts[row['reuse_status']] += 1
                    owners[owner or 'no_old_source'] += 1
                    checked += 1
                    if checked % 10000 == 0:
                        print('Independently checked directed reuse states', checked, flush=True)
        if handle.read():
            raise ValueError('Unexpected directed reuse rows')
    if (checked != 4 * len(full) or record['full_reference_pairs'] != len(full)
            or record['directed_dispositions'] != checked or record['counts'] != dict(counts)
            or record['numerical_counts'] != dict(numerical)
            or record['selected_source_dispositions'] != dict(owners)
            or record['unique_model_mask_source_checks'] != len(model_checks)):
        raise ValueError('Full reuse totals differ')
    verify()
    result = dict(status='passed_full_reference_input_checkpoint_numeric_reuse_readback',
                  producer_receipt_sha256=sha(rp), full_reference_pairs=len(full), directed_dispositions=checked,
                  unique_model_mask_source_checks=len(model_checks), counts=dict(counts), numerical_counts=dict(numerical),
                  selected_source_dispositions=dict(owners), source_hashes=bindings, checker_sha256=sha(__file__),
                  scientific_eligibility=False, scope='Every actual model-mask identity and full reference pair/mask/order '
                  'export independently reconstructed from original catalogs, inputs, checkpoints and complete '
                  'native numerical/geometry proofs. Independent source preference, input projection equality, '
                  'directed correspondence and numerical exclusions; source I/O/proof loader and native-text parser '
                  'shared with existing pipeline. All raw/PDB/checkpoint bytes checked, no flags cleared or native '
                  'results changed. Whole result union with new measurements and biological controls remain.')
    with args.output.open('x') as handle:
        handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['status', 'full_reference_pairs', 'directed_dispositions', 'counts', 'numerical_counts']}), flush=True)


if __name__ == '__main__':
    main()
