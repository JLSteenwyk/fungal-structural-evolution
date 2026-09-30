#!/usr/bin/env python3
"""Qualify exact old reference measurements without discarding missing or numerical states."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path

from screen_duplication_alignment_reuse import sha
from run_cross_clan_alignments import parse_output


MASKS = ['full', 'plddt70']


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def endpoints(row):
    return [(row['model_a'], int(row['version_a'])), (row['model_b'], int(row['version_b']))]


def projection(row):
    # Collection paths/shard numbers may differ; actual PDB content must not.
    fields = ['model_id', 'version', 'mask', 'status', 'source_sha256', 'sequence',
              'original_positions', 'retained_residues', 'original_length', 'reason', 'sha256']
    return {k: row[k] for k in fields if k in row}


def choose_source(candidates):
    for name in ['primary_expanded_completed', 'reference_old']:
        if name in candidates:
            return name
    if candidates:
        raise ValueError('No configured completed source for existing reference pair')
    return ''


def artifact_bindings(root, record):
    return {str(root / 'receipt.json'): sha(root / 'receipt.json'),
            **{str(root / name): value for name, value in record['artifacts'].items()}}


def load_inputs(specs, wanted, bindings, proofs=False):
    output = {}
    for spec in specs:
        root = Path(spec['inputs'])
        rp = root / 'receipt.json'
        receipt = json.loads(rp.read_text())
        ip = json.loads(Path(spec['input_plan']).read_text())
        if (receipt['status'] != spec['status'] or receipt['plan_sha256'] != sha(spec['input_plan'])
                or Path(ip['output']).resolve() != root.resolve()):
            raise ValueError('Unbound materialized inputs')
        bindings.update(artifact_bindings(root, receipt))
        bindings[spec['input_plan']] = sha(spec['input_plan'])
        # Current complete written-input proofs qualify the identical old bytes.
        if proofs:
            proof = json.loads(Path(spec['readback']).read_text())
            if (proof['status'] != spec['readback_status']
                    or proof[spec['readback_receipt_field']] != sha(rp)):
                raise ValueError('Missing full current written-input proof')
            bindings[spec['readback']] = sha(spec['readback'])
        coordinate = Path(ip['coordinates']) / 'receipt.json'
        native_reader = Path(ip['readback']) / 'receipt.json'
        cr, nr = [json.loads(p.read_text()) for p in [coordinate, native_reader]]
        if (receipt['coordinate_receipt_sha256'] != sha(coordinate)
                or receipt['readback_receipt_sha256'] != sha(native_reader)
                or nr['producer_receipt_sha256'] != sha(coordinate)
                or nr['status'] != 'passed_duplication_exported_ca_readback'):
            raise ValueError('Coordinate/native proof lineage differs')
        bindings.update({str(p): sha(p) for p in [coordinate, native_reader]})
        counts, seen = Counter(), set()
        with (root / 'inputs.jsonl').open() as handle:
            for line in handle:
                row = json.loads(line)
                key = row['model_id'], row['version'], row['mask']
                if key in seen or key[2] not in MASKS:
                    raise ValueError('Repeated/invalid materialized input')
                seen.add(key)
                counts[key[2] + ':' + row['status']] += 1
                if key[:2] in wanted:
                    if key in output:
                        raise ValueError('Overlapping materialized partitions')
                    output[key] = row
        if (len(seen) != receipt['input_dispositions'] or len(seen) != 2 * receipt['models']
                or dict(counts) != receipt['counts']):
            raise ValueError('Incomplete materialized disposition universe')
    expected = {(*key, mask) for key in wanted for mask in MASKS}
    if set(output) != expected:
        raise ValueError('Missing requested input masks')
    return output


def table_subset(path, wanted, fields, expected_count=None):
    result, count = {}, 0
    with Path(path).open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            count += 1
            if row['pair_key'] not in wanted:
                continue
            key = row['pair_key'], row['mask'], int(row['order'])
            if key in result:
                raise ValueError('Repeated selected source row')
            result[key] = row
    if expected_count is not None and count != expected_count:
        raise ValueError('Incomplete full source table')
    return result


def load_source(spec, wanted, bindings):
    sp = Path(spec['alignment_plan'])
    plan = json.loads(sp.read_text())
    root = Path(plan['output'])
    rp = root / 'receipt.json'
    receipt = json.loads(rp.read_text())
    if receipt['status'] != spec['alignment_status'] or receipt['plan_sha256'] != sha(sp):
        raise ValueError('Incomplete old native source')
    bindings.update(artifact_bindings(root, receipt))
    bindings[str(sp)] = sha(sp)
    queue = Path(plan.get('inventory', plan.get('queue')))
    qr = json.loads((queue / 'receipt.json').read_text())
    pairpath = queue / 'model_pairs.tsv'
    if sha(pairpath) != qr['artifacts']['model_pairs.tsv']:
        raise ValueError('Changed old pair source')
    all_pairs = {}
    for row in csv.DictReader(pairpath.open(), delimiter='\t'):
        if row['pair_key'] in all_pairs:
            raise ValueError('Repeated old pair')
        all_pairs[row['pair_key']] = endpoints(row)
    if any(key not in all_pairs for key in wanted):
        raise ValueError('Completed source does not contain selected pair')
    bindings.update({str(pairpath): sha(pairpath), str(queue / 'receipt.json'): sha(queue / 'receipt.json')})
    models = {}
    wanted_models = {key for pair in wanted for key in all_pairs[pair]}
    mp = queue / 'models.jsonl'
    if sha(mp) != qr['artifacts'][mp.name]:
        raise ValueError('Changed old model catalog')
    for line in mp.open():
        row = json.loads(line)
        key = row['model_id'], row['version']
        if key in wanted_models:
            if key in models:
                raise ValueError('Repeated selected old model')
            models[key] = row
    if set(models) != wanted_models:
        raise ValueError('Missing selected old model')
    bindings[str(mp)] = sha(mp)
    inputs = load_inputs(spec['input_sources'], wanted_models, bindings)
    mb = {str(Path(s['inputs']) / name): sha(Path(s['inputs']) / name)
          for s in spec['input_sources'] for name in ['receipt.json', 'inputs.jsonl']}
    bundle = (sha(Path(spec['input_sources'][0]['inputs']) / 'inputs.jsonl') if spec['kind'] == 'primary'
              else hashlib.sha256(json.dumps(mb, sort_keys=True).encode()).hexdigest())
    if spec['kind'] == 'primary':
        if receipt['input_manifest_sha256'] != bundle or receipt['input_receipt_sha256'] != mb[str(Path(spec['input_sources'][0]['inputs']) / 'receipt.json')]:
            raise ValueError('Old primary input lineage differs')
    elif receipt['input_bundle_sha256'] != bundle or receipt['input_bindings'] != mb:
        raise ValueError('Old reference input lineage differs')
    proofs = {}
    count = 0
    for row in csv.DictReader((root / 'checkpoint_manifest.tsv').open(), delimiter='\t'):
        count += 1
        name = Path(row['path']).stem
        pair, mask, order = name.rsplit('-', 2)
        if pair in wanted:
            key = pair, mask, int(order)
            if key in proofs:
                raise ValueError('Duplicate old directed checkpoint')
            proofs[key] = row
    if count != receipt['directed_dispositions'] or set(proofs) != {(p, m, o) for p in wanted for m in MASKS for o in [0, 1]}:
        raise ValueError('Incomplete old native directed source')
    dp = Path(spec['diagnostic'])
    gp = Path(spec['geometry'])
    dr, gr = [json.loads((p / 'receipt.json').read_text()) for p in [dp, gp]]
    proof = json.loads(Path(spec['geometry_readback']).read_text())
    if (dr['status'] != spec['diagnostic_status'] or dr['producer_receipt_sha256'] != sha(rp)
            or gr['status'] != spec['geometry_status'] or gr['diagnostic_receipt_sha256'] != sha(dp / 'receipt.json')
            or proof['status'] != spec['geometry_readback_status']
            or proof['producer_receipt_sha256'] != sha(gp / 'receipt.json')
            or proof['alignments_checked'] != gr['alignments']):
        raise ValueError('Old full numerical/geometry proof lineage differs')
    bindings.update(artifact_bindings(dp, dr))
    bindings.update(artifact_bindings(gp, gr))
    bindings[spec['geometry_readback']] = sha(spec['geometry_readback'])
    numeric = table_subset(dp / 'numeric_readback.tsv', wanted, None, dr['numerically_checked_alignments'])
    geometry = table_subset(gp / 'alignment_geometry.tsv', wanted, None, gr['alignments'])
    if set(numeric) != set(geometry):
        raise ValueError('Old numerical/geometry selected universes differ')
    return dict(plan=plan, plan_sha256=sha(sp), root=root, pairs=all_pairs, models=models, inputs=inputs,
                proofs=proofs, numeric=numeric, geometry=geometry, bundle=bundle)


def load_sources(plan):
    bindings = dict(plan['pins'])
    inventory = Path(plan['inventory'])
    rp = inventory / 'receipt.json'
    receipt = json.loads(rp.read_text())
    proof = json.loads(Path(plan['inventory_readback']).read_text())
    if (receipt['status'] != 'complete_provisional_reference_comparison_inventory'
            or proof['status'] != 'passed_full_reference_comparison_ledger_and_native_model_readback'
            or proof['producer_receipt_sha256'] != sha(rp)):
        raise ValueError('Current full reference ledger lacks native proof')
    bindings.update(artifact_bindings(inventory, receipt))
    rows = list(csv.DictReader((inventory / 'model_pairs.tsv').open(), delimiter='\t'))
    if len(rows) != receipt['unique_distinct_model_pairs']:
        raise ValueError('Incomplete current reference pairs')
    full, selected = {}, {}
    wanted = {key for row in rows for key in endpoints(row)}
    curmodels = {}
    for line in (inventory / 'models.jsonl').open():
        row = json.loads(line)
        key = row['model_id'], row['version']
        if key in wanted:
            if key in curmodels:
                raise ValueError('Repeated current model')
            curmodels[key] = row
    if set(curmodels) != wanted:
        raise ValueError('Incomplete current endpoint catalog')
    root = Path(plan['reuse_candidates'])
    rr = json.loads((root / 'receipt.json').read_text())
    readback = json.loads(Path(plan['reuse_readback']).read_text())
    if (rr['status'] != 'complete_full_pair_union_catalog_reuse_screen_not_authorization'
            or readback['status'] != 'passed_full_pair_union_reuse_candidate_sql_readback'
            or readback['producer_receipt_sha256'] != sha(root / 'receipt.json')):
        raise ValueError('Unverified full catalog/source screen')
    bindings.update(artifact_bindings(root, rr))
    for name in ['receipt.json', 'models.jsonl', 'model_pairs.tsv']:
        path = str(inventory / name)
        if rr['source_hashes'].get(path) != sha(path):
            raise ValueError('Source screen binds another reference ledger')
    candidates = {r['pair_key']: r for r in csv.DictReader((root / 'pair_reuse_candidates.tsv').open(), delimiter='\t')
                  if 'reference_expanded' in json.loads(r['new_sources'])}
    for row in rows:
        ends = endpoints(row)
        pair = hashlib.sha256(json.dumps(sorted(ends), separators=(',', ':')).encode()).hexdigest()
        if pair != row['pair_key'] or pair in full or pair not in candidates or sorted(ends) != endpoints(candidates[pair]) or ends[0] == ends[1]:
            raise ValueError('Invalid current reference pair/source identity')
        choice = choose_source(json.loads(candidates[pair]['matching_old_sources']))
        selected[pair] = choice
        full[pair] = row
    if set(full) != set(candidates):
        raise ValueError('Incomplete full reference source screen')
    expected = rr['per_new_source_counts'].get('reference_expanded:matching_catalog_sources_pending_input_and_result_checks', 0)
    if sum(bool(s) for s in selected.values()) != expected:
        raise ValueError('Full source-matched count differs')
    current = load_inputs(plan['current_inputs'], wanted, bindings, proofs=True)
    old = {label: load_source(spec, {p for p, owner in selected.items() if owner == label}, bindings)
           for label, spec in plan['sources'].items()}
    if set(s for s in selected.values() if s) - set(old):
        raise ValueError('Unconfigured selected source')
    return full, selected, current, curmodels, old, bindings


def inspect_checkpoint(source, pair, mask, old_order, bindings):
    key = pair, mask, old_order
    item = source['proofs'][key]
    path = source['root'] / item['path']
    if item['path'] != f'pairs/{pair[:2]}/{pair}-{mask}-{old_order}.json' or sha(path) != item['sha256']:
        raise ValueError('Old checkpoint path/hash differs')
    bindings[str(path)] = item['sha256']
    record = json.loads(path.read_text())
    ends = source['pairs'][pair][::1 if old_order == 0 else -1]
    rows = [source['inputs'][(*end, mask)] for end in ends]
    identities = [dict(model_id=e[0], version=e[1], status=r['status'], path=r.get('path'), sha256=r.get('sha256')) for e, r in zip(ends, rows)]
    ready = all(r['status'] == 'ready' for r in rows)
    command = [source['plan']['usalign'], *[r['path'] for r in rows], *source['plan']['options']] if ready else None
    if (record['pair_key'] != pair or record['mask'] != mask or record['order'] != old_order
            or record['plan_sha256'] != source['plan_sha256']
            or record['input_manifest_sha256'] != source['bundle'] or record['inputs'] != identities
            or record['command'] != command or record['status'] != item['status']):
        raise ValueError('Old directed native checkpoint lineage differs')
    numeric, geometry = source['numeric'].get(key), source['geometry'].get(key)
    if record['status'] == 'aligned':
        if (not ready or record['returncode'] != 0 or numeric is None or geometry is None
                or parse_output(record['stdout'], [r['sequence'] for r in rows]) != record['metrics']
                or numeric['rmsd_status'] != geometry['rmsd_status']
                or numeric['aligned_length'] != geometry['aligned_length']):
            raise ValueError('Old successful checkpoint/proof differs')
    else:
        if numeric is not None or geometry is not None:
            raise ValueError('Nonaligned native disposition has numeric geometry')
        if record['status'] == 'input_unavailable':
            if ready or record['input_statuses'] != [r['status'] for r in rows]:
                raise ValueError('Old unavailable inputs differ')
        elif record['status'] == 'native_error':
            if not ready or record['returncode'] == 0:
                raise ValueError('Invalid old native error')
        elif record['status'] == 'parse_error':
            if not ready or record['returncode'] != 0 or not record.get('error'):
                raise ValueError('Invalid old parse error')
        elif record['status'] == 'timeout':
            if not ready or record['timeout_seconds'] != source['plan']['per_pair_timeout_seconds']:
                raise ValueError('Invalid old timeout')
        else:
            raise ValueError('Unknown old native disposition')
    return path, item['sha256'], record, numeric, geometry


def numeric_exclusions(record, numeric, geometry):
    if record['status'] != 'aligned':
        return [record['status']]
    reasons = []
    if numeric['rmsd_status'] != 'within_printed_rounding':
        reasons.append('rmsd_discrepancy')
    if int(numeric['aligned_length']) < 3:
        reasons.append('fewer_than_three_pairs')
    if geometry['geometry_status'] != 'unique_at_numeric_tolerance':
        reasons.append('nonunique_rotation')
    return reasons


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    initial = {str(args.plan): sha(args.plan), **plan['pins']}
    def verify(bindings):
        for path, value in bindings.items():
            if sha(path) != value:
                raise ValueError('Changed pinned source: ' + path)
    verify(initial)
    full, selected, current, curmodels, old, bindings = load_sources(plan)
    bindings.update(initial)
    verify(bindings)
    target = json.loads(Path(plan['target_alignment_plan']).read_text())
    executable_sha = sha(target['usalign'])
    for source in old.values():
        previous = source['plan']
        if (previous['options'] != target['options'] or sha(previous['usalign']) != executable_sha
                or previous['pins'][previous['usalign']] != executable_sha):
            raise ValueError('Different native executable/settings require separate workload')
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    checked_inputs = {}
    checked_files = set()
    def bind_file(path, expected):
        if (path, expected) not in checked_files:
            if sha(path) != expected:
                raise ValueError('Changed actual source bytes: ' + path)
            checked_files.add((path, expected))
        if path in bindings and bindings[path] != expected:
            raise ValueError('Conflicting actual source binding')
        bindings[path] = expected
    @lru_cache(maxsize=None)
    def model_input_check(owner, key, mask):
        a, b = current[(*key, mask)], old[owner]['inputs'][(*key, mask)]
        model = curmodels[key]
        previous = old[owner]['models'][key]
        compatible = all(model[f] == previous[f] for f in ['sha256', 'sequence_sha256', 'length']) and projection(a) == projection(b)
        for entry in [model, previous]:
            bind_file(entry['path'], entry['sha256'])
        for row in [a, b]:
            if row['source_sha256'] != model['sha256']:
                raise ValueError('Input source hash differs from raw coordinate catalog')
            if row['status'] == 'ready':
                bind_file(row['path'], row['sha256'])
        checked_inputs[(owner, *key, mask)] = dict(source=owner, model_id=key[0], version=key[1], mask=mask,
                    compatible=compatible, current_projection_sha256=digest(projection(a)),
                    old_projection_sha256=digest(projection(b)), current_raw_path=model['path'],
                    old_raw_path=previous['path'], raw_sha256=model['sha256'],
                    current_pdb_path=a.get('path'), old_pdb_path=b.get('path'))
        return compatible
    totals, numerical, sources = Counter(), Counter(), Counter()
    result_path = out / 'reference_reuse_dispositions.jsonl'
    with result_path.open('w') as handle:
        for i, (pair, row) in enumerate(sorted(full.items()), 1):
            ends, owner = endpoints(row), selected[pair]
            for mask in MASKS:
                compatible = bool(owner) and all([model_input_check(owner, key, mask) for key in ends])
                for order in [0, 1]:
                    record = dict(pair_key=pair, model_a=ends[0][0], version_a=ends[0][1],
                                  model_b=ends[1][0], version_b=ends[1][1], mask=mask, order=order,
                                  selected_source=owner, numerical_usable=False,
                                  source_checkpoint=None, source_checkpoint_sha256=None, source_order=None,
                                  source_native_status=None, source_numeric=None, source_geometry=None,
                                  numerical_exclusion_reasons=[])
                    if not owner:
                        record['reuse_status'] = 'new_native_measurement_pending'
                    elif not compatible:
                        record['reuse_status'] = 'incompatible_inputs_require_new_measurement'
                    else:
                        previous = old[owner]
                        target_ends = ends[::1 if order == 0 else -1]
                        if previous['pairs'][pair] == target_ends:
                            source_order = 0
                        elif previous['pairs'][pair][::-1] == target_ends:
                            source_order = 1
                        else:
                            raise ValueError('Different source order/endpoints')
                        path, h, native, numeric, geometry = inspect_checkpoint(previous, pair, mask, source_order, bindings)
                        reasons = numeric_exclusions(native, numeric, geometry)
                        record.update(reuse_status='verified_identical_input_checkpoint_and_retained_disposition',
                                      source_checkpoint=str(path), source_checkpoint_sha256=h, source_order=source_order,
                                      source_native_status=native['status'], source_numeric=numeric, source_geometry=geometry,
                                      numerical_exclusion_reasons=reasons, numerical_usable=not reasons)
                        numerical['usable' if not reasons else ';'.join(reasons)] += 1
                    totals[record['reuse_status']] += 1
                    sources[owner or 'no_old_source'] += 1
                    handle.write(json.dumps(record, separators=(',', ':')) + '\n')
            if i % 1000 == 0:
                (out / 'state.json').write_text(json.dumps(dict(stage='qualifying_reference_reuse', pairs=i,
                                                              total_pairs=len(full), counts=dict(totals))) + '\n')
                print('Reference reuse pairs', i, '/', len(full), flush=True)
    checks = out / 'model_input_identity_checks.jsonl'
    with checks.open('w') as handle:
        for key, record in sorted(checked_inputs.items()):
            handle.write(json.dumps(record, separators=(',', ':')) + '\n')
    verify(bindings)
    receipt = dict(status='complete_full_reference_reuse_qualification_pending_independent_readback',
                   plan_sha256=sha(args.plan), full_reference_pairs=len(full), directed_dispositions=4 * len(full),
                   counts=dict(totals), numerical_counts=dict(numerical), selected_source_dispositions=dict(sources),
                   unique_model_mask_source_checks=len(checked_inputs), source_hashes=bindings,
                   artifacts={p.name: sha(p) for p in [result_path, checks]},
                   scientific_eligibility=False, scope='Complete full-reference pair/two-mask/two-order universe; '
                   'selected completed source is expanded primary first, otherwise old reference. Exact raw '
                   'coordinate, full sequence/length, mask/status/residue positions/sequence/PDB bytes, executable '
                   'settings and directed native checkpoints qualified; all old numeric/geometry flags retained. '
                   'Original checkpoints remain unchanged and referenced with hashes; unavailable/error/degenerate '
                   'dispositions are retained, never zero distances. Source matches alone do not qualify results; '
                   'incompatibilities require native measurement. Full independent reader and union with new '
                   'measurements, order/coverage/domain/PAE/orthology/phylogenetic inference remain required.')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: receipt[k] for k in ['status', 'full_reference_pairs', 'directed_dispositions', 'counts', 'numerical_counts']}), flush=True)


if __name__ == '__main__':
    main()
