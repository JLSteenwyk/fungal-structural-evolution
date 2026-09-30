"""Bind current audited inputs and partition the complete reference ledger by prior sources."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha
from validate_duplication_reference_coordinates_v2 import load_additional_models, model_map


def load_all_inputs(plan):
    inventory = Path(plan['inventory'])
    base = Path(plan['base_queue'])
    qr = json.loads((inventory/'receipt.json').read_text())
    # Recheck the exact model partition and the full comparison ledger binding.
    extras = load_additional_models(dict(inventory=str(inventory), base_queue=str(base),
                                        inventory_readback=plan['inventory_readback']))
    original = model_map(base/'models.jsonl')
    expected_sources = {'primary': original,
                        'additional': {(m['model_id'],m['version']):m for m in extras}}
    inputs, bindings = {}, {}
    for label in ['primary','additional']:
        spec = plan['input_sources'][label]
        folder = Path(spec['inputs'])
        receiptpath = folder/'receipt.json'
        receipt = json.loads(receiptpath.read_text())
        materialization = json.loads(Path(spec['input_plan']).read_text())
        if Path(materialization['output']).resolve() != folder.resolve():
            raise ValueError('Materialization directory differs')
        status = ('complete_duplication_alignment_input_materialization' if label == 'primary'
                  else 'complete_additional_reference_alignment_input_materialization')
        field, upstream = (('queue_receipt_sha256', base) if label == 'primary'
                           else ('inventory_receipt_sha256', inventory))
        if (receipt['status'] != status or receipt['plan_sha256'] != sha(spec['input_plan'])
                or receipt[field] != sha(upstream/'receipt.json')):
            raise ValueError('Incomplete or unbound materialization: '+label)
        manifest = folder/'inputs.jsonl'
        if sha(manifest) != receipt['artifacts']['inputs.jsonl']:
            raise ValueError('Changed input manifest')
        expected = {(m,v,mask) for m,v in expected_sources[label] for mask in ['full','plddt70']}
        seen = set()
        counts = Counter()
        with manifest.open() as handle:
            for line in handle:
                row = json.loads(line)
                key = row['model_id'],row['version'],row['mask']
                if key not in expected or key in seen or key in inputs:
                    raise ValueError('Duplicate or unexpected materialized model')
                if row['source_sha256'] != expected_sources[label][key[:2]]['sha256']:
                    raise ValueError('Materialized source hash differs')
                if row['status'] not in ['ready','too_few_retained_residues','source_rejected']:
                    raise ValueError('Invalid input disposition')
                counts[row['mask']+':'+row['status']] += 1
                seen.add(key)
                inputs[key] = row
        if (seen != expected or dict(counts) != receipt['counts']
                or receipt['models'] != len(expected_sources[label])
                or receipt['input_dispositions'] != len(expected)):
            raise ValueError('Incomplete materialization universe')
        bindings[str(receiptpath)] = sha(receiptpath)
        bindings[str(manifest)] = sha(manifest)
    pairsfile = inventory/'model_pairs.tsv'
    with pairsfile.open() as handle:
        all_pairs = list(csv.DictReader(handle, delimiter='\t'))
    base_receipt = json.loads((base/'receipt.json').read_text())
    if sha(base/'model_pairs.tsv') != base_receipt['artifacts']['model_pairs.tsv']:
        raise ValueError('Changed primary pair table')
    with (base/'model_pairs.tsv').open() as handle:
        base_keys = {row['pair_key'] for row in csv.DictReader(handle, delimiter='\t')}
    pairs, seen, reused = [], set(), 0
    for row in all_pairs:
        endpoints = [(row['model_a'],int(row['version_a'])),(row['model_b'],int(row['version_b']))]
        key = hashlib.sha256(json.dumps(sorted(endpoints),separators=(',',':')).encode()).hexdigest()
        if key != row['pair_key'] or key in seen or endpoints[0] == endpoints[1]:
            raise ValueError('Invalid reference pair identity')
        seen.add(key)
        for endpoint in endpoints:
            if any((*endpoint,mask) not in inputs for mask in ['full','plddt70']):
                raise ValueError('Missing reference pair input')
        if key in base_keys:
            if row['work_disposition'] != 'existing_duplicate_pair':
                raise ValueError('Primary pair reuse differs')
            reused += 1
        else:
            if row['work_disposition'] != 'additional_pair':
                raise ValueError('Additional pair disposition differs')
            pairs.append(row)
    if (len(pairs) != qr['additional_model_pairs'] or reused != qr['existing_duplicate_model_pairs']
            or len(all_pairs) != qr['unique_distinct_model_pairs']):
        raise ValueError('Reference pair scope differs')
    # The job checkpoint binds the full two-source manifest/receipt bundle.
    bundle_hash = hashlib.sha256(json.dumps(bindings,sort_keys=True).encode()).hexdigest()
    return inputs,pairs,bindings,bundle_hash


def load_handoff(plan):
    inputs, _, bindings, _ = load_all_inputs(plan)
    for label in ['primary', 'additional']:
        spec = plan['input_sources'][label]
        proof = json.loads(Path(spec['readback']).read_text())
        status = ('passed_full_duplication_alignment_input_readback' if label == 'primary'
                  else 'passed_full_reference_alignment_input_readback')
        rp = Path(spec['inputs']) / 'receipt.json'
        field = 'source_receipt_sha256' if label == 'primary' else 'producer_receipt_sha256'
        if proof['status'] != status or proof[field] != sha(rp):
            raise ValueError('Written inputs lack matching full readback: ' + label)
        bindings[spec['readback']] = sha(spec['readback'])
    root = Path(plan['reuse_candidates'])
    rp = root / 'receipt.json'
    receipt = json.loads(rp.read_text())
    proof = json.loads(Path(plan['reuse_readback']).read_text())
    path = root / 'pair_reuse_candidates.tsv'
    if (receipt['status'] != 'complete_full_pair_union_catalog_reuse_screen_not_authorization'
            or proof['status'] != 'passed_full_pair_union_reuse_candidate_sql_readback'
            or proof['producer_receipt_sha256'] != sha(rp)
            or sha(path) != receipt['artifacts'][path.name]):
        raise ValueError('Unverified full source partition')
    for name in ['receipt.json', 'models.jsonl', 'model_pairs.tsv']:
        p = str(Path(plan['inventory']) / name)
        if receipt['source_hashes'].get(p) != sha(p):
            raise ValueError('Partition binds another current reference source')
    candidates = {}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if 'reference_expanded' not in json.loads(row['new_sources']):
                continue
            if row['pair_key'] in candidates:
                raise ValueError('Duplicate partition pair')
            candidates[row['pair_key']] = row
    full = []
    pairs = []
    seen = set()
    counts = Counter()
    with (Path(plan['inventory']) / 'model_pairs.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = row['pair_key']
            if key not in candidates or key in seen:
                raise ValueError('Reference partition universe differs')
            seen.add(key)
            prior = candidates[key]
            ends = sorted([(row['model_a'], int(row['version_a'])), (row['model_b'], int(row['version_b']))])
            previous = [(prior['model_a'], int(prior['version_a'])), (prior['model_b'], int(prior['version_b']))]
            if ends != previous:
                raise ValueError('Partition endpoint/version differs')
            matching = json.loads(prior['matching_old_sources'])
            if matching:
                if prior['disposition'] != 'matching_catalog_sources_pending_input_and_result_checks':
                    raise ValueError('Conflicting pending reuse disposition')
                disposition = 'pending_actual_input_result_and_numeric_reuse_checks'
            else:
                if prior['disposition'] not in ['new_pair_requires_alignment', 'changed_catalog_sources_require_alignment']:
                    raise ValueError('Conflicting new measurement disposition')
                disposition = 'native_measurement_required'
                pairs.append(row)
            counts[prior['disposition']] += 1
            full.append({**row, **{f: prior[f] for f in ['new_sources', 'matching_old_sources', 'changed_old_sources']},
                         'catalog_disposition': prior['disposition'], 'measurement_disposition': disposition})
    if seen != set(candidates):
        raise ValueError('Incomplete reference work partition')
    expected_counts = {k.split(':', 1)[1]: v for k, v in receipt['per_new_source_counts'].items()
                       if k.startswith('reference_expanded:')}
    if dict(counts) != expected_counts:
        raise ValueError('Reference partition count differs')
    bindings.update({str(rp): sha(rp), str(path): sha(path), plan['reuse_readback']: sha(plan['reuse_readback'])})
    # Checkpoints bind their two materialized receipt/manifests. Full source/proof
    # bindings are additionally retained in the pinned plan; no old results imported.
    manifest_bindings = {str(Path(s['inputs']) / name): sha(Path(s['inputs']) / name)
                         for s in plan['input_sources'].values() for name in ['receipt.json', 'inputs.jsonl']}
    bundle = hashlib.sha256(json.dumps(manifest_bindings, sort_keys=True).encode()).hexdigest()
    return inputs, pairs, bindings, bundle, sorted(full, key=lambda r: r['pair_key'])


def write_work_partition(path, rows):
    with Path(path).open('w') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
