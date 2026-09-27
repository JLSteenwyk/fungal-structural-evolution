"""Bind primary and supplementary inputs to the additional reference-pair workload."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha
from validate_duplication_reference_coordinates import load_additional_models, model_map


def load_handoff(plan):
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
                inputs[key] = {k:row[k] for k in ['status','path','sha256','sequence'] if k in row}
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
