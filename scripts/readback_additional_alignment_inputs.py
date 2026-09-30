#!/usr/bin/env python3
"""Check every supplemental reference/background PDB and its complete native-model partition."""
import argparse
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

from readback_duplication_alignment_inputs import check_entry
from screen_duplication_alignment_reuse import sha


ROLES = {
    'background': ('complete_background_alignment_input_materialization',
                   'complete_background_coordinate_validation_with_dispositions',
                   'passed_background_exported_ca_readback', 'source_inventory_receipt_sha256'),
    'reference': ('complete_additional_reference_alignment_input_materialization',
                  'complete_duplication_coordinate_validation_with_dispositions',
                  'passed_duplication_exported_ca_readback', 'inventory_receipt_sha256'),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    bindings = {str(args.plan): sha(args.plan), **plan['pins']}

    def bind(path, digest=None):
        path = str(path)
        observed = sha(path)
        if digest is not None and observed != digest:
            raise ValueError('Changed artifact: ' + path)
        if path in bindings and bindings[path] != observed:
            raise ValueError('Changed frozen source: ' + path)
        bindings[path] = observed
        return observed

    def verify():
        for path, digest in bindings.items():
            if sha(path) != digest:
                raise ValueError('Changed source: ' + path)

    verify()
    status, coordinate_status, readback_status, inventory_field = ROLES[plan['role']]
    pp = Path(plan['producer_plan'])
    producer = json.loads(pp.read_text())
    for name, digest in producer['pins'].items():
        bind(name, digest)
    root = Path(producer['output'])
    rp = root / 'receipt.json'
    receipt = json.loads(rp.read_text())
    rh = bind(rp)
    if receipt['status'] != status or receipt['plan_sha256'] != bind(pp):
        raise ValueError('Incomplete materialization')
    source_plan = Path(producer['coordinate_plan'])
    sp = json.loads(source_plan.read_text())
    source = Path(producer['coordinates'])
    audit = Path(producer['readback'])
    if Path(sp['output']).resolve() != source.resolve():
        raise ValueError('Coordinate source differs')
    crp, arp = source / 'receipt.json', audit / 'receipt.json'
    cr, ar = json.loads(crp.read_text()), json.loads(arp.read_text())
    if (cr['status'] != coordinate_status or cr['plan_sha256'] != bind(source_plan)
            or ar['status'] != readback_status or ar['plan_sha256'] != bind(producer['readback_plan'])
            or receipt['coordinate_receipt_sha256'] != bind(crp)
            or ar['producer_receipt_sha256'] != sha(crp)
            or receipt['readback_receipt_sha256'] != bind(arp)):
        raise ValueError('Coordinate/readback lineage differs')
    ap = json.loads(Path(producer['readback_plan']).read_text())
    inventory = Path(sp['inventory'])
    modelpath = inventory / 'additional_models.jsonl'
    if (Path(ap['source']).resolve() != source.resolve()
            or Path(ap['source_plan']).resolve() != source_plan.resolve()
            or Path(ap['models']).resolve() != modelpath.resolve()):
        raise ValueError('Independent coordinate checker used another source')
    irp = inventory / 'receipt.json'
    ir = json.loads(irp.read_text())
    if cr[inventory_field] != bind(irp) or receipt['inventory_receipt_sha256'] != sha(irp):
        raise ValueError('Native inventory lineage differs')
    bind(modelpath, ir['artifacts'][modelpath.name])
    with modelpath.open() as handle:
        models = [json.loads(line) for line in handle]
    expected = {(m['model_id'], m['version']): m for m in models}
    if len(expected) != len(models) or len(expected) != plan['models'] or ir['additional_models'] != len(models):
        raise ValueError('Additional model partition differs')
    if any(record['models'] != len(models) for record in [receipt, cr, ar]):
        raise ValueError('Incomplete model universe')
    manifest = root / 'inputs.jsonl'
    bind(manifest, receipt['artifacts'][manifest.name])
    seen, counts, bytes_checked = set(), Counter(), 0
    n = sp['models_per_shard']
    if len(cr['shards']) != (len(models) + n - 1) // n:
        raise ValueError('Incomplete coordinate shard set')
    with manifest.open() as rows:
        for i, shard in enumerate(cr['shards']):
            subset = models[i*n:(i+1)*n]
            mh = hashlib.sha256(json.dumps(subset, sort_keys=True).encode()).hexdigest()
            if shard['models_sha256'] != mh or shard['models'] != len(subset) or shard['output'] != f'shard-{i:05d}.jsonl.gz':
                raise ValueError('Native coordinate shard partition differs')
            sourcepath = source / shard['output']
            bind(sourcepath, shard['output_sha256'])
            proofpath = audit / (sourcepath.name + '.readback.json')
            bind(proofpath, ar['proofs'][proofpath.name])
            proof = json.loads(proofpath.read_text())
            if proof['source_sha256'] != sha(sourcepath) or proof['models'] != len(subset):
                raise ValueError('Incomplete coordinate shard proof')
            local = Counter()
            with gzip.open(sourcepath, 'rt') as handle:
                for model in subset:
                    line = handle.readline()
                    if not line:
                        raise ValueError('Missing native coordinate record')
                    record = json.loads(line)
                    key = record['model_id'], record['version']
                    if key in seen or key != (model['model_id'], model['version']):
                        raise ValueError('Wrong/repeated model record')
                    for field, sourcefield in [('source_sha256', 'sha256'), ('source_path', 'path'), ('sequence_sha256', 'sequence_sha256')]:
                        if record[field] != model[sourcefield]:
                            raise ValueError('Native source descriptor differs')
                    seen.add(key)
                    local[record['status']] += 1
                    for mask in ['full', 'plddt70']:
                        row = json.loads(next(rows))
                        if row['coordinate_shard'] != sourcepath.name:
                            raise ValueError('Materialization shard differs')
                        blob = None
                        if row['status'] == 'ready':
                            name = hashlib.sha256(json.dumps(key, separators=(',', ':')).encode()).hexdigest()
                            path = root / 'pdb' / name[:2] / (name + '-' + mask + '.pdb')
                            if row['path'] != str(path):
                                raise ValueError('Materialized coordinate path differs')
                            blob = path.read_bytes()
                        bytes_checked += check_entry(row, record, mask, blob)
                        counts[mask + ':' + row['status']] += 1
                if handle.readline():
                    raise ValueError('Extra native coordinate record')
            if dict(local) != shard['counts'] or dict(local) != proof['counts']:
                raise ValueError('Coordinate shard dispositions differ')
            print('Checked supplemental alignment-input shards', i+1, '/', len(cr['shards']), flush=True)
        if rows.read():
            raise ValueError('Extra input disposition')
    if (seen != set(expected) or sum(counts.values()) != 2*len(models)
            or receipt['input_dispositions'] != 2*len(models) or dict(counts) != receipt['counts']
            or bytes_checked != receipt['pdb_bytes']):
        raise ValueError('Incomplete full input partition or summaries')
    verify()
    result = dict(status='passed_full_' + plan['role'] + '_alignment_input_readback',
                  models=len(models), input_dispositions=2*len(models), counts=dict(counts),
                  pdb_bytes=bytes_checked, producer_receipt_sha256=rh,
                  plan_sha256=sha(args.plan), source_hashes=bindings,
                  scientific_eligibility=False,
                  scope='Every additional native model and both confidence masks checked against full coordinate/shard/readback lineage; every ready PDB residue identity, original position, coordinate/confidence rounding and bytes hash independently verified. Short/rejected inputs retained. Production masks remain full and pLDDT70; no input reuse, pair alignment, rejection-cause adjudication, PAE or biological qualification.')
    with Path(plan['output']).open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: result[k] for k in ['status', 'models', 'input_dispositions', 'pdb_bytes']}), flush=True)


if __name__ == '__main__':
    main()
