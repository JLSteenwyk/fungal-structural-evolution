#!/usr/bin/env python3
"""Prepare a validated candidate model inventory for existing PAE retrieval."""
import argparse
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--screen', type=Path, required=True)
    p.add_argument('--coordinate-readback', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    source = a.screen / 'receipt.json'
    candidates = a.screen / 'candidate_models.jsonl'
    receipt = json.loads(source.read_text())
    checked = json.loads(a.coordinate_readback.read_text())
    if checked['status'] != 'passed_coordinate_sequence_and_confidence_readback':
        raise ValueError('Coordinate validation incomplete')
    if checked['screen_receipt_sha256'] != sha(source):
        raise ValueError('Different screening receipt')
    if sha(candidates) != checked['candidate_inventory_sha256'] or sha(candidates) != receipt['artifacts']['candidate_models.jsonl']:
        raise ValueError('Changed candidates')
    models = [json.loads(line)['model'] for line in candidates.read_text().splitlines()]
    indexed = {(r['model_id'], r['version']): r for r in checked['results']}
    if len(indexed) != len(models) or len(models) != checked['models']:
        raise ValueError('Validation scope differs')
    seen = set()
    for model in models:
        key = model['model_id'], model['version']
        if key in seen:
            raise ValueError('Duplicate model')
        seen.add(key)
        row = indexed[key]
        for field in ['path', 'sequence_sha256', 'length']:
            if row[field] != model[field]:
                raise ValueError('Changed validated model')
        if row['coordinate_sha256'] != model['sha256'] or sha(Path(model['path'])) != model['sha256']:
            raise ValueError('Changed coordinates')
    a.output.mkdir(parents=True, exist_ok=False)
    target = a.output / 'model_provenance.json'
    target.write_text(json.dumps(models, indent=2) + '\n')
    result = dict(status='complete_recovered_candidate_provenance_only', models=len(models),
                  screen_receipt_sha256=sha(source), coordinate_readback_sha256=sha(a.coordinate_readback),
                  script_sha256=sha(Path(__file__)), artifacts={target.name: sha(target)},
                  scope='Validated candidate inventory for PAE retrieval; not a marker mapping or a confidence-qualified analysis snapshot.')
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
