#!/usr/bin/env python3
"""Check supplementary inventory binding, exact partition and native validation handoff."""
import argparse
import copy
import json
import subprocess
import sys
import tempfile
from pathlib import Path

from run_ortholog_pair_guide_comparison import sha
from validate_duplication_reference_coordinates import load_additional_models
from readback_duplication_coordinates import check_shard


def write(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--models', type=Path, required=True)
    args = parser.parse_args()
    with args.models.open() as handle:
        models = [json.loads(next(handle)) for _ in range(2)]
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        inventory, base = root / 'inventory', root / 'base'
        inventory.mkdir()
        base.mkdir()
        for i, model in enumerate(models):
            source = root / f'{i}.cif'
            source.write_bytes(Path(model['path']).read_bytes())
            model['path'] = str(source)
        (base / 'models.jsonl').write_text(json.dumps(models[0]) + '\n')
        write(base / 'receipt.json', dict(
            status='complete_reviewed_duplication_model_pair_queue', unique_models=1,
            artifacts={'models.jsonl': sha(base / 'models.jsonl')}))
        (inventory / 'models.jsonl').write_text(''.join(json.dumps(m)+'\n' for m in models))
        (inventory / 'additional_models.jsonl').write_text(json.dumps(models[1])+'\n')
        receipt = dict(status='complete_provisional_reference_comparison_inventory',
                       additional_models=1, all_reference_comparison_models=2,
                       source_pins={str(base/name): sha(base/name)
                                    for name in ['receipt.json', 'models.jsonl']},
                       artifacts={name: sha(inventory/name)
                                  for name in ['models.jsonl', 'additional_models.jsonl']})
        audit = dict(status='passed_full_reference_comparison_ledger_readback',
                     additional_models=1, unique_models=2)

        def bind():
            receipt['artifacts']['additional_models.jsonl'] = sha(inventory/'additional_models.jsonl')
            write(inventory/'receipt.json', receipt)
            audit['producer_receipt_sha256'] = sha(inventory/'receipt.json')
            write(root/'audit.json', audit)

        bind()
        plan = dict(inventory=str(inventory), base_queue=str(base),
                    inventory_readback=str(root/'audit.json'), output=str(root/'out'),
                    workers=1, models_per_shard=1, minimum_free_disk_gib=0, pins={})
        assert load_additional_models(plan) == [models[1]]
        # Even rehashed output with a valid-looking count cannot substitute a base model.
        (inventory/'additional_models.jsonl').write_text(json.dumps(models[0])+'\n')
        bind()
        try:
            load_additional_models(plan)
        except ValueError:
            pass
        else:
            raise AssertionError('Wrong partition accepted')
        (inventory/'additional_models.jsonl').write_text(json.dumps(models[1])+'\n')
        bind()
        bad = copy.deepcopy(audit)
        bad['producer_receipt_sha256'] = 'wrong'
        write(root/'audit.json', bad)
        try:
            load_additional_models(plan)
        except ValueError:
            pass
        else:
            raise AssertionError('Unbound audit accepted')
        bind()
        write(root/'plan.json', plan)
        command = [sys.executable, 'scripts/validate_duplication_reference_coordinates.py',
                   '--plan', str(root/'plan.json')]
        subprocess.run(command, check=True)
        result = json.loads((root/'out/receipt.json').read_text())
        assert result['models'] == 1 and result['counts'] == {'validated': 1}
        proofs = root/'proofs'
        proofs.mkdir()
        proof = check_shard((str(root/'out'), result['shards'][0], [models[1]],
                             str(proofs), 'fixture'))
        assert proof['validated_residues'] == models[1]['length']
        subprocess.run(command, check=True)
    print('Passed exact partition, changed audit, producer/readback handoff and checkpoint reuse.')


if __name__ == '__main__':
    main()
