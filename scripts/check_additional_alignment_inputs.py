#!/usr/bin/env python3
"""Exercise both supplemental materializers and full input checker with rehashed corruptions."""
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

from screen_duplication_alignment_reuse import sha


def save(path, record):
    Path(path).write_text(json.dumps(record) + '\n')


def main():
    for role in ['background', 'reference']:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            coords, audit, inventory, base = [root / n for n in ['coords', 'audit', 'inventory', 'base']]
            for folder in [coords, audit, inventory, base]:
                folder.mkdir()
            records = []
            for version, confidence in [(6, [70, 20, 90, 80, 60]), (10, [0]*5)]:
                records.append(dict(model_id='M', version=version, status='validated', sequence='ACDEF',
                                    length=5, ca_xyz=[[i, i/2, 0] for i in range(5)], ca_plddt=confidence,
                                    source_sha256='fixture-' + str(version), source_path='fixture.cif',
                                    sequence_sha256=hashlib.sha256(b'ACDEF').hexdigest()))
            records.append(dict(model_id='Z', version=1, status='rejected_content', reason='fixture rejection',
                                source_sha256='fixture-Z', source_path='fixture-Z.cif', sequence_sha256='fixture'))
            models = [dict(model_id=r['model_id'], version=r['version'], path=r['source_path'],
                           sha256=r['source_sha256'], sequence_sha256=r['sequence_sha256']) for r in records]
            blob = ''.join(json.dumps(row) + '\n' for row in models)
            (inventory / 'additional_models.jsonl').write_text(blob)
            (inventory / 'models.jsonl').write_text(blob)
            (base / 'models.jsonl').write_text('')
            save(base / 'receipt.json', dict(status='complete_reviewed_duplication_model_pair_queue', unique_models=0,
                                             artifacts={'models.jsonl': sha(base / 'models.jsonl')}))
            ir = dict(status=('complete_background_measurement_inventory_pending_readback' if role == 'background'
                              else 'complete_provisional_reference_comparison_inventory'), additional_models=3,
                      all_reference_comparison_models=3,
                      source_pins={str(base/name): sha(base/name) for name in ['receipt.json', 'models.jsonl']},
                      artifacts={name: sha(inventory/name) for name in ['models.jsonl', 'additional_models.jsonl']})
            save(inventory / 'receipt.json', ir)
            ia = root / 'inventory-proof.json'
            save(ia, dict(status=('passed_full_background_measurement_inventory_readback' if role == 'background'
                                  else 'passed_full_reference_comparison_ledger_and_native_model_readback'),
                          additional_models=3, unique_models=3, producer_receipt_sha256=sha(inventory / 'receipt.json')))
            cp, ap = root / 'coordinate-plan.json', root / 'audit-plan.json'
            save(cp, dict(inventory=str(inventory), base_queue=str(base), readback=str(ia), inventory_readback=str(ia),
                          output=str(coords), models_per_shard=3))
            save(ap, dict(source_plan=str(cp), source=str(coords), models=str(inventory / 'additional_models.jsonl')))
            shard = coords / 'shard-00000.jsonl.gz'
            with gzip.open(shard, 'wt') as handle:
                for row in records:
                    handle.write(json.dumps(row) + '\n')
            counts = {'validated': 2, 'rejected_content': 1}
            sr = dict(status=('complete_background_coordinate_validation_with_dispositions' if role == 'background'
                              else 'complete_duplication_coordinate_validation_with_dispositions'),
                      plan_sha256=sha(cp), models=3,
                      shards=[dict(output=shard.name, output_sha256=sha(shard), models=3, counts=counts,
                                   models_sha256=hashlib.sha256(json.dumps(models, sort_keys=True).encode()).hexdigest())])
            sr['source_inventory_receipt_sha256' if role == 'background' else 'inventory_receipt_sha256'] = sha(inventory/'receipt.json')
            save(coords / 'receipt.json', sr)
            sp = audit / (shard.name + '.readback.json')
            save(sp, dict(source_sha256=sha(shard), models=3, counts=counts))
            save(audit / 'receipt.json', dict(status=('passed_background_exported_ca_readback' if role == 'background'
                                                      else 'passed_duplication_exported_ca_readback'),
                                              plan_sha256=sha(ap), models=3, producer_receipt_sha256=sha(coords/'receipt.json'),
                                              proofs={sp.name: sha(sp)}))
            output = root / 'output'
            pp = root / 'producer-plan.json'
            save(pp, dict(producer=dict(pid=99999999, created=0, cmdline=[]),
                          coordinates=str(coords), coordinate_plan=str(cp), readback=str(audit), readback_plan=str(ap),
                          output=str(output), minimum_free_disk_gib=0,
                          pins={str(cp): sha(cp), str(ap): sha(ap)}))
            script = 'scripts/materialize_background_inputs.py' if role == 'background' else 'scripts/materialize_duplication_reference_inputs_v2.py'
            proc = subprocess.run([sys.executable, script, '--plan', str(pp)], capture_output=True, text=True)
            assert proc.returncode == 0, proc.stderr
            plan = root / 'readback-plan.json'
            proof = root / 'proof.json'
            save(plan, dict(role=role, models=3, producer_plan=str(pp), output=str(proof), pins={str(pp): sha(pp)}))
            command = [sys.executable, 'scripts/readback_additional_alignment_inputs.py', '--plan', str(plan)]
            proc = subprocess.run(command, capture_output=True, text=True)
            assert proc.returncode == 0, proc.stderr
            readback = json.loads(proof.read_text())
            assert readback['models'] == 3 and readback['input_dispositions'] == 6
            manifest, receipt = output / 'inputs.jsonl', output / 'receipt.json'
            pristine = manifest.read_bytes(), receipt.read_bytes()
            rows = [json.loads(line) for line in manifest.read_text().splitlines()]
            masked = next(row for row in rows if row['version'] == 6 and row['mask'] == 'plddt70')
            assert masked['original_positions'] == [1, 3, 4] and masked['sequence'] == 'ADE'

            def reject(edited):
                proof.unlink()
                manifest.write_text(''.join(json.dumps(row) + '\n' for row in edited))
                r = json.loads(receipt.read_text())
                r['artifacts'][manifest.name] = sha(manifest)
                save(receipt, r)
                proc = subprocess.run(command, capture_output=True, text=True)
                assert proc.returncode != 0 and not proof.exists()
                manifest.write_bytes(pristine[0])
                receipt.write_bytes(pristine[1])
                save(proof, {'fixture_reset': True})

            bad = json.loads(json.dumps(rows))
            next(r for r in bad if r['version'] == 6 and r['mask'] == 'plddt70')['original_positions'] = [1, 2, 4]
            reject(bad)
            bad = json.loads(json.dumps(rows));bad[0]['version'] = 10
            reject(bad)
            reject(rows[:-1])
            # Rehash both PDB and manifest: a wrong coordinate must still fail native reconstruction.
            bad = json.loads(json.dumps(rows))
            ready = next(r for r in bad if r['status'] == 'ready')
            pdb = Path(ready['path']);data = pdb.read_text()
            pdb.write_text(data[:30] + '  10.000' + data[38:])
            ready['sha256'] = sha(pdb)
            reject(bad)
    print('PASS: both full supplemental handoffs, numeric versions 6/10, six mask/disposition outcomes and sparse positions; reject rehashed wrong residues, versions, missing rows and PDB coordinates.')


if __name__ == '__main__':
    main()
