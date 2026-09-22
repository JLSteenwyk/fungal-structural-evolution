#!/usr/bin/env python3
"""Combine disjoint converted ESMFold inventories without erasing batch settings."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import shutil
import time


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8388608),b''):
            h.update(block)
    return h.hexdigest()


def identity(row,sequences,models):
    if row['status']!='verified' or len(row['models'])!=1:
        raise ValueError('Expected one verified model per sequence')
    model=row['models'][0]
    sid=row['record_id']
    if sid!='S'+model['sequence_sha256'] or model['provider']!='local' or model['tool']!='ESMFold v1':
        raise ValueError('Unexpected model identity or source')
    if sid in sequences or model['model_id'] in models:
        raise ValueError('Overlapping sequence/model requires explicit reconciliation')
    sequences.add(sid)
    models.add(model['model_id'])
    return model


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan',required=True,type=Path)
    a=ap.parse_args()
    plan=json.loads(a.plan.read_text())
    plan_sha=sha(a.plan)
    def verify():
        if sha(a.plan)!=plan_sha:
            raise ValueError('Plan changed')
        for p,digest in plan['pins'].items():
            if sha(p)!=digest:
                raise ValueError('Changed input: '+p)
    verify()
    out=Path(plan['output'])
    if out.exists():
        raise FileExistsError(out)
    if shutil.disk_usage(out.parent).free<plan['resources']['minimum_free_disk_gib']*2**30:
        raise ValueError('Insufficient disk')
    out.mkdir()
    start=time.time()
    configs={}
    for path in plan['configurations']:
        digest=sha(path)
        if digest in configs:
            raise ValueError('Duplicate configuration')
        configs[digest]=dict(path=path,configuration=json.loads(Path(path).read_text()))
    records=[]; bindings=[]; sources=[]; sequences=set(); models=set(); batches=Counter()
    for cohort in plan['cohorts']:
        folder=Path(cohort['inventory'])
        receipt=json.loads((folder/'receipt.json').read_text())
        if receipt['status']!=cohort['expected_status'] or receipt['models']!=cohort['models']:
            raise ValueError('Unexpected source conversion scope/status')
        for name,digest in receipt['artifacts'].items():
            if sha(folder/name)!=digest:
                raise ValueError('Source artifact changed')
        count=0
        with (folder/'inventory.jsonl').open() as f:
            for line in f:
                row=json.loads(line)
                model=identity(row,sequences,models)
                digest=model['prediction_config_sha256']
                if digest not in configs:
                    raise ValueError('Missing prediction configuration')
                prediction=Path(model['prediction_receipt_path'])
                if (sha(prediction)!=model['prediction_receipt_sha256']
                        or sha(prediction.parent/'config.json')!=digest
                        or sha(model['path'])!=model['sha256']):
                    raise ValueError('Changed prediction provenance or coordinates')
                records.append(row)
                bindings.append(dict(cohort=cohort['name'],sequence_id=row['record_id'],model_id=model['model_id'],
                                     prediction_config_sha256=digest,length=model['length'],
                                     source_inventory=str(folder/'inventory.jsonl')))
                batches[digest]+=1
                count+=1
                if count%1000==0:
                    (out/'state.json').write_text(json.dumps(dict(status='verifying_and_collecting',cohort=cohort['name'],cohort_models=count,total_models=len(records)))+'\n')
        if count!=cohort['models']:
            raise ValueError('Inventory model count differs')
        sources.append(dict(cohort=cohort['name'],inventory=str(folder),models=count,
                            receipt_sha256=sha(folder/'receipt.json'),inventory_sha256=sha(folder/'inventory.jsonl')))
    if len(records)!=plan['expected_models'] or set(batches)!=set(configs):
        raise ValueError('Wrong combined scope or unused configuration')
    with (out/'inventory.jsonl').open('w') as f:
        for row in sorted(records,key=lambda x:x['record_id']):
            f.write(json.dumps(row,sort_keys=True)+'\n')
    with (out/'model_cohorts.tsv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(bindings[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(sorted(bindings,key=lambda r:r['sequence_id']))
    config_keys=set().union(*(set(c['configuration']) for c in configs.values()))
    varying=[key for key in sorted(config_keys) if len({json.dumps(c['configuration'].get(key),sort_keys=True) for c in configs.values()})>1]
    (out/'prediction_configurations.json').write_text(json.dumps(dict(configurations=configs,models_by_configuration=dict(batches),varying_fields=varying),indent=2)+'\n')
    verify()
    result=dict(status='complete_disjoint_prediction_inventory_union_with_batch_provenance',models=len(records),
                source_cohorts=sources,configurations=len(configs),models_by_configuration=dict(batches),varying_configuration_fields=varying,
                plan_sha256=plan_sha,script_sha256=sha(__file__),elapsed_seconds=time.time()-start,
                artifacts={p.name:sha(p) for p in [out/'inventory.jsonl',out/'model_cohorts.tsv',out/'prediction_configurations.json']},
                scope='Unchanged source model records, disjoint exact sequences and model IDs, all coordinate and prediction receipt hashes checked. Configurations remain explicit. PAE files remain bound by original hashes; no new PAE or coordinate-accuracy audit. Mapping, feature qualification and evolutionary analyses remain separate.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    (out/'state.json').write_text(json.dumps(dict(status=result['status'],models=len(records)))+'\n')


if __name__=='__main__':
    main()
