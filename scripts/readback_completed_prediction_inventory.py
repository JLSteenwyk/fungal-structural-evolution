#!/usr/bin/env python3
"""Compare every combined record to sources and inventory exact marker links."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):
            h.update(b)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan',type=Path,required=True)
    ap.add_argument('--global-links',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    plan=json.loads(a.plan.read_text())
    folder=Path(plan['output'])
    receipt=json.loads((folder/'receipt.json').read_text())
    pins={str(a.plan):sha(a.plan),str(a.global_links):sha(a.global_links),str(folder/'receipt.json'):sha(folder/'receipt.json'),**plan['pins']}
    for name,digest in receipt['artifacts'].items():
        pins[str(folder/name)]=digest
    def verify():
        for path,digest in pins.items():
            if sha(path)!=digest:
                raise ValueError('Changed input: '+path)
    verify()
    if receipt['status']!='complete_disjoint_prediction_inventory_union_with_batch_provenance' or receipt['plan_sha256']!=sha(a.plan):
        raise ValueError('Wrong union receipt')
    expected={}; cohorts={}; models=set(); batch_counts=Counter()
    for source in plan['cohorts']:
        count=0
        with (Path(source['inventory'])/'inventory.jsonl').open() as f:
            for line in f:
                row=json.loads(line); sid=row['record_id']; model=row['models'][0]
                if sid in expected or model['model_id'] in models:
                    raise ValueError('Overlapping source universe')
                expected[sid]=row;cohorts[sid]=source;models.add(model['model_id']);count+=1
                batch_counts[model['prediction_config_sha256']]+=1
        if count!=source['models']:
            raise ValueError('Source scope differs')
    observed=set()
    with (folder/'inventory.jsonl').open() as f:
        for line in f:
            row=json.loads(line);sid=row['record_id']
            if sid in observed or row!=expected.get(sid):
                raise ValueError('Changed, repeated or additional source record')
            observed.add(sid)
    if observed!=set(expected) or len(observed)!=receipt['models'] or len(observed)!=plan['expected_models']:
        raise ValueError('Missing source models')
    seen=set()
    with (folder/'model_cohorts.tsv').open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            sid=row['sequence_id']
            if sid in seen or sid not in expected:
                raise ValueError('Repeated or additional batch binding')
            model=expected[sid]['models'][0];source=cohorts[sid]
            correct=dict(cohort=source['name'],sequence_id=sid,model_id=model['model_id'],
                         prediction_config_sha256=model['prediction_config_sha256'],length=str(model['length']),
                         source_inventory=str(Path(source['inventory'])/'inventory.jsonl'))
            if row!=correct:
                raise ValueError('Changed batch binding')
            seen.add(sid)
    if seen!=observed:
        raise ValueError('Missing batch bindings')
    configs=json.loads((folder/'prediction_configurations.json').read_text())
    correct={sha(p):dict(path=p,configuration=json.loads(Path(p).read_text())) for p in plan['configurations']}
    if configs['configurations']!=correct or configs['models_by_configuration']!=dict(batch_counts):
        raise ValueError('Changed configuration catalog')
    taxa=Counter();markers=Counter();cohort_links=Counter();matched=[];universe=set();total=0
    with a.global_links.open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            total+=1
            sid='S'+row['sequence_sha256'];universe.add(sid)
            if sid in expected:
                matched.append(row);taxa[row['taxon_id']]+=1;markers[row['marker']]+=1
                cohort_links[cohorts[sid]['name']]+=1
    if not observed<=universe:
        raise ValueError('Models outside global marker inventory')
    a.output.mkdir()
    for name,data,key in [('taxon_model_links.tsv',taxa,'taxon_id'),('marker_model_links.tsv',markers,'marker'),('cohort_model_links.tsv',cohort_links,'cohort')]:
        with (a.output/name).open('w') as f:
            writer=csv.writer(f,delimiter='\t',lineterminator='\n');writer.writerow([key,'exact_marker_protein_links'])
            writer.writerows(sorted(data.items()))
    with (a.output/'expected_global_links.tsv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(matched[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(matched)
    verify()
    result=dict(status='passed_complete_disjoint_prediction_inventory_readback',models=len(observed),configurations=len(correct),
                global_marker_sequences=len(universe),global_marker_links=total,matched_marker_protein_links=len(matched),
                represented_taxa=len(taxa),represented_markers=len(markers),models_by_configuration=dict(batch_counts),
                producer_receipt_sha256=sha(folder/'receipt.json'),input_hashes=pins,script_sha256=sha(__file__),
                artifacts={p.name:sha(p) for p in a.output.iterdir()},
                scope='Every combined model record and batch binding matched exactly to source; full configurations preserved. Coverage is exact sequence availability among marker records, not confidence-qualified structural coverage, all-proteome coverage or validated biological accuracy.')
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['input_hashes','artifacts','models_by_configuration']},indent=2))


if __name__=='__main__':
    main()
