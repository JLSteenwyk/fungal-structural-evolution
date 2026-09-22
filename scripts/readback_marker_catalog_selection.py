#!/usr/bin/env python3
"""Reconstruct source-specific marker catalog selection from frozen JSONL."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from collections import defaultdict


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8388608),b''):h.update(block)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--catalog',type=Path,required=True)
    p.add_argument('--markers',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    receipt_path=a.catalog/'receipt.json';receipt=json.loads(receipt_path.read_text())
    pins={str(receipt_path):sha(receipt_path),str(a.markers):sha(a.markers)}
    if receipt['status']!='complete_source_specific_model_catalog' or pins[str(a.markers)]!=receipt['marker_mapping_sha256']:
        raise ValueError('Invalid catalog or marker frame')
    for name,digest in receipt['artifacts'].items():
        if sha(a.catalog/name)!=digest:raise ValueError('Changed catalog artifact')
        pins[str(a.catalog/name)]=digest
    markers=list(csv.DictReader(a.markers.open(),delimiter='\t'))
    sequences={r['sequence_sha256'] for r in markers}
    first_seen={};latest_relevant={};digest=hashlib.sha256();n=0
    # Preserve accession insertion order even when its first record was irrelevant.
    with Path(receipt['inventory_path']).open('rb') as f:
        for raw in f:
            if not raw.endswith(b'\n'):raise ValueError('Inventory is not a complete-line snapshot')
            digest.update(raw);row=json.loads(raw);accession=row['uniprot_accession']
            first_seen.setdefault(accession,n);n+=1
            relevant=[m for m in row.get('models',[]) if m['sequence_sha256'] in sequences] if row['status']=='verified' else []
            if relevant:latest_relevant[accession]=relevant
            else:latest_relevant.pop(accession,None)
    if digest.hexdigest()!=receipt['inventory_sha256']:raise ValueError('Frozen inventory mismatch')
    candidates=defaultdict(list);policy=receipt['source_policy']
    for accession in sorted(latest_relevant,key=first_seen.get):
        for model in latest_relevant[accession]:
            if model.get('provider')==policy['provider'] and model.get('tool')==policy['tool']:
                candidates[model['sequence_sha256']].append(dict(model,uniprot_accession=accession))
    expected_models={};expected_links=[]
    for marker in markers:
        options=candidates[marker['sequence_sha256']]
        if not options:continue
        best=sorted(options,key=lambda m:(-m['mean_ca_plddt'],-m['version']),reverse=False)
        # Model ID breaks ties in descending lexical order; preserve accession
        # insertion order among otherwise identical candidates.
        score=max((m['mean_ca_plddt'],m['version'],m['model_id']) for m in best)
        model=next(m for m in options if (m['mean_ca_plddt'],m['version'],m['model_id'])==score)
        expected_models.setdefault(model['path'],model)
        expected_links.append(dict(marker=marker['marker'],taxon_id=marker['taxon_id'],protein_id=marker['protein_id'],
            sequence_sha256=marker['sequence_sha256'],model_id=model['model_id'],version=str(model['version']),model_path=model['path']))
    models=json.loads((a.catalog/'model_provenance.json').read_text())
    actual={m['path']:m for m in models}
    links=list(csv.DictReader((a.catalog/'marker_model_links.tsv').open(),delimiter='\t'))
    if len(actual)!=len(models) or actual!=expected_models or links!=expected_links:
        raise ValueError('Selected records or marker links differ from inventory reconstruction')
    if len(models)!=receipt['distinct_models'] or len(links)!=receipt['marker_proteins_linked']:
        raise ValueError('Catalog totals differ')
    for path,digest in pins.items():
        if sha(path)!=digest:raise ValueError('Source changed during readback')
    result=dict(status='passed_full_frozen_inventory_marker_catalog_selection_readback',
        inventory_records=n,models=len(models),marker_links=len(links),
        inventory_sha256=receipt['inventory_sha256'],catalog_receipt_sha256=pins[str(receipt_path)],
        script_sha256=sha(__file__),input_hashes=pins,
        scope='Full frozen log replay, latest-accession semantics, source policy, tie handling and exact selected model/link records checked without importing catalog producer functions. Does not repeat coordinate parsing or qualify structural confidence.')
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
