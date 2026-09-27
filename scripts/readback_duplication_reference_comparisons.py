#!/usr/bin/env python3
"""Independently enumerate the complete reference ledger and model-work partitions."""
import argparse,csv,hashlib,json
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha

p=argparse.ArgumentParser(description=__doc__)
for n in ['inventory','references','base-queue','output']:p.add_argument('--'+n,type=Path,required=True)
a=p.parse_args()
if a.output.exists():raise FileExistsError(a.output)
r=json.loads((a.inventory/'receipt.json').read_text());refs=json.loads((a.references/'receipt.json').read_text());base=json.loads((a.base_queue/'receipt.json').read_text())
for name,h in r['artifacts'].items():
    if sha(a.inventory/name)!=h:raise ValueError('Changed output artifact')
for path,h in r['source_pins'].items():
    if sha(path)!=h:raise ValueError('Changed source pin')
expected=set()
for guide in ['profile','mafft']:
    path=a.references/(guide+'_sister_references.tsv')
    if sha(path)!=refs['artifacts'][path.name] or r['source_pins'][str(path)]!=sha(path):raise ValueError('Reference table mismatch')
    with path.open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            if row['status']!='provisional_reference_available':continue
            for gene in json.loads(row['nearest_reference_genes']):
                for side in ['a','b']:expected.add((guide,row['family'],row['gene_node'],row['gene_a'],row['gene_b'],gene,side,row['gene_'+side],str(int(gene==row['chosen_reference_gene']))))
path=a.base_queue/'models.jsonl'
if sha(path)!=base['artifacts'][path.name] or r['source_pins'][str(path)]!=sha(path):raise ValueError('Base models mismatch')
with path.open() as f:base_models={(m['model_id'],str(m['version'])) for line in f for m in [json.loads(line)]}
path=a.base_queue/'model_pairs.tsv'
if sha(path)!=base['artifacts'][path.name]:raise ValueError('Base pairs changed')
with path.open() as f:base_pairs={tuple(sorted([(x['model_a'],x['version_a']),(x['model_b'],x['version_b'])])) for x in csv.DictReader(f,delimiter='\t')}
observed=set();models=set();pairs=set();counts={}
with (a.inventory/'event_reference_comparisons.tsv').open() as f:
    for row in csv.DictReader(f,delimiter='\t'):
        key=tuple(row[k] for k in ['guide','family','gene_node','gene_a','gene_b','reference_gene','focal_side','focal_gene','lexical_representative'])
        if key in observed:raise ValueError('Repeated event reference link')
        observed.add(key);left=(row['focal_model'],row['focal_version']);right=(row['reference_model'],row['reference_version']);models.update([left,right]);pair=tuple(sorted([left,right]))
        digest=hashlib.sha256(json.dumps([[m,int(v)] for m,v in pair],separators=(',',':')).encode()).hexdigest()
        if digest!=row['pair_key']:raise ValueError('Pair digest mismatch')
        status='identical_model' if left==right else ('existing_duplicate_pair' if pair in base_pairs else 'additional_pair')
        if row['work_disposition']!=status:raise ValueError('Work partition differs')
        counts[status]=counts.get(status,0)+1
        if left!=right:pairs.add(pair)
if observed!=expected:raise ValueError('Full tied-reference ledger differs')
with (a.inventory/'model_pairs.tsv').open() as f:
    rows=list(csv.DictReader(f,delimiter='\t'));written={tuple(sorted([(x['model_a'],x['version_a']),(x['model_b'],x['version_b'])])) for x in rows}
if written!=pairs or len(rows)!=len(pairs):raise ValueError('Distinct pair set differs')
for name,expected_models in [('models.jsonl',models),('additional_models.jsonl',models-base_models)]:
    with (a.inventory/name).open() as f:records=[json.loads(l) for l in f]
    found={(m['model_id'],str(m['version'])) for m in records}
    if found!=expected_models or len(records)!=len(found):raise ValueError('Model partition differs')
if len(models-base_models)!=r['additional_models'] or len(pairs-base_pairs)!=r['additional_model_pairs']:raise ValueError('Receipt count differs')
result=dict(status='passed_full_reference_comparison_ledger_readback',producer_receipt_sha256=sha(a.inventory/'receipt.json'),event_reference_comparisons=len(observed),unique_models=len(models),additional_models=len(models-base_models),distinct_model_pairs=len(pairs),additional_model_pairs=len(pairs-base_pairs),work_dispositions=counts,script_sha256=sha(__file__),scope='All eligible event x tied-reference x duplicate-side links independently enumerated from source tables. Exact model/pair work partitions, hashes and output scope checked. Does not independently redo source gene-to-model mapping, reference selection or raw coordinates.')
a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
