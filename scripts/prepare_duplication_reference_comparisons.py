#!/usr/bin/env python3
"""Inventory every tied provisional reference and additional structural-comparison work."""
import argparse,csv,hashlib,json,sqlite3
from collections import Counter
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha


def pair_identity(a,b):
    key=tuple(sorted([a,b]));return key,hashlib.sha256(json.dumps(key,separators=(',',':')).encode()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['inventory','inventory-plan','base-queue','catalog','output']:p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    rp=a.inventory/'receipt.json';r=json.loads(rp.read_text());plan=json.loads(a.inventory_plan.read_text());base=json.loads((a.base_queue/'receipt.json').read_text());cat=json.loads((a.catalog/'receipt.json').read_text())
    if r['status']!='complete_duplication_sister_reference_inventory' or r['plan_sha256']!=sha(a.inventory_plan):raise ValueError('Reference inventory mismatch')
    pins={str(rp):sha(rp),str(a.inventory_plan):sha(a.inventory_plan),str(a.base_queue/'receipt.json'):sha(a.base_queue/'receipt.json'),str(a.catalog/'receipt.json'):sha(a.catalog/'receipt.json')}
    bridge=Path(plan['bridge']);pins[str(bridge)]=sha(bridge)
    if pins[str(bridge)]!=plan['pins'][str(bridge)]:raise ValueError('Reference bridge changed')
    brp=bridge.parent/'receipt.json';br=json.loads(brp.read_text());pins[str(brp)]=sha(brp)
    if br['catalog_receipt_sha256']!=sha(a.catalog/'receipt.json'):raise ValueError('Different frozen catalog')
    events=[];needed=set();counts=Counter()
    for guide in ['profile','mafft']:
        path=a.inventory/(guide+'_sister_references.tsv');pins[str(path)]=sha(path)
        if pins[str(path)]!=r['artifacts'][path.name]:raise ValueError('Changed reference table')
        with path.open() as f:
            for row in csv.DictReader(f,delimiter='\t'):
                if row['status']!='provisional_reference_available':continue
                ties=json.loads(row['nearest_reference_genes'])
                if not ties or ties!=sorted(set(ties)) or row['chosen_reference_gene']!=ties[0]:raise ValueError('Invalid reference ties')
                events.append((guide,row,ties));needed.update([row['gene_a'],row['gene_b'],*ties]);counts[guide+'_eligible_events']+=1;counts[guide+'_event_reference_links']+=len(ties)
    con=sqlite3.connect('file:'+str(bridge.resolve())+'?mode=ro',uri=True);links={}
    for taxon,protein,seq,model,version,path in con.execute('SELECT taxon_id,protein_id,sequence_sha256,model_id,version,model_path FROM structures'):
        gene=taxon+'_'+protein
        if gene in needed:links[gene]=(model,version,seq,path)
    con.close()
    if set(links)!=needed:raise ValueError('Unmapped reference/focal protein')
    base_models={};path=a.base_queue/'models.jsonl';pins[str(path)]=sha(path)
    if pins[str(path)]!=base['artifacts'][path.name]:raise ValueError('Changed base model catalog')
    with path.open() as f:
        for line in f:
            m=json.loads(line);base_models[m['model_id'],m['version']]=m
    path=a.base_queue/'model_pairs.tsv';pins[str(path)]=sha(path)
    if pins[str(path)]!=base['artifacts'][path.name]:raise ValueError('Changed base pairs')
    with path.open() as f:base_pairs={tuple(sorted([(x['model_a'],int(x['version_a'])),(x['model_b'],int(x['version_b']))])) for x in csv.DictReader(f,delimiter='\t')}
    a.output.mkdir(parents=True);pairs={};models_needed=set()
    with (a.output/'event_reference_comparisons.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['guide','family','gene_node','gene_a','gene_b','reference_gene','lexical_representative','focal_side','focal_gene','focal_model','focal_version','reference_model','reference_version','pair_key','work_disposition'])
        for guide,row,ties in events:
            chosen=links[row['chosen_reference_gene']]
            if (chosen[0],str(chosen[1]))!=(row['reference_model'],row['reference_version']):raise ValueError('Chosen reference model differs')
            for g,k in [(row['gene_a'],'model_a'),(row['gene_b'],'model_b')]:
                if links[g][0]!=row[k]:raise ValueError('Focal model differs')
            for ref in ties:
                rm=links[ref][:2];models_needed.add(rm)
                for side in ['a','b']:
                    gene=row['gene_'+side];fm=links[gene][:2];models_needed.add(fm);key,digest=pair_identity(fm,rm)
                    status='identical_model' if fm==rm else ('existing_duplicate_pair' if key in base_pairs else 'additional_pair')
                    counts[status+'_event_links']+=1
                    if fm!=rm:pairs[key]=digest
                    w.writerow([guide,row['family'],row['gene_node'],row['gene_a'],row['gene_b'],ref,int(ref==row['chosen_reference_gene']),side,gene,*fm,*rm,digest,status])
    selected={};path=a.catalog/'models.jsonl';pins[str(path)]=sha(path)
    if pins[str(path)]!=cat['artifacts'][path.name]:raise ValueError('Changed original catalog')
    with path.open() as f:
        for line in f:
            m=json.loads(line);key=(m['model_id'],m['version'])
            if key in models_needed:
                if key in selected:raise ValueError('Repeated catalog model')
                selected[key]=m
    if set(selected)!=models_needed:raise ValueError('Missing frozen model')
    for model,version,seq,path in links.values():
        m=selected[model,version]
        if (m['sequence_sha256'],m['path'])!=(seq,path):raise ValueError('Sequence/path identity differs')
        if (model,version) in base_models and m!=base_models[model,version]:raise ValueError('Base catalog descriptor differs')
    extras=set(selected)-set(base_models);additional=set(pairs)-base_pairs
    for name,keys in [('models.jsonl',set(selected)),('additional_models.jsonl',extras)]:
        with (a.output/name).open('w') as f:
            for key in sorted(keys):f.write(json.dumps(selected[key],separators=(',',':'))+'\n')
    with (a.output/'model_pairs.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['pair_key','model_a','version_a','model_b','version_b','work_disposition'])
        for key,digest in sorted(pairs.items()):w.writerow([digest,*key[0],*key[1],'additional_pair' if key in additional else 'existing_duplicate_pair'])
    byte_count=sum(Path(selected[key]['path']).stat().st_size for key in extras)
    for path,h in pins.items():
        if sha(path)!=h:raise ValueError('Source changed during inventory')
    result=dict(status='complete_provisional_reference_comparison_inventory',source_pins=pins,counts=dict(counts),all_reference_comparison_models=len(selected),additional_models=len(extras),additional_coordinate_bytes=byte_count,unique_distinct_model_pairs=len(pairs),additional_model_pairs=len(additional),existing_duplicate_model_pairs=len(set(pairs)&base_pairs),additional_directed_mask_dispositions=4*len(additional),script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in a.output.iterdir()},scope='Both duplicate copies compared with every tied nearest provisional reference under each guide; all event links retained, model computations deduplicated. Existing frozen queue unchanged. Model/sequence/version/path provenance bound; extra coordinate bytes measured by stat, not validated. Reference orthology, complete sister-choice/path readback, additional coordinate validation, alignments and asymmetry tests remain pending.')
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['source_pins','artifacts']},indent=2))

if __name__=='__main__':main()
