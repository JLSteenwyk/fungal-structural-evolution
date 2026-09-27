#!/usr/bin/env python3
"""Resolve qualified source-overlap cells to full encoded protein/model identities."""
import argparse,csv,json,hashlib
from collections import Counter
from pathlib import Path
import numpy as np
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['comparison','reference','local','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();cr=json.loads((a.comparison/'receipt.json').read_text())
    if cr['status']!='complete_paired_source_state_comparison':raise ValueError('Incomplete source comparison')
    pins=dict(cr['source_sha256']);pins[str(a.comparison/'receipt.json')]=sha(a.comparison/'receipt.json')
    table=a.comparison/'taxon_marker_state_comparison.tsv'
    if sha(table)!=cr['artifacts'][table.name]:raise ValueError('Changed comparison table')
    pins[str(table)]=sha(table)
    with table.open() as f:cells=list(csv.DictReader(f,delimiter='\t'))
    targets={(r['marker'],r['taxon']) for r in cells}
    if len(targets)!=len(cells):raise ValueError('Repeated comparison cell')
    indexes=[]
    for root in [a.reference,a.local]:
        pr=json.loads((root/'receipt.json').read_text())
        if str(root/'receipt.json') not in pins or sha(root/'receipt.json')!=pins[str(root/'receipt.json')]:raise ValueError('Different paired source')
        source={}
        for label in ['snapshot','encodings']:
            d=pr['source_receipts'][label];folder=Path(d['path']);rp=folder/'receipt.json'
            if sha(rp)!=d['sha256']:raise ValueError('Changed upstream receipt')
            pins[str(rp)]=sha(rp);source[label]=(folder,json.loads(rp.read_text()))
        folder,receipt=source['snapshot'];lp=folder/'marker_structure_links.tsv'
        if sha(lp)!=receipt['artifacts'][lp.name]:raise ValueError('Changed model links')
        pins[str(lp)]=sha(lp)
        links={}
        for row in csv.DictReader(lp.open(),delimiter='\t'):
            key=row['marker'],row['taxon_id']
            if key in targets:
                if key in links:raise ValueError('Ambiguous source model')
                links[key]=row
        if set(links)!=targets:raise ValueError('Missing overlap model')
        folder,receipt=source['encodings'];ep=folder/'model_summary.tsv'
        if sha(ep)!=receipt['artifacts'][ep.name]:raise ValueError('Changed encoding index')
        pins[str(ep)]=sha(ep);encodings={}
        for row in csv.DictReader(ep.open(),delimiter='\t'):
            key=row['model_id'],row['version']
            if key in encodings:raise ValueError('Repeated encoding identity')
            encodings[key]=row
        selected={}
        for key,link in links.items():
            encoding=encodings[link['model_id'],link['model_version']];path=Path(encoding['encoding_path'])
            if sha(path)!=encoding['encoding_sha256']:raise ValueError('Changed encoded model')
            with np.load(path,allow_pickle=False) as data:sequence=str(data['sequence'])
            digest=hashlib.sha256(sequence.encode()).hexdigest()
            if digest!=encoding['sequence_sha256'] or digest!=link['sequence_sha256'] or len(sequence)!=int(encoding['length']):raise ValueError('Full sequence identity mismatch')
            mp=Path(link['model_path'])
            if sha(mp)!=link['model_sha256']:raise ValueError('Changed coordinate model')
            pins[str(path)]=sha(path);pins[str(mp)]=sha(mp)
            selected[key]=(link,encoding,sequence)
        indexes.append(selected)
    a.output.mkdir(parents=True,exist_ok=False);rows=[];counts=Counter();sites=Counter();mismatches=Counter();modelpairs=set()
    for cell in cells:
        key=cell['marker'],cell['taxon'];left,right=[x[key] for x in indexes]
        same=left[2]==right[2];status='identical_complete_encoded_sequence' if same else 'different_complete_encoded_sequence'
        row=dict(cell,full_sequence_status=status,same_protein_identifier=int(left[0]['protein_id']==right[0]['protein_id']))
        pair=[]
        for prefix,(link,encoding,sequence) in zip(['reference','local'],[left,right]):
            for name in ['protein_id','sequence_sha256','model_id','model_version','model_provider','model_tool','model_path','model_sha256']:row[prefix+'_'+name]=link[name]
            row[prefix+'_length']=len(sequence);row[prefix+'_encoding_path']=encoding['encoding_path'];row[prefix+'_encoding_sha256']=encoding['encoding_sha256'];pair.append((link['model_id'],link['model_version']))
        modelpairs.add(tuple(pair));counts[status]+=1;sites[status]+=int(cell['same_aa_positions']);mismatches[status]+=int(cell['state_mismatches']);rows.append(row)
    path=a.output/'source_model_context.tsv'
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    for path,h in pins.items():
        if sha(path)!=h:raise ValueError('Changed provenance during comparison')
    result=dict(status='complete_paired_source_full_sequence_context',source_sha256=pins,script_sha256=sha(__file__),cells=len(rows),distinct_model_pairs=len(modelpairs),counts=dict(counts),observed_positions_by_status=dict(sites),state_mismatches_by_status=dict(mismatches),artifacts={'source_model_context.tsv':sha(a.output/'source_model_context.tsv')},scope='All overlap cells joined to selected source models and full native-encoding sequences; exact sequences/lengths and encoding/coordinate file hashes checked. Complete encoded sequence identity does not prove biological gene-copy orthology or equal prediction settings. State disagreement remains context/confidence-selected, not experimental error or a calibrated evolutionary effect.')
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['source_sha256','artifacts']},indent=2))

if __name__=='__main__':main()
