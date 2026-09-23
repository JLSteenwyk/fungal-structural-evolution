#!/usr/bin/env python3
"""Rebuild a complete PAE manifest from unchanged verified records and repaired cache entries."""
import argparse,json
from pathlib import Path
from catalog_whole_proteome_structures import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    for name,h in plan['pins'].items():
        if sha(name)!=h:raise ValueError('Changed input: '+name)
    old=Path(plan['prior']);prior=json.loads((old/'receipt.json').read_text());rows=json.loads((old/'pae_manifest.json').read_text());catalog=Path(plan['catalog']);cr=json.loads((catalog/'receipt.json').read_text())
    if sha(old/'pae_manifest.json')!=prior['artifacts']['pae_manifest.json'] or prior['mapping_receipt_sha256']!=sha(catalog/'receipt.json'):raise ValueError('Prior binding differs')
    models=json.loads((catalog/'model_provenance.json').read_text());expected={(m['model_id'],m['version']):m for m in models}
    if len(expected)!=len(models) or len(rows)!=len(models) or {(r['model_id'],r['version']) for r in rows}!=set(expected):raise ValueError('Incomplete source grid')
    records=[];repaired=[]
    for previous in rows:
        key=previous['model_id'],previous['version'];m=expected[key];cache=Path('data/structures/pae')/f'{key[0]}-v{key[1]}.receipt.json';r=json.loads(cache.read_text())
        if r['status']!='verified' or any(r[k]!=m[k] for k in ['model_id','version','length','sequence_sha256']) or r['url']!=m['pae_url']:raise ValueError('Cache provenance differs')
        if previous['status']=='verified':
            if r!=previous:raise ValueError('Prior verified record changed')
        else:repaired.append(key)
        data=Path(r['path'])
        if data.stat().st_size!=r['compressed_bytes'] or sha(data)!=r['gzip_sha256']:raise ValueError('Changed PAE bytes')
        records.append(r)
    if len(repaired)!=prior['models_failed'] or len(records)!=cr['distinct_models']:raise ValueError('Recovery scope differs')
    target=Path(plan['output']);target.mkdir(parents=True,exist_ok=False);records.sort(key=lambda r:(r['model_id'],r['version']));manifest=target/'pae_manifest.json';manifest.write_text(json.dumps(records,indent=2)+'\n')
    result={'mapping_snapshot':str(catalog),'mapping_receipt_sha256':sha(catalog/'receipt.json'),'models_requested':len(models),'models_verified':len(records),'models_failed':0,'compressed_bytes':sum(r['compressed_bytes'] for r in records),'json_bytes':sum(r['json_bytes'] for r in records),'script_sha256':sha(__file__),'prior_receipt_sha256':sha(old/'receipt.json'),'repaired_models':repaired,'artifacts':{'pae_manifest.json':sha(manifest)},'scope':'Every compressed matrix rehashed and exact model/version/sequence/length/URL/cache receipt verified. Previous verified records preserved unchanged; repaired record passed original retrieval matrix validation. Prior dimensional/numerical validation reused for unchanged matrices; complete mapping-bound context validation remains downstream.'}
    (target/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    for name,h in plan['pins'].items():
        if sha(name)!=h:raise ValueError('Input changed during recovery')
    if sha(a.plan)!=ph:raise ValueError('Plan changed')
    control=Path(plan['control']);control.mkdir(parents=True,exist_ok=False)
    receipt={'status':'complete_refreshed_catalog_pae_prefetch_pending_mapping_binding','plan_sha256':ph,'stage':'pae','result_receipt':str(target/'receipt.json'),'result_receipt_sha256':sha(target/'receipt.json'),'recovered_models':len(repaired)}
    (control/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
