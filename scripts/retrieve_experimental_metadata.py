#!/usr/bin/env python3
"""Freeze complete RCSB entity and entry metadata for the experimental candidate set."""
import argparse,concurrent.futures,datetime,fcntl,json,time,urllib.request
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--inventory',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    inv=checked_receipt(a.inventory);ids=(a.inventory/'polymer_entity_ids.txt').read_text().splitlines()
    if len(ids)!=len(set(ids)) or len(ids)!=inv['unique_polymer_entities']:raise ValueError('Inventory identities differ')
    entries=sorted({i.rsplit('_',1)[0] for i in ids})
    a.output.mkdir(parents=True,exist_ok=True);lock=(a.output/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    config={'inventory_receipt_sha256':sha(a.inventory/'receipt.json'),'script_sha256':sha(Path(__file__)),'endpoint':'https://data.rcsb.org/rest/v1/core','workers':2,'per_worker_pause_seconds':0.25,'entity_count':len(ids),'entry_count':len(entries)}
    cp=a.output/'config.json'
    if cp.exists() and json.loads(cp.read_text())!=config:raise ValueError('Changed retrieval configuration')
    cp.write_text(json.dumps(config,indent=2)+'\n')
    for kind in ['polymer_entity','entry']:(a.output/kind).mkdir(exist_ok=True)
    def fetch(job):
        kind,identifier=job;path=a.output/kind/(identifier+'.json');receipt=path.with_suffix('.receipt.json')
        suffix=identifier.replace('_','/') if kind=='polymer_entity' else identifier
        url=config['endpoint']+'/'+kind+'/'+suffix
        if receipt.exists():
            r=json.loads(receipt.read_text())
            if r['url']!=url or r['config_sha256']!=sha(cp) or r['response_sha256']!=sha(path):raise ValueError('Changed cached response')
        else:
            for attempt in range(4):
                try:
                    with urllib.request.urlopen(url,timeout=30) as response:
                        if response.status!=200:raise ValueError('Unexpected status')
                        raw=response.read()
                    data=json.loads(raw)
                    if data.get('rcsb_id')!=identifier:raise ValueError('Response identity differs')
                    break
                except Exception:
                    if attempt==3:raise
                    time.sleep(2**attempt)
            tmp=path.with_suffix('.partial');tmp.write_bytes(raw);tmp.replace(path)
            r={'status':'verified_identity_metadata_response','kind':kind,'identifier':identifier,'url':url,'config_sha256':sha(cp),'response_sha256':sha(path),'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
            tmp=receipt.with_suffix('.partial');tmp.write_text(json.dumps(r,indent=2)+'\n');tmp.replace(receipt);time.sleep(.25)
        if json.loads(path.read_text()).get('rcsb_id')!=identifier:raise ValueError('Cached identity differs')
        return {'kind':kind,'identifier':identifier,'receipt_sha256':sha(receipt),'response_sha256':sha(path)}
    jobs=[('polymer_entity',i) for i in ids]+[('entry',i) for i in entries];results=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        for index,row in enumerate(pool.map(fetch,jobs),1):
            results.append(row)
            if index%100==0:print(index,'/',len(jobs),'verified',flush=True)
    result={'status':'complete_experimental_candidate_metadata','config_sha256':sha(cp),'entity_count':len(ids),'entry_count':len(entries),'responses':results,'interpretation':'Metadata response identities and bytes verified. Sequence correspondence, observed coordinates, experimental quality, biological context and training overlap still require analysis.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print('Complete',len(results),flush=True)


if __name__=='__main__':main()
