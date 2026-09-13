#!/usr/bin/env python3
"""Retrieve all PDB coordinate archives nominated by exact deposited-sequence matches."""
import argparse,datetime,fcntl,gzip,json,re,shutil,time,urllib.request
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha,read_table


def verify_gzip(path,entry):
    total=0;header=b''
    with gzip.open(path,'rb') as f:
        while True:
            data=f.read(1024*1024)
            if not data:break
            if not header:header=data[:4096]
            total+=len(data)
    match=re.search(rb'(?m)^data_([^\s]+)',header)
    if not match or match.group(1).decode().upper()!=entry.upper():raise ValueError('Coordinate data-block identity differs')
    return total


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--screen',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    checked_receipt(a.screen);rows=read_table(a.screen/'sequence_correspondence.tsv')
    entries=sorted({r['entity_id'].rsplit('_',1)[0] for r in rows if r['sequence_class']=='exact_full_sequence'})
    if not entries or any(not re.fullmatch('[A-Za-z0-9]+',x) for x in entries):raise ValueError('Invalid entry inventory')
    a.output.mkdir(parents=True,exist_ok=True);lock=(a.output/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    config={'screen_receipt_sha256':sha(a.screen/'receipt.json'),'script_sha256':sha(Path(__file__)),'entries':entries,'endpoint':'https://files.rcsb.org/download/','selection':'Every entry with an exact full canonical deposited-sequence match in the frozen screen; no quality-based benchmark eligibility inferred.'}
    cp=a.output/'config.json'
    if cp.exists() and json.loads(cp.read_text())!=config:raise ValueError('Changed retrieval configuration')
    cp.write_text(json.dumps(config,indent=2)+'\n');results=[]
    for index,entry in enumerate(entries,1):
        path=a.output/(entry+'.cif.gz');rp=a.output/(entry+'.receipt.json');url=config['endpoint']+entry+'.cif.gz'
        if rp.exists():
            r=json.loads(rp.read_text())
            if r['config_sha256']!=sha(cp) or r['gzip_sha256']!=sha(path):raise ValueError('Changed cached coordinates')
            if verify_gzip(path,entry)!=r['uncompressed_bytes']:raise ValueError('Cached archive content changed')
        else:
            if shutil.disk_usage(a.output).free<20000000000:raise ValueError('Insufficient planned disk headroom')
            partial=path.with_suffix('.partial')
            for attempt in range(4):
                try:
                    with urllib.request.urlopen(url,timeout=30) as response,partial.open('wb') as f:
                        if response.status!=200:raise ValueError('Unexpected response status')
                        shutil.copyfileobj(response,f,1024*1024)
                    size=verify_gzip(partial,entry);partial.replace(path);break
                except Exception:
                    if attempt==3:raise
                    time.sleep(2**attempt)
            r={'status':'downloaded_gzip_integrity_and_block_identity_checked','entry_id':entry,'url':url,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'config_sha256':sha(cp),'gzip_sha256':sha(path),'compressed_bytes':path.stat().st_size,'uncompressed_bytes':size,'interpretation':'Archive CRC/readability and data-block ID verified. Atomic residues, entity/chain correspondence, experimental quality and benchmark eligibility not yet validated.'}
            temp=rp.with_suffix('.partial');temp.write_text(json.dumps(r,indent=2)+'\n');temp.replace(rp);time.sleep(.25)
        results.append({'entry_id':entry,'receipt_sha256':sha(rp),'gzip_sha256':r['gzip_sha256'],'compressed_bytes':r['compressed_bytes'],'uncompressed_bytes':r['uncompressed_bytes']});print(index,'/',len(entries),entry,r['compressed_bytes'],flush=True)
    receipt={'status':'complete_frozen_experimental_coordinate_download','config_sha256':sha(cp),'entries':len(entries),'compressed_bytes':sum(r['compressed_bytes'] for r in results),'uncompressed_bytes':sum(r['uncompressed_bytes'] for r in results),'results':results,'interpretation':'All nominated archives downloaded; candidate coverage follows the frozen screen, which may be partial. No experimental accuracy benchmark or coordinate-residue validation is completed by this stage.'}
    (a.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':main()
