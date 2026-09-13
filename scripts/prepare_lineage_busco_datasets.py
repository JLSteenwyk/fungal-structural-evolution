#!/usr/bin/env python3
"""Acquire the explicitly selected BUSCO 12.2 datasets from a frozen publisher index."""
import argparse
import csv
import hashlib
import json
import tarfile
import urllib.request
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

ROOT=Path(__file__).resolve().parents[1]
NAMES=['fungi','ascomycota','basidiomycota','chytridiomycota','microsporidia','mucoromycota']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--resume',action='store_true');a=p.parse_args()
    a.output=a.output.resolve()
    if a.output.exists() and not a.resume:raise FileExistsError('Use new immutable dataset collection or explicit resume')
    if (a.output/'receipt.json').exists():raise FileExistsError('Completed collection is immutable')
    index=ROOT/'data/busco_downloads/file_versions.tsv'
    with index.open() as f:versions={r[0]:r for r in csv.reader(f,delimiter='\t')}
    a.output.mkdir(parents=True,exist_ok=a.resume)
    frozen=a.output/'publisher_index.tsv'
    if frozen.exists() and sha(frozen)!=sha(index):raise ValueError('Publisher index changed since partial acquisition')
    if not frozen.exists():frozen.write_bytes(index.read_bytes())
    def get(short):
        name=short+'_odb12.2';v=versions[name]
        if v[4]!='lineages':raise ValueError('Unexpected dataset category')
        url=f'https://busco-data.ezlab.org/v5/data/lineages/{name}.{v[1]}.tar.gz';archive=a.output/(name+'.tar.gz')
        if archive.exists():
            if hashlib.md5(archive.read_bytes()).hexdigest()!=v[2]:raise ValueError('Existing archive is incomplete or changed; preserve and review before replacement')
        else:
            with urllib.request.urlopen(url,timeout=120) as response,archive.open('xb') as out:
                while block:=response.read(1024*1024):out.write(block)
            if hashlib.md5(archive.read_bytes()).hexdigest()!=v[2]:raise ValueError('Publisher archive MD5 mismatch')
        with tarfile.open(archive,'r:gz') as tar:
            for m in tar.getmembers():
                if Path(m.name).parts[0]!=name or m.issym() or m.islnk() or not (m.isdir() or m.isfile()):raise ValueError('Unexpected archive member')
            tar.extractall(a.output,filter='data')
        base=a.output/name;cfg={k:v for k,v in (line.split('=',1) for line in (base/'dataset.cfg').read_text().splitlines() if '=' in line)}
        if cfg['name']!=name or cfg['creation_date']!=v[1]:raise ValueError('Archive version differs')
        files={str(x.relative_to(base)):sha(x) for x in sorted(base.rglob('*')) if x.is_file()}
        if len(list((base/'hmms').glob('*.hmm')))!=int(cfg['number_of_BUSCOs']):raise ValueError('HMM marker count differs')
        return {'dataset':name,'url':url,'publisher_md5':v[2],'archive_sha256':sha(archive),'archive_bytes':archive.stat().st_size,'path':str(base.relative_to(ROOT)),'config':cfg,'files_sha256':files}
    with ThreadPoolExecutor(max_workers=2) as pool:datasets=list(pool.map(get,NAMES))
    r={'status':'complete_pinned_lineage_busco_dataset_collection','publisher_index_sha256':sha(a.output/'publisher_index.tsv'),'script_sha256':sha(Path(__file__)),'datasets':datasets,'interpretation':'Explicit taxonomic QC panels from the same OrthoDB release. Dataset recovery percentages have different marker denominators and are not interchangeable completeness measurements.'}
    (a.output/'receipt.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps([{k:x[k] for k in ['dataset','archive_bytes','config']} for x in datasets],indent=2))


if __name__=='__main__':main()
