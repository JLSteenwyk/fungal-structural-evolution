#!/usr/bin/env python3
"""Catalog disjoint audited annotation shards after complete full-proteome search."""
import argparse,csv,json
from pathlib import Path
from prepare_pfam import digest


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['search','pfam-receipt','output']:ap.add_argument('--'+name,type=Path,required=True)
    ap.add_argument('--snapshots',nargs='+',type=Path,required=True);ap.add_argument('--readbacks',nargs='+',type=Path,required=True);a=ap.parse_args()
    if a.output.exists() or len(a.snapshots)!=len(a.readbacks):raise ValueError('Fresh output and paired sources required')
    sr=json.loads((a.search/'receipt.json').read_text());cfg=json.loads((a.search/'config.json').read_text());pfam=json.loads(a.pfam_receipt.read_text())
    if sr['status']!='complete_raw_domain_search' or sr['config_sha256']!=digest(a.search/'config.json') or cfg['pfam_receipt_sha256']!=digest(a.pfam_receipt):raise ValueError('Complete compatible raw search required')
    expected={Path(x['path']).stem:x for x in pfam['chunks']};raw={x['chunk']:x for x in sr['chunks']}
    if set(raw)!=set(expected) or len(sr['chunks'])!=len(expected):raise ValueError('Incomplete raw chunk universe')
    catalog=[];seen=set();sources=[]
    for folder,audit in zip(a.snapshots,a.readbacks):
        r=json.loads((folder/'receipt.json').read_text());ar=json.loads((audit/'receipt.json').read_text())
        if ar['status']!='passed_full_pfam_annotation_raw_field_readback' or ar['annotation_receipt_sha256']!=digest(folder/'receipt.json') or ar['rows_checked']!=r['annotated_hit_rows'] or ar['chunks_checked']!=r['completed_chunks']:raise ValueError('Incomplete/mismatched annotation readback')
        for item in r['shards']:
            key=item['chunk'];rp=folder/(key+'.annotation.receipt.json');shard=json.loads(rp.read_text());path=folder/(key+'.annotated.tsv.gz')
            if key in seen:raise ValueError('Overlapping annotation snapshots')
            seen.add(key)
            if digest(rp)!=item['receipt_sha256'] or digest(path)!=item['annotation_sha256']:raise ValueError('Changed annotation shard')
            if (shard['input_receipt_sha256']!=cfg['input_receipt_sha256'] or shard['pfam_receipt_sha256']!=cfg['pfam_receipt_sha256'] or shard['raw_table_sha256']!=raw[key]['table_sha256'] or raw[key]['profile_sha256']!=expected[key]['sha256'] or shard['rows']!=raw[key]['domain_hit_rows']):raise ValueError('Shard/raw lineage mismatch')
            catalog.append({'chunk':key,'profiles':expected[key]['profiles'],'hit_rows':shard['rows'],'annotation_path':str(path),'annotation_sha256':digest(path),'shard_receipt':str(rp),'shard_receipt_sha256':digest(rp),'raw_table_sha256':raw[key]['table_sha256']})
        sources.append({'snapshot':str(folder),'receipt_sha256':digest(folder/'receipt.json'),'readback':str(audit),'readback_sha256':digest(audit/'receipt.json')})
    if seen!=set(expected) or sum(x['hit_rows'] for x in catalog)!=sr['domain_hit_rows'] or sum(x['profiles'] for x in catalog)!=pfam['families']:raise ValueError('Incomplete full annotation universe')
    a.output.mkdir(parents=True)
    with (a.output/'annotation_shards.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(catalog[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(sorted(catalog,key=lambda x:x['chunk']))
    result={'status':'complete_audited_additional_full_proteome_annotation_catalog','chunks':len(catalog),'pfam_profiles':pfam['families'],'annotated_hit_rows':sr['domain_hit_rows'],'source_snapshots':sources,'raw_search_receipt_sha256':digest(a.search/'receipt.json'),'script_sha256':digest(Path(__file__)),'interpretation':'Complete disjoint catalog of additional full-proteome annotation shards. Marker query partition remains separate. Not overlap-resolved architectures or proof of biological absence. Source-partition E-values retained.','artifacts':{'annotation_shards.tsv':digest(a.output/'annotation_shards.tsv')}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
