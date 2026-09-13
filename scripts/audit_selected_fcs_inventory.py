#!/usr/bin/env python3
"""Audit all selected FCS reports, separately validating reports lacking publisher MD5."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
from retrieve_selected_fcs_reports import parse_report, fetch, sha, FIELDS


def interval_union(records, actions):
    # Independent event-sweep count in zero-based half-open coordinates.
    events=defaultdict(lambda:defaultdict(int))
    for r in records:
        if r['action'] in actions:
            events[r['seq_id']][r['start_pos']-1]+=1;events[r['seq_id']][r['end_pos']]-=1
    result=0
    for positions in events.values():
        previous=None;active=0
        for position,change in sorted(positions.items()):
            if active>0:result+=position-previous
            active+=change;previous=position
        if active!=0:raise ValueError('Unbalanced interval events')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inventory',type=Path,required=True);p.add_argument('--manifest',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable output')
    r=json.loads((a.inventory/'receipt.json').read_text())
    if r['manifest_sha256']!=sha(a.manifest):raise ValueError('Manifest changed')
    for name,digest in r['artifacts'].items():
        if sha(a.inventory/name)!=digest:raise ValueError('Inventory table changed')
    with a.manifest.open() as f:manifest={x['taxon_id']:x for x in csv.DictReader(f,delimiter='\t')}
    with (a.inventory/'flagged_regions.tsv').open() as f:old=list(csv.DictReader(f,delimiter='\t'))
    with (a.inventory/'taxon_summary.tsv').open() as f:old_summary={x['taxon_id']:x for x in csv.DictReader(f,delimiter='\t')}
    sources={x['taxon_id']:x for x in r['source_receipts']}
    if len(sources)!=len(r['source_receipts']) or set(sources)!=set(manifest) or set(old_summary)!=set(manifest):raise ValueError('Incomplete source grid')
    a.output.mkdir(parents=True);summaries=[];allrows=[];verified_original=[];repeats=[]
    for taxon in sorted(manifest):
        source=sources[taxon];row=manifest[taxon];acc=row['assembly_accession'];status=source['status'];records=[];header=None;report_sha='';path=''
        if source['assembly_accession']!=acc:raise ValueError('Assembly mismatch')
        if status!='external_source_no_ncbi_report_requested':
            folder=a.inventory/'reports'/acc
            for name,digest in source['saved_artifacts'].items():
                if sha(folder/name)!=digest:raise ValueError('Saved raw artifact changed')
            url=row['proteome_url'].removesuffix('_protein.faa.gz')+'_fcs_report.txt'
            if source['url']!=url or not url.rsplit('/',1)[1].startswith(acc+'_'):raise ValueError('Source URL mismatch')
            report=folder/url.rsplit('/',1)[1]
            if status=='publisher_verified_report':
                header,records=parse_report(report.read_text());path=str(report);report_sha=sha(report)
                if header!=source['header'] or report_sha!=source['report_sha256'] or hashlib.md5(report.read_bytes()).hexdigest()!=source['publisher_md5']:raise ValueError('Publisher report provenance mismatch')
                verified_original.extend(dict(taxon_id=taxon,assembly_accession=acc,**{k:str(v) for k,v in item.items()}) for item in records)
            elif report.exists() and (folder/'md5checksums.txt').exists():
                checks=[s.split()[0] for s in (folder/'md5checksums.txt').read_text().splitlines() if len(s.split())==2 and s.split()[1].removeprefix('./')==report.name]
                if not checks:
                    # Lack of publisher checksum is not a checksum mismatch.
                    try:
                        raw=fetch(url);target=a.output/(acc+'_repeat_fcs_report.txt');target.write_bytes(raw)
                        if raw!=report.read_bytes():raise ValueError('Repeated retrieval differs from original')
                        header,records=parse_report(raw.decode());status='repeat_identical_https_report_without_publisher_md5';path=str(report);report_sha=sha(report)
                        repeats.append({'taxon_id':taxon,'url':url,'original_path':str(report),'repeat_path':str(target),'sha256':sha(target),'publisher_md5_available':False})
                    except Exception as exc:
                        status='repeat_or_parse_unresolved';repeats.append({'taxon_id':taxon,'url':url,'error':str(exc)});records=[]
                elif len(checks)!=1 or hashlib.md5(report.read_bytes()).hexdigest()!=checks[0]:status='publisher_checksum_mismatch_unresolved'
        validated=status in ['publisher_verified_report','repeat_identical_https_report_without_publisher_md5']
        info={'exclude_fix_trim_union_bp':interval_union(records,{'EXCLUDE','FIX','TRIM'}),'review_rare_union_bp':interval_union(records,{'REVIEW_RARE'}),'review_union_bp':interval_union(records,{'REVIEW'})}
        if source['status']=='publisher_verified_report':
            for key,value in info.items():
                if int(old_summary[taxon][key])!=value:raise ValueError('Independent interval union differs')
            if len(records)!=int(old_summary[taxon]['report_rows']) or len({x['seq_id'] for x in records})!=int(old_summary[taxon]['flagged_sequence_count']):raise ValueError('Count readback differs')
        summaries.append({'taxon_id':taxon,'assembly_accession':acc,'validation_status':status,'report_rows':len(records) if validated else '', 'flagged_sequence_count':len({x['seq_id'] for x in records}) if validated else '',**{k:v if validated else '' for k,v in info.items()},'reported_actions':json.dumps(dict(Counter(x['action'] for x in records)),sort_keys=True),'report_sha256':report_sha,'report_path':path,'fcs_run_metadata':json.dumps(header) if validated else ''})
        allrows.extend(dict(taxon_id=taxon,assembly_accession=acc,validation_status=status,**item) for item in records)
    sortkey=lambda x:json.dumps(x,sort_keys=True)
    if sorted(old,key=sortkey)!=sorted(verified_original,key=sortkey):raise ValueError('Original full flagged-region table differs from raw report readback')
    for name,data,fields in [('taxon_summary.tsv',summaries,list(summaries[0])),('flagged_regions.tsv',allrows,['taxon_id','assembly_accession','validation_status']+FIELDS)]:
        with (a.output/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    result={'status':'complete_fcs_inventory_raw_readback_and_missing_checksum_review','source_receipt_sha256':sha(a.inventory/'receipt.json'),'manifest_sha256':sha(a.manifest),'script_sha256':sha(Path(__file__)),'parser_script_sha256':sha(Path(__file__).with_name('retrieve_selected_fcs_reports.py')),'taxa':len(summaries),'validation_counts':dict(Counter(x['validation_status'] for x in summaries)),'validated_reports_with_rows':sum(isinstance(x['report_rows'],int) and x['report_rows']>0 for x in summaries),'taxa_with_exclude_fix_trim_regions':sum(isinstance(x['exclude_fix_trim_union_bp'],int) and x['exclude_fix_trim_union_bp']>0 for x in summaries),'flagged_region_rows':len(allrows),'original_raw_rows_compared':len(old),'repeat_retrievals':repeats,'artifacts':{p.name:sha(p) for p in a.output.iterdir() if p.is_file()},'interpretation':'Exact-version publisher reports; missing MD5 separately labeled and validated by identical repeated HTTPS retrieval. Full original region-table readback and independent interval unions passed. Review-only and organelle records remain distinct; no contamination-free certificate, new FCS run, biological contamination confirmation or sequence removal. Gene overlaps remain pending.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['repeat_retrievals','artifacts']},indent=2))


if __name__=='__main__':main()
