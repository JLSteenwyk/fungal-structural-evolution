#!/usr/bin/env python3
"""Annotate a pinned snapshot of completed full-proteome Pfam chunks by streaming."""
import argparse,csv,gzip,json,time
from pathlib import Path
from Bio import SeqIO
from prepare_pfam import digest,ROOT
from summarize_marker_domains import read_metadata,parse_hit


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args()
    plan=json.loads(a.plan.read_text());out=Path(plan['output'])
    if out.exists():raise FileExistsError('Fresh immutable snapshot required')
    for file,expected in plan['pins'].items():
        if digest(Path(file))!=expected:raise ValueError('Changed source: '+file)
    inputs=Path(plan['inputs']);search=Path(plan['search']);ir=json.loads((inputs/'receipt.json').read_text());cfg=json.loads((search/'config.json').read_text());pfam=json.loads(Path(plan['pfam_receipt']).read_text())
    if cfg['input_receipt_sha256']!=digest(inputs/'receipt.json') or cfg['pfam_receipt_sha256']!=digest(Path(plan['pfam_receipt'])):raise ValueError('Search input/Pfam lineage differs')
    for name,expected in ir['artifacts'].items():
        if digest(inputs/name)!=expected:raise ValueError('Changed input artifact '+name)
    metadata_path=ROOT/'data/pfam/38.2/Pfam-A.hmm.dat';meta_record=next(x for x in pfam['files'] if x['uncompressed_file']==metadata_path.name)
    if digest(metadata_path)!=meta_record['uncompressed_sha256']:raise ValueError('Changed Pfam metadata')
    metadata=read_metadata(metadata_path)
    if len(metadata)!=pfam['families']:raise ValueError('Wrong metadata family universe')
    lengths={};residues=0
    for record in SeqIO.parse(inputs/'additional_sequences.faa','fasta'):
        if record.id in lengths:raise ValueError('Duplicate sequence')
        import hashlib
        sequence=str(record.seq)
        if record.id!='S'+hashlib.sha256(sequence.encode()).hexdigest():raise ValueError('Query sequence identity mismatch')
        lengths[record.id]=len(sequence);residues+=len(sequence)
    if len(lengths)!=ir['additional_unique_sequences'] or residues!=ir['additional_unique_residues']:raise ValueError('Query sequence/residue counts differ')
    expected={Path(c['path']).stem:c for c in pfam['chunks']}
    if len({x['chunk'] for x in plan['chunks']})!=len(plan['chunks']):raise ValueError('Duplicate planned chunk')
    out.mkdir(parents=True);records=[]
    for item in plan['chunks']:
        key=item['chunk'];rp=search/(key+'.receipt.json');raw=search/(key+'.domtblout');r=json.loads(rp.read_text())
        if digest(rp)!=item['receipt_sha256'] or r['profile_sha256']!=expected[key]['sha256'] or r['config_sha256']!=digest(search/'config.json') or digest(raw)!=r['table_sha256']:raise ValueError('Changed/incomplete raw chunk '+key)
        if '--cut_ga' not in r['command']:raise ValueError('Gathering thresholds required')
        dest=out/(key+'.annotated.tsv.gz');count=0;writer=None;start=time.monotonic()
        with gzip.open(dest,'wt',newline='') as target,raw.open() as source:
            for line in source:
                if not line.strip() or line.startswith('#'):continue
                hit=parse_hit(line,metadata,lengths)
                hit={'hit_id':f'FULL-{key}-{count+1:09d}','search_partition':'additional_full_proteome',**hit}
                if writer is None:writer=csv.DictWriter(target,fieldnames=list(hit),delimiter='\t',lineterminator='\n');writer.writeheader()
                writer.writerow(hit);count+=1
        if count!=r['domain_hit_rows']:raise ValueError('Raw hit count differs')
        result={'status':'complete_full_pfam_annotation_shard','chunk':key,'rows':count,'raw_chunk_receipt_sha256':digest(rp),'raw_table_sha256':digest(raw),'input_receipt_sha256':digest(inputs/'receipt.json'),'pfam_receipt_sha256':digest(Path(plan['pfam_receipt'])),'script_sha256':digest(Path(__file__)),'parser_sha256':digest(ROOT/'scripts/summarize_marker_domains.py'),'artifacts':{dest.name:digest(dest)}}
        shard=out/(key+'.annotation.receipt.json');shard.write_text(json.dumps(result,indent=2)+'\n');records.append({'chunk':key,'rows':count,'receipt_sha256':digest(shard),'annotation_sha256':digest(dest)})
        print(key,count,round(time.monotonic()-start,1),flush=True)
    result={'status':'complete_pinned_full_pfam_annotation_snapshot','full_search_complete':len(records)==len(expected),'completed_chunks':len(records),'planned_search_chunks':len(expected),'missing_chunks':sorted(set(expected)-{x['chunk'] for x in records}),'annotated_hit_rows':sum(x['rows'] for x in records),'additional_query_sequences_checked':len(lengths),'additional_query_residues_checked':residues,'plan_sha256':digest(a.plan),'script_sha256':digest(Path(__file__)),'shards':records,'interpretation':'Streaming annotation of explicitly pinned completed full-proteome raw search chunks. Search-partition E-values retained without global rescaling. Pfam types, GA scores, HMM/sequence boundaries and partial hits preserved. No overlap resolution, complete-proteome coverage, confirmed domain absence, architecture or evolutionary event claim. Marker query partition remains separate.','artifacts':{p.name:digest(p) for p in out.glob('*.annotation.receipt.json')}}
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print('snapshot complete',result['annotated_hit_rows'],flush=True)

if __name__=='__main__':main()
