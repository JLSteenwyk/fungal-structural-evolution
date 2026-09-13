#!/usr/bin/env python3
"""Independently enumerate region/CDS overlaps and verify protein/marker joins."""
import argparse
from collections import defaultdict, Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path
from urllib.parse import unquote
import numpy as np


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def table(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def union(intervals):
    events=defaultdict(int)
    for lo,hi in intervals:events[lo-1]+=1;events[hi]-=1
    active=0;previous=0;total=0
    for position in sorted(events):
        if active:total+=position-previous
        active+=events[position];previous=position
    return total


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['mapping','fcs','annotations','representatives','markers','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use new immutable audit output')
    r=json.loads((a.mapping/'receipt.json').read_text())
    for key,path in [('fcs',a.fcs/'receipt.json'),('annotations',a.annotations),('representatives',a.representatives)]:
        if r['source_receipts'][key]!=sha(path):raise ValueError('Changed source')
    if sha(a.markers)!=r['marker_mapping_sha256']:raise ValueError('Changed markers')
    for name,digest in r['artifacts'].items():
        if sha(a.mapping/name)!=digest:raise ValueError('Changed mapping artifact')
    regions=table(a.mapping/'region_annotation_review.tsv');original=table(a.fcs/'flagged_regions.tsv')
    if len(regions)!=len(original):raise ValueError('Region grid differs')
    bytax=defaultdict(list)
    original_ids={'FCS'+str(i+1).zfill(6):x for i,x in enumerate(original)}
    for reg in regions:
        if any(reg[k]!=str(v) for k,v in original_ids[reg['region_id']].items()):raise ValueError('Region fields changed')
        bytax[reg['taxon_id']].append(reg)
    features=table(a.mapping/'cds_feature_overlaps.tsv');actual={}
    for x in features:
        key=x['taxon_id'],x['region_id'],int(x['gff_line'])
        if key in actual:raise ValueError('Duplicate feature intersection')
        actual[key]=x
    sources={x['taxon_id']:x for x in json.loads(a.annotations.read_text())};seen=set();scan_count=0
    for taxon,regs in sorted(bytax.items()):
        source=sources[taxon];path=Path(source['path'])
        if sha(path)!=source['sha256']:raise ValueError('Changed annotation')
        seqids=set();byseq=defaultdict(list);lines={}
        with gzip.open(path,'rt') as f:
            for number,line in enumerate(f,1):
                if line.startswith('##FASTA'):break
                if line.startswith('##sequence-region '):seqids.add(line.split()[1])
                if line.startswith('#') or not line.strip():continue
                fields=line.rstrip('\n').split('\t');seqids.add(fields[0])
                if fields[2]=='CDS':byseq[fields[0]].append((int(fields[3]),int(fields[4]),number));lines[number]=fields;scan_count+=1
        for reg in regs:
            expected_state='represented_in_gff' if reg['seq_id'] in seqids else 'not_represented_in_gff'
            if reg['gff_sequence_status']!=expected_state:raise ValueError('GFF sequence presence differs')
            arr=np.array(byseq[reg['seq_id']],dtype=np.int64).reshape((-1,3));lo=int(reg['start_pos']);hi=int(reg['end_pos'])
            lengths=np.maximum(0,np.minimum(arr[:,1],hi)-np.maximum(arr[:,0],lo)+1);selected=np.flatnonzero(lengths>0);intervals=[];pids=set()
            for j in selected:
                start,end,number=map(int,arr[j]);key=taxon,reg['region_id'],number
                if key not in actual:raise ValueError('Missed region/CDS intersection')
                x=actual[key];fields=lines[number];seen.add(key);expected=(max(start,lo),min(end,hi));intervals.append(expected)
                if int(x['overlap_bp'])!=int(lengths[j]) or (int(x['intersection_start']),int(x['intersection_end']))!=expected or int(x['cds_start'])!=start or int(x['cds_end'])!=end or x['strand']!=fields[6] or x['phase']!=fields[7]:raise ValueError('CDS intersection fields differ')
                attrs={}
                for field in fields[8].split(';'):
                    if '=' in field:
                        k,v=field.split('=',1);attrs[k]=[unquote(z) for z in v.split(',')]
                if json.loads(x['protein_ids_json'])!=attrs.get('protein_id',[]):raise ValueError('Protein identity differs')
                pids.update(attrs.get('protein_id',[]))
            if len(selected)!=int(reg['overlapping_cds_feature_rows']) or union(intervals)!=int(reg['overlapping_cds_union_bp']) or sorted(pids)!=json.loads(reg['overlapping_protein_ids_json']):raise ValueError('Region summary differs')
    if seen!=set(actual):raise ValueError('Spurious CDS overlap')
    byprotein=defaultdict(list)
    for x in features:
        for pid in json.loads(x['protein_ids_json']):byprotein[x['taxon_id'],pid].append(x)
    proteins=table(a.mapping/'protein_overlap_review.tsv');reps={x['taxon_id']:x for x in json.loads(a.representatives.read_text())['taxa']};decisions={}
    if {(x['taxon_id'],x['protein_id']) for x in proteins}!=set(byprotein) or len(proteins)!=len(byprotein):raise ValueError('Protein grid differs')
    for x in proteins:
        t=x['taxon_id'];pid=x['protein_id'];hits=byprotein[t,pid]
        if t not in decisions:
            path=Path(reps[t]['path']).with_suffix('.decisions.tsv')
            if sha(path)!=reps[t]['decisions_sha256']:raise ValueError('Decisions changed')
            decisions[t]={r['protein_id']:r for r in table(path)}
        d=decisions[t].get(pid)
        if d and any(x[k]!=d[k2] for k,k2 in [('gene_ids_json','gene_ids_json'),('gene_mapping_status','status'),('representative_decision','decision')]):raise ValueError('Gene/representative join differs')
        expected={}
        for action in {h['fcs_action'] for h in hits}:
            expected[action]=sum(union([(int(h['intersection_start']),int(h['intersection_end'])) for h in hits if h['fcs_action']==action and h['seq_id']==seq]) for seq in {h['seq_id'] for h in hits})
        if json.loads(x['overlap_union_bp_by_action_json'])!=expected:raise ValueError('Protein union differs')
    markers=table(a.markers);review=table(a.mapping/'marker_overlap_review.tsv')
    if len(markers)!=len(review):raise ValueError('Marker count differs')
    for orig,x in zip(markers,review):
        if any(x[k]!=v for k,v in orig.items()):raise ValueError('Marker fields changed')
        hits=byprotein[x['taxon_id'],x['protein_id']]
        if (x['fcs_overlap_status']=='cds_overlaps_fcs_region')!=bool(hits) or x['fcs_actions']!=';'.join(sorted({h['fcs_action'] for h in hits})):raise ValueError('Marker overlap differs')
    a.output.mkdir(parents=True)
    result={'status':'passed_full_region_cds_intersection_and_protein_marker_join_audit','mapping_receipt_sha256':sha(a.mapping/'receipt.json'),'script_sha256':sha(Path(__file__)),'taxa_scanned':len(bytax),'cds_feature_rows_scanned':scan_count,'regions':len(regions),'overlap_feature_rows':len(features),'overlapping_proteins':len(proteins),'marker_links':len(markers),'interpretation':'All exact reported intervals independently tested against all annotated CDS segments on the matching sequence; full positive/negative overlap enumeration, source field identities, region/protein union counts and marker overlap joins checked. Does not establish biological contamination or resolve sequences absent from GFF.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
