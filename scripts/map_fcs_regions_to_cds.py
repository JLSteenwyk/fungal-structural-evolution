#!/usr/bin/env python3
"""Intersect FCS intervals with exact-version annotated CDS segments and marker IDs."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
from map_proteins_to_genes import attributes


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def table(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def union_length(intervals):
    total=0;end=-1
    for start,stop in sorted(intervals):
        total+=max(0,stop-max(end,start)+1);end=max(end,stop+1)
    return total


def overlap(a,b,c,d):return max(0,min(b,d)-max(a,c)+1)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['fcs','annotations','representatives','markers','marker-index','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use immutable output')
    fr=json.loads((a.fcs/'receipt.json').read_text())
    for name,digest in fr['artifacts'].items():
        if sha(a.fcs/name)!=digest:raise ValueError('Changed audited FCS artifact')
    index=json.loads((a.marker_index/'receipt.json').read_text())
    if sha(a.markers)!=index['mapping_sha256']:raise ValueError('Marker mapping differs from full CDS identity index')
    sources={x['taxon_id']:x for x in json.loads(a.annotations.read_text())};reps={x['taxon_id']:x for x in json.loads(a.representatives.read_text())['taxa']}
    markers=table(a.markers);byprotein=defaultdict(list)
    for r in markers:byprotein[r['taxon_id'],r['protein_id']].append(r['marker'])
    summaries={r['taxon_id']:r for r in table(a.fcs/'taxon_summary.tsv')};regions=table(a.fcs/'flagged_regions.tsv');bytax=defaultdict(list)
    for i,r in enumerate(regions):
        r['region_id']='FCS'+str(i+1).zfill(6)
        for key in ['start_pos','end_pos','seq_len']:r[key]=int(r[key])
        bytax[r['taxon_id']].append(r)
    a.output.mkdir(parents=True);features=[];region_summary=[];protein_rows=[];taxon_rows=[];proofs=[];hits_byprotein=defaultdict(list);unresolved_taxa=set()
    for taxon,local in sorted(bytax.items()):
        source=sources[taxon];path=Path(source['path']);acc=summaries[taxon]['assembly_accession']
        if source['assembly_accession']!=acc or sha(path)!=source['sha256']:raise ValueError('Annotation assembly/hash mismatch')
        decision=Path(reps[taxon]['path']).with_suffix('.decisions.tsv')
        if sha(decision)!=reps[taxon]['decisions_sha256']:raise ValueError('Representative decision changed')
        decisions={r['protein_id']:r for r in table(decision)}
        byseq=defaultdict(list)
        for r in local:byseq[r['seq_id']].append(r)
        seqids=set();bounds={};hits=defaultdict(list);cds_total=0;missing_pid=0;header_accessions=[]
        with gzip.open(path,'rt') as f:
            for lineno,line in enumerate(f,1):
                if line.startswith('##FASTA'):break
                if line.startswith('#!genome-build-accession '):header_accessions.append(line.strip().split(':')[-1])
                if line.startswith('##sequence-region '):
                    _,seq,start,end=line.split();seqids.add(seq);bounds[seq]=(int(start),int(end))
                if line.startswith('#') or not line.strip():continue
                fields=line.rstrip('\n').split('\t')
                if len(fields)!=9:raise ValueError('Malformed GFF')
                seqids.add(fields[0])
                if fields[2]!='CDS':continue
                cds_total+=1
                if fields[0] not in byseq:continue
                start,end=int(fields[3]),int(fields[4]);attrs=attributes(fields[8]);pids=attrs.get('protein_id',[])
                if not 1<=start<=end:raise ValueError('Invalid GFF CDS coordinates')
                for region in byseq[fields[0]]:
                    n=overlap(start,end,region['start_pos'],region['end_pos'])
                    if not n:continue
                    hit={'taxon_id':taxon,'assembly_accession':acc,'region_id':region['region_id'],'seq_id':fields[0],'gff_line':lineno,'cds_start':start,'cds_end':end,'strand':fields[6],'phase':fields[7],'cds_ids_json':json.dumps(attrs.get('ID',[])),'protein_ids_json':json.dumps(pids),'parent_ids_json':json.dumps(attrs.get('Parent',[])),'fcs_action':region['action'],'intersection_start':max(start,region['start_pos']),'intersection_end':min(end,region['end_pos']),'overlap_bp':n}
                    features.append(hit);hits[region['region_id']].append(hit)
                    if not pids:missing_pid+=1
                    for pid in pids:hits_byprotein[taxon,pid].append(hit)
        if header_accessions!=[acc]:raise ValueError('GFF header does not establish exact assembly version')
        unknown=0
        for region in local:
            matched=region['seq_id'] in seqids
            if not matched:unknown+=1;unresolved_taxa.add(taxon)
            b=bounds.get(region['seq_id'])
            selected=hits[region['region_id']]
            # Independent direct overlap count is retained per CDS feature, not intron-spanning gene interval.
            covered=union_length([(x['intersection_start'],x['intersection_end']) for x in selected])
            if covered>region['end_pos']-region['start_pos']+1:raise ValueError('Region union exceeds region')
            region_summary.append(dict(region,gff_sequence_status='represented_in_gff' if matched else 'not_represented_in_gff',gff_sequence_bounds_json=json.dumps(b),reported_sequence_length_matches_gff_bounds='' if b is None else b==(1,region['seq_len']),overlapping_cds_feature_rows=len(selected),overlapping_cds_union_bp=covered,overlapping_protein_ids_json=json.dumps(sorted({pid for x in selected for pid in json.loads(x['protein_ids_json'])}))))
        local_proteins=sorted(pid for t,pid in hits_byprotein if t==taxon)
        for pid in local_proteins:
            selected=hits_byprotein[taxon,pid];d=decisions.get(pid);actions=sorted({r['fcs_action'] for r in selected})
            psummary={'taxon_id':taxon,'assembly_accession':acc,'protein_id':pid,'gene_ids_json':d['gene_ids_json'] if d else '[]','gene_mapping_status':d['status'] if d else 'protein_absent_from_representative_decisions','representative_decision':d['decision'] if d else 'unresolved','fcs_actions':';'.join(actions),'fcs_region_ids':';'.join(sorted({r['region_id'] for r in selected})),'markers':';'.join(sorted(byprotein[taxon,pid]))}
            action_lengths={action:sum(union_length([(x['intersection_start'],x['intersection_end']) for x in selected if x['fcs_action']==action and x['seq_id']==seq]) for seq in {x['seq_id'] for x in selected}) for action in actions}
            psummary['overlap_union_bp_by_action_json']=json.dumps(action_lengths,sort_keys=True);protein_rows.append(psummary)
        taxon_rows.append({'taxon_id':taxon,'assembly_accession':acc,'fcs_regions':len(local),'gff_cds_feature_rows_scanned':cds_total,'regions_without_gff_sequence':unknown,'overlapping_proteins':len(local_proteins),'overlap_feature_rows_without_protein_id':missing_pid})
        proofs.append({'taxon_id':taxon,'annotation_path':str(path),'annotation_sha256':sha(path),'representative_decision_sha256':sha(decision)})
        if len(proofs)%10==0:print('Mapped',len(proofs),'of',len(bytax),'flagged taxa',flush=True)
    marker_rows=[]
    for r in markers:
        t=r['taxon_id'];hits=hits_byprotein[t,r['protein_id']];s=summaries[t];validated=s['validation_status'] in ['publisher_verified_report','repeat_identical_https_report_without_publisher_md5']
        status='cds_overlaps_fcs_region' if hits else 'fcs_report_unresolved_or_external' if not validated else 'no_regions_reported_by_fcs' if int(s['report_rows'])==0 else 'no_direct_cds_overlap_but_fcs_sequence_mapping_incomplete' if t in unresolved_taxa else 'no_annotated_cds_overlap_with_reported_regions'
        marker_rows.append(dict(r,fcs_overlap_status=status,fcs_actions=';'.join(sorted({x['fcs_action'] for x in hits})),fcs_region_ids=';'.join(sorted({x['region_id'] for x in hits}))))
    artifacts={}
    fallback=['taxon_id','assembly_accession','protein_id']
    for name,data in [('cds_feature_overlaps.tsv',features),('region_annotation_review.tsv',region_summary),('protein_overlap_review.tsv',protein_rows),('taxon_mapping_summary.tsv',taxon_rows),('marker_overlap_review.tsv',marker_rows)]:
        with (a.output/name).open('w',newline='') as f:
            w=csv.DictWriter(f,list(data[0]) if data else fallback,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
        artifacts[name]=sha(a.output/name)
    result={'status':'complete_fcs_cds_marker_overlap_mapping_pending_independent_audit','source_receipts':{'fcs':sha(a.fcs/'receipt.json'),'annotations':sha(a.annotations),'representatives':sha(a.representatives),'marker_index':sha(a.marker_index/'receipt.json')},'marker_mapping_sha256':sha(a.markers),'script_sha256':sha(Path(__file__)),'gff_attribute_parser_sha256':sha(Path(__file__).with_name('map_proteins_to_genes.py')),'flagged_taxa_scanned':len(proofs),'regions':len(regions),'cds_feature_overlap_rows':len(features),'overlapping_proteins':len(protein_rows),'overlapping_selected_representatives':sum(r['representative_decision'] in ['longest_per_gene_lexical_tiebreak','retained_unresolved_gene'] for r in protein_rows),'marker_links':len(marker_rows),'marker_overlap_status_counts':dict(Counter(r['fcs_overlap_status'] for r in marker_rows)),'annotation_proofs':proofs,'artifacts':artifacts,'interpretation':'Exact CDS-segment intersections with FCS report intervals. Overlap is a review flag, not confirmed foreign origin. No intron-only overlap is counted as coding overlap; duplicate/overlapping CDS features are unioned for per-protein action summaries. Unrepresented GFF sequences and missing reports stay unresolved. No deletion, structural residue localization or recomputed phylogeny.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['annotation_proofs','artifacts']},indent=2))


if __name__=='__main__':main()
