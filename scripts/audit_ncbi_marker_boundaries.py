#!/usr/bin/env python3
"""Audit complete NCBI marker set for code provenance and partial CDS boundaries."""
import argparse,csv,gzip,json
from collections import Counter,defaultdict
from pathlib import Path
from Bio import SeqIO
from map_proteins_to_genes import ROOT,attributes
from audit_busco_gene_copies import sha,read_table
from audit_ncbi_cds_translation import compare


def describe(parts,regions):
    if not parts:return {'annotation_status':'missing_cds_features'}
    if len({(p['seqid'],p['strand']) for p in parts})!=1:return {'annotation_status':'multiple_sequences_or_strands'}
    strand=parts[0]['strand']
    if strand not in ['+','-']:return {'annotation_status':'unknown_strand'}
    ordered=sorted(parts,key=lambda x:x['start'],reverse=strand=='-')
    codes=set();provenance=set();flags=set();length=0;issues=[]
    initial=ordered[0]['phase']
    if initial not in ['0','1','2']:issues.append('missing_initial_phase')
    for n,p in enumerate(ordered):
        attr=p['attributes'];explicit=set(attr.get('transl_table',[]));regional=set(regions.get(p['seqid'],[]))
        if explicit:
            codes.update(explicit);provenance.add('explicit_cds')
            if regional and regional!=explicit:issues.append('cds_region_code_disagreement')
        elif regional:codes.update(regional);provenance.add('source_region')
        else:codes.add('1');provenance.add('ncbi_documented_default')
        for key in ['partial','start_range','end_range','exception','transl_except','pseudo','pseudogene']:
            if key in attr:flags.add(key)
        if p['end']<p['start'] or p['start']<1:issues.append('invalid_interval')
        if initial in ['0','1','2']:
            expected=int(initial) if n==0 else (3-(length-int(initial))%3)%3
            if p['phase']!=str(expected):issues.append('inconsistent_internal_phase')
        length+=p['end']-p['start']+1
    intervals=sorted((x['start'],x['end']) for x in parts)
    if any(x[1]>=y[0] for x,y in zip(intervals,intervals[1:])):issues.append('overlapping_cds_intervals')
    if len(codes)!=1:issues.append('conflicting_translation_codes')
    firstkey='start_range' if strand=='+' else 'end_range';lastkey='end_range' if strand=='+' else 'start_range'
    five=firstkey in ordered[0]['attributes'];three=lastkey in ordered[-1]['attributes']
    internal=any(('start_range' in p['attributes'] and not ((strand=='+' and i==0) or (strand=='-' and i==len(ordered)-1))) or ('end_range' in p['attributes'] and not ((strand=='-' and i==0) or (strand=='+' and i==len(ordered)-1))) for i,p in enumerate(ordered))
    return {'annotation_status':'annotation_flags' if issues else 'consistent_ordered_cds_features','translation_table':next(iter(codes)) if len(codes)==1 else '', 'code_provenance':';'.join(sorted(provenance)),'annotation_issues':';'.join(sorted(set(issues))),'initial_phase':initial,'genomic_cds_span_bases':length,'cds_parts':len(parts),'five_prime_partial':five,'three_prime_partial':three,'internal_partial_boundary':internal,'annotation_flags':';'.join(sorted(flags)),'ordered_segments_json':json.dumps(ordered)}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable boundary audit')
    base=ROOT/'results/cds/marker-source-index-v1';receipt=json.loads((base/'receipt.json').read_text())
    assert receipt['status']=='complete_full_marker_cds_identity_index'
    for name,h in receipt['artifacts'].items():
        if sha(base/name)!=h:raise ValueError('Changed marker index')
    rows=read_table(base/'marker_cds_source_index.tsv');bytax=defaultdict(list)
    for r in rows:
        if r['cds_source_kind']=='ncbi_unmodified_cds':bytax[r['taxon_id']].append(r)
    if len(bytax)!=519:raise ValueError('Expected full 519 NCBI taxa')
    cds={r.id:str(r.seq) for r in SeqIO.parse(base/'marker_source_cds.fna','fasta')}
    ap=ROOT/'metadata/annotation_download_receipts.json';annotations={r['taxon_id']:r for r in json.loads(ap.read_text())}
    pp=ROOT/'metadata/qc_input_receipts.json';proteins={r['taxon_id']:r for r in json.loads(pp.read_text())}
    output=[];sources=[];summaries=[];a.output.mkdir(parents=True)
    for t,links in bytax.items():
        ann=annotations[t];path=ROOT/ann['path'];inp=proteins[t];protein_path=ROOT/inp['input_path']
        if sha(path)!=ann['sha256'] or sha(protein_path)!=inp['sha256']:raise ValueError('Changed GFF/protein source')
        wanted={r['protein_id'] for r in links};parts=defaultdict(list);regions=defaultdict(set)
        with gzip.open(path,'rt') as f:
            for line in f:
                if line.startswith('##FASTA'):break
                if line.startswith('#') or not line.strip():continue
                x=line.rstrip('\n').split('\t')
                if x[2] not in ['CDS','region']:continue
                attr=attributes(x[8])
                if x[2]=='region':regions[x[0]].update(attr.get('transl_table',[]));continue
                for pid in set(attr.get('protein_id',[]))&wanted:
                    parts[pid].append({'seqid':x[0],'start':int(x[3]),'end':int(x[4]),'strand':x[6],'phase':x[7],'attributes':attr})
        seqs={r.id:str(r.seq) for r in SeqIO.parse(protein_path,'fasta') if r.id in wanted};local=[]
        for link in links:
            pid=link['protein_id'];row={'taxon_id':t,'marker':link['marker'],'protein_id':pid,**describe(parts[pid],regions),'translation_status':'unresolved_annotation_code','terminal_stop':''}
            code=row.get('translation_table','');dna=cds.get(link['marker']+'|'+t)
            if dna is None:row['translation_status']='no_unique_indexed_cds'
            elif code:
                try:row['translation_status'],row['terminal_stop']=compare(dna,seqs[pid],int(code))
                except (ValueError,KeyError) as e:row['translation_status']='translation_error:'+str(e)
            row['indexed_cds_length']=len(dna) if dna else ''
            local.append(row);output.append(row)
        summaries.append({'taxon_id':t,'marker_links':len(local),'translation_status_counts':dict(Counter(r['translation_status'] for r in local)),'annotation_status_counts':dict(Counter(r['annotation_status'] for r in local))})
        sources.append({'taxon_id':t,'gff_sha256':ann['sha256'],'protein_sha256':inp['sha256']});print(t,len(summaries),summaries[-1]['translation_status_counts'],flush=True)
    table=a.output/'marker_annotation_translation.tsv';fields=list(dict.fromkeys(k for r in output for k in r))
    with table.open('w') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(output)
    result={'status':'complete_full_ncbi_marker_boundary_code_audit','marker_links':len(output),'taxa':summaries,'sources':sources,'source_index_receipt_sha256':sha(base/'receipt.json'),'annotation_inventory_sha256':sha(ap),'protein_inventory_sha256':sha(pp),'script_sha256':sha(Path(__file__)),'translation_helper_sha256':sha(ROOT/'scripts/audit_ncbi_cds_translation.py'),'translation_status_counts':dict(Counter(r['translation_status'] for r in output)),'code_provenance_counts':dict(Counter(r.get('code_provenance','') for r in output)),'interpretation':'NCBI CDS code, source-region code, or documented table-1 default, with disagreements flagged. Strict unmodified indexed DNA translation; no frame repair. Strand-aware terminal and internal partial boundaries retained. This is annotation/translation review, not genome sequence reconstruction or codon-selection eligibility.','policy_sources':['https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/','https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/taxonomy-processing/genetic-codes/'],'artifacts':{table.name:sha(table)}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
