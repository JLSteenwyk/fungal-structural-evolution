#!/usr/bin/env python3
"""Strict full-design CDS/protein translation audit with annotation code provenance."""
import argparse,csv,gzip,json,re,fcntl
from collections import Counter,defaultdict
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from map_proteins_to_genes import attributes,ROOT
from download_gene_annotations import digest


def compare(dna,protein,table):
    if len(dna)%3:return 'non_triplet_length',False
    aa=str(Seq(dna).translate(table=table));stop=aa.endswith('*')
    if stop:aa=aa[:-1]
    return ('exact_translation' if aa==protein else 'translation_mismatch'),stop


def code_choice(values):
    if len(values)>1:return None,'conflicting_annotation_codes'
    if not values:return 1,'table_1_assumption_no_explicit_code'
    return int(next(iter(values))),'explicit_gff_transl_table'


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    a.output.mkdir(parents=True,exist_ok=True)
    lock=(a.output/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    paths={k:ROOT/'metadata'/v for k,v in [('cds','cds_download_receipts.json'),('gff','annotation_download_receipts.json'),('protein','qc_input_receipts.json')]}
    config={'script_sha256':digest(Path(__file__)),'input_receipts':{k:digest(p) for k,p in paths.items()},'policy':'Strict unmodified CDS; annotated translation code or explicitly flagged table-1 assumption; no phase trimming, initiation-residue repair, recoding, frame search or selection eligibility claim.'}
    cp=a.output/'config.json'
    if cp.exists():
        if json.loads(cp.read_text())!=config:raise ValueError('Changed audit configuration')
    else:cp.write_text(json.dumps(config,indent=2)+'\n')
    sources=json.loads(paths['cds'].read_text());assert sources['status']=='complete_ncbi_cds_acquisition'
    annotations={r['taxon_id']:r for r in json.loads(paths['gff'].read_text())}
    inputs={r['taxon_id']:r for r in json.loads(paths['protein'].read_text())}
    summaries=[]
    for src in sources['taxa']:
        t=src['taxon_id'];ann=annotations[t];inp=inputs[t]
        assert ann['assembly_accession']==src['assembly_accession']
        source_paths={'cds':ROOT/src['path'],'gff':ROOT/ann['path'],'protein':ROOT/inp['input_path']}
        expected={'cds':src['sha256'],'gff':ann['sha256'],'protein':inp['sha256']}
        for k,p in source_paths.items():
            if digest(p)!=expected[k]:raise ValueError('Changed '+k+' source '+t)
        rp=a.output/(t+'.receipt.json')
        if rp.exists():
            r=json.loads(rp.read_text())
            if r['config_sha256']!=digest(cp) or r['source_hashes']!=expected:raise ValueError('Changed taxon audit')
            for name,h in r['artifacts'].items():
                if digest(a.output/name)!=h:raise ValueError('Changed taxon artifact')
            summaries.append(r);continue
        codes=defaultdict(set);flags=defaultdict(set);annotated=set()
        with gzip.open(source_paths['gff'],'rt') as f:
            for line in f:
                if line.startswith('##FASTA'):break
                if line.startswith('#') or not line.strip():continue
                fields=line.rstrip('\n').split('\t')
                if fields[2]!='CDS':continue
                attrs=attributes(fields[8])
                for pid in attrs.get('protein_id',[]):
                    annotated.add(pid);codes[pid].update(attrs.get('transl_table',[]))
                    for key in ['exception','transl_except','partial','pseudo','pseudogene']:
                        for value in attrs.get(key,[]):flags[pid].add(key+'='+value)
        proteins={r.id:str(r.seq) for r in SeqIO.parse(source_paths['protein'],'fasta')}
        records=[];multiplicity=Counter()
        with gzip.open(source_paths['cds'],'rt') as f:
            for rec in SeqIO.parse(f,'fasta'):
                fields=re.findall(r'\[protein_id=([^\]]+)\]',rec.description)
                pid=fields[0] if len(fields)==1 else ''
                records.append((rec,pid));multiplicity[pid]+=1
        table=a.output/(t+'.audit.tsv');tmp=table.with_suffix('.partial')
        counts=Counter();code_counts=Counter();seen=set()
        columns=['taxon_id','cds_id','protein_id','status','translation_table','code_source','terminal_stop','gff_flags','cds_header']
        with tmp.open('w') as f:
            w=csv.DictWriter(f,columns,delimiter='\t',lineterminator='\n');w.writeheader()
            for rec,pid in records:
                row=dict(taxon_id=t,cds_id=rec.id,protein_id=pid,status='',translation_table='',code_source='',terminal_stop='',gff_flags=';'.join(sorted(flags[pid])),cds_header=rec.description)
                if not pid:status='missing_or_ambiguous_protein_id'
                elif pid not in proteins:status='protein_id_absent_from_normalized_proteome'
                elif multiplicity[pid]!=1:status='multiple_cds_records_for_protein'
                elif pid not in annotated:status='protein_id_absent_from_gff_cds'
                else:
                    code,provenance=code_choice(codes[pid]);row.update(translation_table=code,code_source=provenance)
                    code_counts[provenance+':'+str(code)]+=1
                    if code is None:status='conflicting_annotation_codes'
                    elif any(x.startswith(('exception=','transl_except=','pseudo=','pseudogene=')) for x in flags[pid]):status='annotation_exception_requires_review'
                    else:
                        try:status,row['terminal_stop']=compare(str(rec.seq).upper(),proteins[pid],code)
                        except (ValueError,KeyError) as e:status='translation_error:'+str(e)
                if pid in proteins:seen.add(pid)
                row['status']=status;counts[status]+=1;w.writerow(row)
            for pid in sorted(set(proteins)-seen):
                w.writerow(dict(taxon_id=t,protein_id=pid,status='protein_without_cds_record'))
        tmp.replace(table)
        r={'taxon_id':t,'status':'complete_strict_cds_translation_audit','config_sha256':digest(cp),'source_hashes':expected,'cds_records':len(records),'proteins':len(proteins),'proteins_without_cds':len(set(proteins)-seen),'cds_status_counts':dict(counts),'code_source_counts':dict(code_counts),'artifacts':{table.name:digest(table)}}
        rp.write_text(json.dumps(r,indent=2)+'\n');summaries.append(r)
        print(t,len(summaries),dict(counts),flush=True)
    result={'status':'complete_full_ncbi_strict_cds_translation_audit','config_sha256':digest(cp),'taxa':summaries,'taxa_count':len(summaries),'interpretation':config['policy']}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
