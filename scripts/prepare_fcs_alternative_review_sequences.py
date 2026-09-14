#!/usr/bin/env python3
"""Materialize score-passing alternative proteins for review, without replacing markers."""
import argparse,csv,json,hashlib
from collections import Counter,defaultdict
from pathlib import Path
from Bio import SeqIO


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def table(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['review','output']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use new immutable output')
    r=json.loads((a.review/'receipt.json').read_text());cp=a.review/'score_passing_unflagged_alternatives.tsv'
    if sha(cp)!=r['artifacts'][cp.name]:raise ValueError('Candidate review changed')
    candidates=table(cp);source_path=Path('metadata/qc_input_receipts.json');sources={x['taxon_id']:x for x in json.loads(source_path.read_text())}
    if sha(source_path)!=r['source_hashes']['proteins']:raise ValueError('Protein receipts changed')
    cache={};cds={};proofs=[];outrows=[];a.output.mkdir(parents=True)
    with (a.output/'candidate_proteins.faa').open('w') as fasta:
        for row in candidates:
            t=row['taxon_id'];pid=row['protein_id']
            if t not in cache:
                path=Path(sources[t]['input_path'])
                if sha(path)!=sources[t]['sha256']:raise ValueError('Protein source changed')
                cache[t]={x.id:str(x.seq) for x in SeqIO.parse(path,'fasta')}
                audit=Path('results/cds/ncbi-strict-translation-v1')/(t+'.receipt.json');ar=json.loads(audit.read_text());ap=audit.with_name(t+'.audit.tsv')
                if sha(ap)!=ar['artifacts'][ap.name] or ar['source_hashes']['protein']!=sources[t]['sha256']:raise ValueError('CDS translation audit provenance mismatch')
                cds[t]=defaultdict(list)
                for x in table(ap):cds[t][x['protein_id']].append(x)
                proofs.append({'taxon_id':t,'protein_sha256':sha(path),'cds_audit_receipt_sha256':sha(audit),'cds_audit_table_sha256':sha(ap)})
            seq=cache[t][pid]
            if hashlib.sha256(seq.encode()).hexdigest()!=row['sequence_sha256']:raise ValueError('Candidate sequence differs')
            ident='|'.join([t,row['marker'],pid]);fasta.write('>'+ident+'\n'+seq+'\n');cr=cds[t].get(pid,[])
            outrows.append(dict(row,review_sequence_id=ident,cds_audit_records=len(cr),cds_audit_statuses=';'.join(sorted({x['status'] for x in cr})),translation_tables=';'.join(sorted({x['translation_table'] for x in cr if x['translation_table']})),review_disposition='candidate_only_phylogenetic_copy_and_alignment_review_required'))
    with (a.output/'candidate_review.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,list(outrows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(outrows)
    reread={x.id:str(x.seq) for x in SeqIO.parse(a.output/'candidate_proteins.faa','fasta')}
    if len(reread)!=len(outrows) or any(hashlib.sha256(reread[x['review_sequence_id']].encode()).hexdigest()!=x['sequence_sha256'] for x in outrows):raise ValueError('Full candidate FASTA readback differs')
    result={'status':'complete_alternative_review_sequence_materialization','source_review_receipt_sha256':sha(a.review/'receipt.json'),'script_sha256':sha(Path(__file__)),'candidate_marker_protein_records':len(outrows),'unique_proteins':len({(x['taxon_id'],x['protein_id']) for x in outrows}),'marker_taxon_observations':len({(x['taxon_id'],x['marker']) for x in outrows}),'cds_status_counts':dict(Counter(x['cds_audit_statuses'] for x in outrows)),'source_proofs':proofs,'artifacts':{p.name:sha(p) for p in a.output.iterdir()},'interpretation':'Exact original alternative-protein sequences and existing CDS translation evidence only. Raw HMM score passage and lack of recorded FCS overlap do not establish orthology or cleanliness. No marker replacements, new predictions or inferred evolutionary results.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['source_proofs','artifacts']},indent=2))


if __name__=='__main__':main()
