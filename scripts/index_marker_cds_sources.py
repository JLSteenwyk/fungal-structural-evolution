#!/usr/bin/env python3
"""Index exact CDS source identities for every marker; preserve ambiguity and gaps."""
import argparse,csv,gzip,json,re
from collections import defaultdict,Counter
from pathlib import Path
from Bio import SeqIO
from audit_busco_gene_copies import ROOT,sha,read_table


def unique_record(records):
    if not records:return None,'no_source_cds'
    if len(records)!=1:return None,'multiple_source_cds_records'
    return records[0],'unique_source_cds'


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable source index')
    mapping_path=ROOT/'results/phylogeny/markers-full-v1/protein_mapping.tsv';mapping=read_table(mapping_path)
    if len({(r['taxon_id'],r['marker']) for r in mapping})!=len(mapping):raise ValueError('Repeated marker/taxon')
    bytax=defaultdict(list)
    for r in mapping:bytax[r['taxon_id']].append(r)
    paths={k:ROOT/'metadata'/v for k,v in [('ncbi','cds_download_receipts.json'),('external','external_cds_source_inventory.json'),('representatives','gene_representatives_receipt.json'),('proteins','qc_input_receipts.json')]}
    ncbi=json.loads(paths['ncbi'].read_text());assert ncbi['status']=='complete_ncbi_cds_acquisition'
    sources={r['taxon_id']:dict(path=r['path'],sha256=r['sha256'],kind='ncbi_unmodified_cds',receipt_sha256=sha(paths['ncbi'])) for r in ncbi['taxa']}
    ext=json.loads(paths['external'].read_text())
    projection=ROOT/'results/cds/published-outgroup-codon-projection-v2';pr=json.loads((projection/'receipt.json').read_text())
    assert pr['status']=='complete_published_outgroup_genome_codon_projection_audit'
    projected={r['taxon_id'] for r in pr['taxa']}
    for r in ext['taxa']:
        t=r['taxon_id']
        if t in projected:
            name=t+'.verified_codon_spans.fna';sources[t]=dict(path=str((projection/name).relative_to(ROOT)),sha256=pr['artifacts'][name],kind='published_genome_projected_codon_span',receipt_sha256=sha(projection/'receipt.json'))
        else:sources[t]=dict(path=r['path'],sha256=r['sha256'],kind='external_verified_cds_subset',receipt_sha256=sha(paths['external']))
    reps={r['taxon_id']:r for r in json.loads(paths['representatives'].read_text())['taxa']}
    proteins={r['taxon_id']:r for r in json.loads(paths['proteins'].read_text())}
    if set(bytax)-set(sources):raise ValueError('Missing source inventory')
    a.output.mkdir(parents=True);fasta=a.output/'marker_source_cds.fna';rows=[];summaries=[]
    with fasta.open('w') as out:
        for t,links in bytax.items():
            src=sources[t];path=ROOT/src['path']
            if sha(path)!=src['sha256']:raise ValueError('Changed CDS source')
            protein_path=ROOT/proteins[t]['input_path']
            if sha(protein_path)!=proteins[t]['sha256']:raise ValueError('Changed protein source')
            seqs={r.id:str(r.seq) for r in SeqIO.parse(protein_path,'fasta')}
            decpath=(ROOT/reps[t]['path']).with_suffix('.decisions.tsv')
            if sha(decpath)!=reps[t]['decisions_sha256']:raise ValueError('Changed representative decisions')
            decisions={r['protein_id']:r for r in read_table(decpath)}
            wanted={r['protein_id'] for r in links};hits=defaultdict(list)
            opener=gzip.open if path.suffix=='.gz' else open
            with opener(path,'rt') as f:
                for rec in SeqIO.parse(f,'fasta'):
                    if src['kind']=='ncbi_unmodified_cds':
                        ids=re.findall(r'\[protein_id=([^\]]+)\]',rec.description);pid=ids[0] if len(ids)==1 else ''
                    else:pid=rec.id
                    if pid in wanted:hits[pid].append(rec)
            local=[]
            for link in links:
                pid=link['protein_id'];aa=seqs[pid]
                import hashlib
                if hashlib.sha256(aa.encode()).hexdigest()!=link['sequence_sha256']:raise ValueError('Marker protein mismatch')
                rec,status=unique_record(hits[pid]);decision=decisions[pid]
                row=dict(link,cds_source_kind=src['kind'],source_status=status,source_record_ids_json=json.dumps([r.id for r in hits[pid]]),cds_sequence_sha256='',cds_length='',gene_ids_json=decision['gene_ids_json'],gene_mapping_status=decision['status'],representative_decision=decision['decision'])
                if rec is not None:
                    dna=str(rec.seq).upper();row.update(cds_length=len(dna),cds_sequence_sha256=hashlib.sha256(dna.encode()).hexdigest())
                    out.write('>'+link['marker']+'|'+t+'\n'+dna+'\n')
                rows.append(row);local.append(row)
            summaries.append(dict(taxon_id=t,marker_links=len(local),source_status_counts=dict(Counter(r['source_status'] for r in local)),source=src))
            print(t,len(summaries),summaries[-1]['source_status_counts'],flush=True)
    table=a.output/'marker_cds_source_index.tsv'
    with table.open('w') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    result={'status':'complete_full_marker_cds_identity_index','marker_links':len(rows),'taxa':summaries,'source_status_counts':dict(Counter(r['source_status'] for r in rows)),'mapping_sha256':sha(mapping_path),'source_receipt_hashes':{k:sha(v) for k,v in paths.items()},'projection_receipt_sha256':sha(projection/'receipt.json'),'script_sha256':sha(Path(__file__)),'artifacts':{p.name:sha(p) for p in [fasta,table]},'interpretation':'Exact protein identity and source CDS indexing only. NCBI sequences are unmodified and not yet translation-qualified here; external sets inherit explicit extraction/partial-boundary policies. No selection eligibility claim. Representative alternatives and unresolved gene mappings remain flagged; no silent representative replacement. Join completed translation and boundary audits before codon alignment.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
