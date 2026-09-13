#!/usr/bin/env python3
"""Independently read back marker CDS identities and summarize missingness by lineage."""
import argparse,csv,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path
from Bio import SeqIO
from audit_busco_gene_copies import ROOT,sha,read_table


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--index',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable readback directory')
    r=json.loads((a.index/'receipt.json').read_text())
    if r['status']!='complete_full_marker_cds_identity_index':raise ValueError('Completed full source index required')
    for name,h in r['artifacts'].items():
        if sha(a.index/name)!=h:raise ValueError('Changed index artifact')
    mapping_path=ROOT/'results/phylogeny/markers-full-v1/protein_mapping.tsv'
    if sha(mapping_path)!=r['mapping_sha256']:raise ValueError('Changed marker mapping')
    mapping=read_table(mapping_path);expected={(x['marker'],x['taxon_id']):x for x in mapping}
    rows=read_table(a.index/'marker_cds_source_index.tsv');observed={(x['marker'],x['taxon_id']):x for x in rows}
    if len(observed)!=len(rows) or len(expected)!=len(mapping) or set(expected)!=set(observed):raise ValueError('Marker identity universe differs')
    for key,row in observed.items():
        if any(row[k]!=expected[key][k] for k in expected[key]):raise ValueError('Marker provenance differs')
        ids=json.loads(row['source_record_ids_json'])
        n=len(ids);status='no_source_cds' if n==0 else 'unique_source_cds' if n==1 else 'multiple_source_cds_records'
        if row['source_status']!=status:raise ValueError('Source multiplicity/status differs')
    wanted={key:row for key,row in observed.items() if row['source_status']=='unique_source_cds'}
    seen=set();nucleotides=0
    for rec in SeqIO.parse(a.index/'marker_source_cds.fna','fasta'):
        key=tuple(rec.id.split('|'))
        if key in seen or key not in wanted:raise ValueError('Unexpected/repeated exported CDS')
        seen.add(key);dna=str(rec.seq);row=wanted[key]
        if hashlib.sha256(dna.encode()).hexdigest()!=row['cds_sequence_sha256'] or len(dna)!=int(row['cds_length']):raise ValueError('Changed exported nucleotide sequence')
        if not dna or not set(dna)<=set('ACGTRYSWKMBDHVN'):raise ValueError('Invalid exported CDS alphabet')
        nucleotides+=len(dna)
    if seen!=set(wanted):raise ValueError('Missing exported CDS')
    manifest_path=ROOT/'metadata/analysis_manifest.tsv';manifest=read_table(manifest_path);taxa={x['taxon_id']:x for x in manifest}
    local=defaultdict(list)
    for row in rows:local[row['taxon_id']].append(row)
    if set(local)!=set(taxa):raise ValueError('Full design taxon coverage differs')
    summaries=[]
    for t,m in taxa.items():
        counts=Counter(x['source_status'] for x in local[t])
        summaries.append({'taxon_id':t,'species_name':m['species_name'],'study_role':m['study_role'],'lineage':m['lineage'].split(';')[0],'marker_links':len(local[t]),'unique_source_cds':counts['unique_source_cds'],'no_source_cds':counts['no_source_cds'],'multiple_source_cds_records':counts['multiple_source_cds_records'],'alternative_product_markers':sum(x['representative_decision']=='alternative_product_retained_in_source' for x in local[t]),'unresolved_gene_markers':sum(x['gene_mapping_status']!='unique_gene' for x in local[t])})
    totals=dict(Counter(x['source_status'] for x in rows))
    if totals!=r['source_status_counts'] or len(rows)!=r['marker_links']:raise ValueError('Receipt totals differ')
    a.output.mkdir(parents=True)
    for name,data,fields in [('taxon_cds_coverage.tsv',summaries,list(summaries[0])),('marker_cds_exceptions.tsv',[x for x in rows if x['source_status']!='unique_source_cds' or x['representative_decision']=='alternative_product_retained_in_source' or x['gene_mapping_status']!='unique_gene'],list(rows[0]))]:
        with (a.output/name).open('w') as f:
            w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    result={'status':'complete_full_marker_cds_index_readback','source_index_receipt_sha256':sha(a.index/'receipt.json'),'manifest_sha256':sha(manifest_path),'marker_links':len(rows),'taxa':len(taxa),'exported_sequences':len(seen),'exported_nucleotides':nucleotides,'source_status_counts':totals,'representative_decision_counts':dict(Counter(x['representative_decision'] for x in rows)),'gene_mapping_status_counts':dict(Counter(x['gene_mapping_status'] for x in rows)),'script_sha256':sha(Path(__file__)),'interpretation':'All index artifacts, exact marker identity coverage, record multiplicity and exported sequence hashes verified independently. This is not translation verification, complete-gene validation or selection eligibility.','artifacts':{p.name:sha(p) for p in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
