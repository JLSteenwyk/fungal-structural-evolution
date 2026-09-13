#!/usr/bin/env python3
"""Combine complete focused/marker domain hits and restore every gene identity."""
import argparse,csv,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path
from Bio import SeqIO
from audit_busco_gene_copies import ROOT,sha,read_table


def checked(folder,status=None):
    r=json.loads((folder/'receipt.json').read_text())
    if status and r.get('status')!=status:raise ValueError('Incomplete prerequisite')
    for name,h in r['artifacts'].items():
        if sha(folder/name)!=h:raise ValueError('Changed prerequisite '+name)
    return r


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable candidate collection')
    focus=ROOT/'results/domains/tfiib-family-search-v1';fr=checked(focus,'complete_focused_tfiib_domain_search');fc=json.loads((focus/'config.json').read_text())
    full=ROOT/'data/domains/full-inputs-v1';ir=checked(full,'complete_search_input_preparation')
    marker_inputs=ROOT/'data/domains/marker-inputs-v1';mi=checked(marker_inputs)
    marker_hits=ROOT/'results/domains/marker-annotations-v1';mr=checked(marker_hits)
    if fc['full_input_receipt_sha256']!=sha(full/'receipt.json') or fc['marker_annotation_receipt_sha256']!=sha(marker_hits/'receipt.json') or ir['marker_input_receipt_sha256']!=sha(marker_inputs/'receipt.json'):raise ValueError('Incompatible candidate sources')
    wanted=set(fr['profiles']);hits=[]
    with (focus/'hits.domtblout').open() as f:
        for line in f:
            if not line.strip() or line.startswith('#'):continue
            x=line.split(maxsplit=22)
            hits.append({'sequence_id':x[0],'protein_length':int(x[2]),'pfam_accession':x[4],'pfam_name':x[3],'hmm_length':int(x[5]),'hmm_start':int(x[15]),'hmm_end':int(x[16]),'alignment_start':int(x[17]),'alignment_end':int(x[18]),'envelope_start':int(x[19]),'envelope_end':int(x[20]),'domain_score':float(x[13]),'source':'focused_additional'})
    with (marker_hits/'raw_annotated_hits.tsv').open() as f:
        for x in csv.DictReader(f,delimiter='\t'):
            if x['pfam_accession'] not in wanted:continue
            row={k:x[k] for k in hits[0] if k!='source'}
            for k in ['protein_length','hmm_length','hmm_start','hmm_end','alignment_start','alignment_end','envelope_start','envelope_end']:row[k]=int(row[k])
            row['domain_score']=float(row['domain_score']);row['source']='reused_marker';hits.append(row)
    hit_by_sequence=defaultdict(list)
    for h in hits:
        if h['pfam_accession'] not in wanted:raise ValueError('Unexpected profile')
        hit_by_sequence[h['sequence_id']].append(h)
    candidates=[];seen=set()
    with (full/'protein_links.tsv').open() as f:
        for x in csv.DictReader(f,delimiter='\t'):
            if x['sequence_id'] not in hit_by_sequence:continue
            key=(x['taxon_id'],x['protein_id'])
            if key in seen:raise ValueError('Repeated representative protein identity')
            seen.add(key);candidates.append(x)
    used={x['sequence_id'] for x in candidates};seqs={}
    for path in [full/'additional_sequences.faa',marker_inputs/'sequences.faa']:
        for rec in SeqIO.parse(path,'fasta'):
            if rec.id not in used:continue
            seq=str(rec.seq)
            if rec.id in seqs or rec.id!='S'+hashlib.sha256(seq.encode()).hexdigest():raise ValueError('Repeated or changed sequence identity')
            if any(x['protein_length']!=len(seq) or not 1<=x['alignment_start']<=x['alignment_end']<=len(seq) for x in hit_by_sequence[rec.id]):raise ValueError('Hit coordinates/length differ')
            seqs[rec.id]=seq
    if set(seqs)!=used:raise ValueError('Missing candidate sequences')
    rep_path=ROOT/'metadata/gene_representatives_receipt.json';rr=json.loads(rep_path.read_text())
    if sha(rep_path)!=ir['source_receipt_sha256']:raise ValueError('Changed representative source')
    bytax=defaultdict(list)
    for x in candidates:bytax[x['taxon_id']].append(x)
    decisions={}
    for x in rr['taxa']:
        t=x['taxon_id']
        if t not in bytax:continue
        path=(ROOT/x['path']).with_suffix('.decisions.tsv')
        if sha(path)!=x['decisions_sha256']:raise ValueError('Changed gene decisions')
        ids={z['protein_id'] for z in bytax[t]}
        for d in read_table(path):
            if d['protein_id'] in ids:decisions[t,d['protein_id']]=d
    marker_map=read_table(ROOT/'results/phylogeny/markers-full-v1/protein_mapping.tsv');busco=defaultdict(list)
    for x in marker_map:busco[x['taxon_id'],x['protein_id']].append(x['marker'])
    a.output.mkdir(parents=True);rows=[]
    with (a.output/'candidate_proteins.faa').open('w') as f:
        for n,x in enumerate(sorted(candidates,key=lambda x:(x['taxon_id'],x['protein_id'])),1):
            t=x['taxon_id'];pid=x['protein_id'];sid=x['sequence_id'];d=decisions[t,pid];entry='G'+str(n).zfill(7)
            if d['decision']=='alternative_product_retained_in_source':raise ValueError('Nonrepresentative protein entered candidate set')
            labels=Counter(z['pfam_accession'] for z in hit_by_sequence[sid])
            row={'entry_id':entry,**x,'protein_length':len(seqs[sid]),'gene_ids_json':d['gene_ids_json'],'gene_mapping_status':d['status'],'representative_decision':d['decision'],'TFIIB_hits':labels['PF00382.25'],'BRF1_hits':labels['PF07741.19'],'Zn_Ribbon_TF_hits':labels['PF08271.18'],'selected_busco_markers_json':json.dumps(busco[t,pid])};rows.append(row)
            f.write('>'+entry+'\n'+seqs[sid]+'\n')
    mp=ROOT/'metadata/analysis_manifest.tsv';manifest=read_table(mp);summaries=[]
    for m in manifest:
        local=[x for x in rows if x['taxon_id']==m['taxon_id']]
        summaries.append({'taxon_id':m['taxon_id'],'species_name':m['species_name'],'study_role':m['study_role'],'lineage':m['lineage'],'candidate_proteins':len(local),'TFIIB_hit_proteins':sum(x['TFIIB_hits']>0 for x in local),'BRF1_hit_proteins':sum(x['BRF1_hits']>0 for x in local),'Zn_Ribbon_TF_hit_proteins':sum(x['Zn_Ribbon_TF_hits']>0 for x in local),'unresolved_gene_candidates':sum(x['gene_mapping_status']!='unique_gene' for x in local)})
    for name,data in [('candidate_protein_mapping.tsv',rows),('domain_hits.tsv',hits),('taxon_candidate_coverage.tsv',summaries)]:
        with (a.output/name).open('w') as f:
            w=csv.DictWriter(f,list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    unused=sorted(set(hit_by_sequence)-used);(a.output/'hit_sequences_without_representative_links.json').write_text(json.dumps(unused,indent=2)+'\n')
    result={'status':'complete_full_tfiib_domain_candidate_collection','taxa_screened':len(manifest),'taxa_with_candidates':len(bytax),'candidate_proteins':len(rows),'unique_candidate_sequences':len(seqs),'raw_domain_hit_rows':len(hits),'source_hit_rows':dict(Counter(x['source'] for x in hits)),'hit_sequences_without_representative_links':len(unused),'gene_mapping_status_counts':dict(Counter(x['gene_mapping_status'] for x in rows)),'focus_receipt_sha256':sha(focus/'receipt.json'),'full_input_receipt_sha256':sha(full/'receipt.json'),'marker_input_receipt_sha256':sha(marker_inputs/'receipt.json'),'marker_annotation_receipt_sha256':sha(marker_hits/'receipt.json'),'representative_receipt_sha256':sha(rep_path),'marker_mapping_sha256':sha(ROOT/'results/phylogeny/markers-full-v1/protein_mapping.tsv'),'manifest_sha256':sha(mp),'script_sha256':sha(Path(__file__)),'interpretation':'Union of three selected Pfam domain-hit candidates across complete representative-protein inputs, with marker-hit reuse. Exact-sequence search deduplication is expanded back to every gene/taxon entry. This union is not a resolved homologous full-length family: shared zinc-ribbon or TFIIB domains can occur in different architectures. Gene/domain correspondence and phylogenetic reconciliation remain required. Zero candidates means no selected GA hit, not demonstrated biological absence. Nonrepresentative marker-only hits remain listed separately.','artifacts':{p.name:sha(p) for p in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='artifacts'},indent=2))

if __name__=='__main__':main()
