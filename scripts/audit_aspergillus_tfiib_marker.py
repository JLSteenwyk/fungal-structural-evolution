#!/usr/bin/env python3
"""Audit domain-hit heterogeneity in the leading Aspergillus divergence case."""
import argparse,csv,hashlib,json,statistics
from collections import defaultdict,Counter
from pathlib import Path
from Bio import SeqIO
from audit_busco_gene_copies import ROOT,sha,read_table


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable case audit')
    marker='4986044at2759';genus='Aspergillus';base=ROOT/'results/cds/outlier-realignment-v1';r=json.loads((base/'receipt.json').read_text())
    if r['status']!='complete_targeted_codon_alignment_review':raise ValueError('Completed realignment required')
    for name,h in r['artifacts'].items():
        if sha(base/name)!=h:raise ValueError('Changed realignment artifact')
    seqs={x.id:str(x.seq) for x in SeqIO.parse(base/(marker+'.'+genus+'.input.faa'),'fasta')}
    identities={'S'+hashlib.sha256(s.encode()).hexdigest():t for t,s in seqs.items()}
    if len(identities)!=len(seqs):raise ValueError('Repeated protein sequence requires explicit many-to-one handling')
    ann=ROOT/'results/domains/marker-annotations-v1';ap=ROOT/'metadata/marker_domain_annotation_receipt.json';ar=json.loads(ap.read_text())
    for name,h in ar['artifacts'].items():
        if sha(ann/name)!=h:raise ValueError('Changed Pfam annotation source')
    hits=defaultdict(list)
    with (ann/'raw_annotated_hits.tsv').open() as f:
        for x in csv.DictReader(f,delimiter='\t'):
            if x['sequence_id'] in identities:hits[identities[x['sequence_id']]].append(x)
    manifest_path=ROOT/'metadata/analysis_manifest.tsv';manifest={x['taxon_id']:x for x in read_table(manifest_path)}
    mappings=read_table(ROOT/'results/phylogeny/markers-full-v1/protein_mapping.tsv');proteins={x['taxon_id']:x['protein_id'] for x in mappings if x['marker']==marker}
    rows=[];partition={}
    for t,s in seqs.items():
        labels=Counter(x['pfam_name'] for x in hits[t]);group='BRF1_hit' if labels['BRF1'] else 'no_BRF1_hit';partition[t]=group
        rows.append({'taxon_id':t,'species_name':manifest[t]['species_name'],'protein_id':proteins[t],'protein_length':len(s),'sequence_sha256':hashlib.sha256(s.encode()).hexdigest(),'hit_group':group,'pfam_hits_json':json.dumps([{k:x[k] for k in ['hit_id','pfam_accession','pfam_name','pfam_type','alignment_start','alignment_end','hmm_coverage','domain_score']} for x in hits[t]])})
    pairs=[x for x in read_table(base/'pair_correspondence.tsv') if x['marker']==marker and x['genus_label']==genus];strata=defaultdict(list)
    for x in pairs:
        left=partition[x['taxon_a']];right=partition[x['taxon_b']];key='between_hit_groups' if left!=right else 'within_'+left;strata[key].append(x)
    summaries=[]
    for key,rr in sorted(strata.items()):
        out={'pair_stratum':key,'pairs':len(rr)}
        for method in ['global','local']:
            eligible=[x for x in rr if int(x[method+'_shared_codons'])>=100]
            out[method+'_pairs_at_least_100_codons']=len(eligible)
            out[method+'_median_aa_difference']=statistics.median(float(x[method+'_aa_difference_fraction']) for x in eligible) if eligible else ''
        summaries.append(out)
    if len(pairs)!=len(seqs)*(len(seqs)-1)//2:raise ValueError('Incomplete pair partition')
    a.output.mkdir(parents=True)
    for name,data in [('protein_domain_hits.tsv',rows),('domain_hit_pair_strata.tsv',summaries)]:
        with (a.output/name).open('w') as f:
            w=csv.DictWriter(f,list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    result={'status':'complete_aspergillus_marker_domain_hit_review','marker':marker,'genus_label':genus,'taxa':len(rows),'hit_groups':dict(Counter(partition.values())),'pair_strata':summaries,'realignment_receipt_sha256':sha(base/'receipt.json'),'domain_annotation_receipt_sha256':sha(ap),'manifest_sha256':sha(manifest_path),'script_sha256':sha(Path(__file__)),'interpretation':'Source Pfam hits distinguish BRF1-hit and no-BRF1-hit subsets; these are observed annotation patterns, not validated functions or resolved orthogroups. Pair divergence is stratified after seeing domain heterogeneity, so this is exploratory case diagnosis, not a confirmatory test. Mixed homologous protein types/hidden paralogy is a competing explanation requiring gene-family phylogeny and reconciliation. No domain gain/loss or structural acceleration is inferred; original BUSCO selections remain unchanged.','artifacts':{p.name:sha(p) for p in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
