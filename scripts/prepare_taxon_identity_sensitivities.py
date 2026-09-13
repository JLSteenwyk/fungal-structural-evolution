#!/usr/bin/env python3
"""Prepare hybrid/uncertain-label exclusion matrices without redefining species identity."""
import argparse,csv,json
from collections import Counter
from pathlib import Path
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT,sha,read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable sensitivity directory')
    mp=ROOT/'metadata/analysis_manifest.tsv';manifest=read_table(mp);taxa={r['taxon_id']:r for r in manifest}
    if len(taxa)!=len(manifest):raise ValueError('Duplicate taxon identity')
    lp=ROOT/'metadata/taxon_label_review.tsv';lr=ROOT/'metadata/taxon_label_review_receipt.json';lrec=json.loads(lr.read_text())
    if lrec['manifest_sha256']!=sha(mp) or lrec['artifact_sha256']!=sha(lp):raise ValueError('Changed label screen')
    labels={r['taxon_id']:r for r in read_table(lp)}
    cp=ROOT/'config/curated_hybrid_evidence.json';curation=json.loads(cp.read_text());ap=ROOT/'metadata/species_assembly_candidates.tsv';assemblies={r['assembly_accession']:r for r in read_table(ap)}
    evidence=[];hybrids=set()
    for r in curation['records']:
        t=taxa[r['taxon_id']];catalog=assemblies[r['assembly_accession']]
        if (t['assembly_accession']!=r['assembly_accession'] or t['species_name']!=r['species_name'] or catalog['infraspecific_name']!='strain='+r['strain'] or r['pubmed_id'] not in catalog['pubmed_id'].split(';')):raise ValueError('Assembly/strain/publication linkage differs')
        if r['taxon_id'] in hybrids:raise ValueError('Duplicate hybrid curation')
        hybrids.add(r['taxon_id']);evidence.append(r|{'catalog_assembly_type':catalog['assembly_type'],'catalog_bioproject':catalog['bioproject'],'catalog_biosample':catalog['biosample'],'interpretation':'Assembly representation field is not a biological ploidy or nonhybrid certificate. Parental assignment of individual genes remains pending.'})
    policies={'exclude_curated_hybrids':hybrids,'exclude_hybrids_and_uncertain_labels':hybrids|set(labels)}
    a.output.mkdir(parents=True);write_table(a.output/'hybrid_evidence.tsv',evidence);membership=[];counts={}
    for policy,excluded in policies.items():
        keep=set(taxa)-excluded
        if {t for t in keep if taxa[t]['study_role']=='outgroup'}!={t for t in taxa if taxa[t]['study_role']=='outgroup'}:raise ValueError('Unexpected outgroup exclusion')
        counts[policy]={'retained':len(keep),'excluded':len(excluded),'roles':dict(Counter(taxa[t]['study_role'] for t in keep))}
        (a.output/(policy+'.tips.txt')).write_text('\n'.join(sorted(keep))+'\n')
        for t in sorted(taxa):membership.append({'policy':policy,'taxon_id':t,'species_name':taxa[t]['species_name'],'study_role':taxa[t]['study_role'],'retained':t in keep,'curated_hybrid':t in hybrids,'label_flags':labels.get(t,{}).get('flags',''),'interpretation':'Sensitivity membership only; retention does not establish unique accepted species or absence of unrecognized hybridization.'})
    write_table(a.output/'taxon_membership.tsv',membership);sources={};matrices=[]
    for name,folder in [('profile','profile-matrix-50-v1'),('mafft','mafft-matrix-50-v1')]:
        base=ROOT/'results/phylogeny'/folder;r=checked_receipt(base);sources[name]=sha(base/'receipt.json')
        records=list(SeqIO.parse(base/'matrix.faa','fasta'));ids=[r.id for r in records]
        if len(ids)!=len(set(ids)) or set(ids)!=set(taxa):raise ValueError('Source matrix taxon universe differs')
        for policy,excluded in policies.items():
            target=a.output/(name+'-'+policy+'.faa');selected=[r for r in records if r.id not in excluded];SeqIO.write(selected,target,'fasta')
            back={r.id:str(r.seq) for r in SeqIO.parse(target,'fasta')}
            if back!={r.id:str(r.seq) for r in selected}:raise ValueError('Sensitivity matrix readback differs')
            matrices.append({'alignment':name,'policy':policy,'taxa':len(selected),'columns':len(selected[0].seq),'path':target.name,'sha256':sha(target)})
    result={'status':'complete_taxon_identity_sensitivity_inputs','manifest_sha256':sha(mp),'label_receipt_sha256':sha(lr),'curation_sha256':sha(cp),'assembly_catalog_sha256':sha(ap),'source_matrix_receipt_sha256':sources,'script_sha256':sha(Path(__file__)),'curated_hybrids':len(hybrids),'policies':counts,'matrices':matrices,'interpretation':'Input preparation only. Sequences/columns remain exactly as in the full matrices; excluded taxa are retained in original data. No topology test, accepted species count, homeolog assignment or comprehensive hybrid screen is claimed.','artifacts':{f.name:sha(f) for f in a.output.iterdir() if f.is_file()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='artifacts'},indent=2))


if __name__=='__main__':main()
