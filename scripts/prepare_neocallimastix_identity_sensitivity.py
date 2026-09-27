#!/usr/bin/env python3
"""Prepare all single-representative sensitivities for a literature-supported complex."""
import csv,json,xml.etree.ElementTree as ET
from datetime import datetime,timezone
from pathlib import Path
from Bio import SeqIO
from screen_duplication_domain_alignment_coverage import sha

GROUP={'F1754190':('GCA_002104975.1','G1'),'F2767002':('GCA_016946835.1','sp3'),'F3108450':('GCA_050613775.1','G3')}
OUTPUT=Path('results/phylogeny/neocallimastix-identity-inputs-20260927-v1')
DATA=Path('data/taxonomy/neocallimastix-review-20260927-v1')


def rows(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    mp=Path('metadata/analysis_manifest.tsv');catalog_path=Path('metadata/species_assembly_candidates.tsv')
    manifest=rows(mp);taxa={r['taxon_id']:r for r in manifest};catalog={r['assembly_accession']:r for r in rows(catalog_path)}
    assert len(taxa)==len(manifest)==526
    sources={str(mp):sha(mp),str(catalog_path):sha(catalog_path)};articles={}
    for name,pmc,doi in [('article.xml','PMC12341934','10.1093/g3journal/jkaf137'),('phylogeny_article.xml','PMC8399178','10.3390/microorganisms9081655')]:
        path=DATA/name;root=ET.parse(path).getroot();ids={x.get('pub-id-type'):x.text for x in root.findall('.//article-id')};assert ids['pmcid']==pmc and ids['doi']==doi
        articles[name]=' '.join(root.itertext());sources[str(path)]=sha(path)
    assert 'PRJNA1052201' in articles['article.xml'] and 'JBODTK000000000' in articles['article.xml'] and 'Neocon1' in articles['article.xml']
    assert 'conspecific' in articles['phylogeny_article.xml'] and all(word in articles['phylogeny_article.xml'] for word in ['californiae','lanati','cameroonii'])
    evidence=[]
    for tid,(accession,strain) in GROUP.items():
        row=taxa[tid];c=catalog[accession];assert row['assembly_accession']==accession and c['infraspecific_name']=='strain='+strain and row['study_role']=='ingroup'
        evidence.append(dict(taxon_id=tid,manifest_name=row['species_name'],assembly_accession=accession,strain=strain,bioproject=c['bioproject'],biosample=c['biosample'],interpretation='Literature-supported potential species-complex membership; not a formal synonymy decision or verified species collapse.'))
    assert catalog[GROUP['F3108450'][0]]['bioproject']=='PRJNA1052201' and catalog[GROUP['F3108450'][0]]['wgs_master'].split('.')[0]=='JBODTK000000000'
    OUTPUT.mkdir(parents=True,exist_ok=False);members=[];matrices=[];outgroups={t for t,r in taxa.items() if r['study_role']=='outgroup'};assert len(outgroups)==25
    for representative in sorted(GROUP):
        excluded=set(GROUP)-{representative};keep=set(taxa)-excluded
        assert len(keep)==524 and outgroups<=keep
        for tid in sorted(taxa):members.append(dict(policy='retain_'+representative,representative=representative,taxon_id=tid,retained=int(tid in keep),species_complex_member=int(tid in GROUP),study_role=taxa[tid]['study_role']))
    for alignment,folder in [('profile','profile-matrix-50-v1'),('mafft','mafft-matrix-50-v1')]:
        root=Path('results/phylogeny')/folder;rp=root/'receipt.json';r=json.loads(rp.read_text());path=root/'matrix.faa'
        assert sha(path)==r['artifacts']['matrix.faa'];sources[str(rp)]=sha(rp);sources[str(path)]=sha(path)
        records=list(SeqIO.parse(path,'fasta'));original={r.id:str(r.seq) for r in records};assert len(original)==len(records)==526 and set(original)==set(taxa)
        for representative in sorted(GROUP):
            excluded=set(GROUP)-{representative};selected=[r for r in records if r.id not in excluded];p=OUTPUT/(alignment+'-retain_'+representative+'.faa');SeqIO.write(selected,p,'fasta')
            observed={r.id:str(r.seq) for r in SeqIO.parse(p,'fasta')};assert observed=={k:v for k,v in original.items() if k not in excluded} and outgroups<=set(observed)
            matrices.append(dict(alignment=alignment,representative=representative,taxa=524,fungal_entries=499,outgroups=25,columns=len(selected[0].seq),path=p.name,sha256=sha(p)))
    p=OUTPUT/'membership.tsv'
    with p.open('w') as f:w=csv.DictWriter(f,list(members[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(members)
    assert rows(p)==[{k:str(v) for k,v in r.items()} for r in members]
    (OUTPUT/'assembly_evidence.json').write_text(json.dumps(evidence,indent=2)+'\n')
    for p,digest in sources.items():assert sha(p)==digest
    result=dict(status='complete_neocallimastix_identity_sensitivity_inputs',source_hashes=sources,script_sha256=sha(__file__),created_utc=datetime.now(timezone.utc).isoformat(),group_taxa=sorted(GROUP),matrices=matrices,articles=[dict(pmc='PMC12341934',doi='10.1093/g3journal/jkaf137',url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12341934/fullTextXML',path=str(DATA/'article.xml')),dict(pmc='PMC8399178',doi='10.3390/microorganisms9081655',url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8399178/fullTextXML',path=str(DATA/'phylogeny_article.xml'))],scope='Six exact aligned-sequence subsets retain each of three potential conspecific representatives in turn. All 25 outgroups and all other taxa retained; original data and active analyses unchanged. Names and taxids do not prove unique species. Literature motivates sensitivity, not automatic synonymy; genome-wide delimitation and phylogenetic comparisons remain pending. Other uncertain labels and hybrids remain unresolved.',artifacts={p.name:sha(p) for p in OUTPUT.iterdir()})
    (OUTPUT/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']},indent=2))


if __name__=='__main__':main()
