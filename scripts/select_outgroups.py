#!/usr/bin/env python3
"""Build explicit 25-taxon outgroup QC draft, retaining source distinctions."""
import csv,json,hashlib
from pathlib import Path
from catalog_inventory import parse_catalog
ROOT=Path(__file__).resolve().parents[1]
config=json.loads((ROOT/'config/outgroups.json').read_text())
rows=[];receipts=json.loads((ROOT/'metadata/source_receipts.json').read_text())
for source in ('refseq','genbank'):
 for group in ('protozoa','invertebrate'):
  path=ROOT/f'data/raw/{source}_{group}_assembly_summary.txt';content=path.read_bytes()
  receipts[f'{source}_{group}']={'url':f'https://ftp.ncbi.nlm.nih.gov/genomes/{source}/{group}/assembly_summary.txt','path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(content).hexdigest(),'bytes':len(content)}
  for r in parse_catalog(content.decode()):r['catalog_source']=source;rows.append(r)
fields=next(csv.reader((ROOT/'metadata/sampling_manifest.tsv').open(),delimiter='\t'))
result=[]
def blank(name,lineage,reason):
 m={k:'' for k in fields};m.update(species_name=name,study_role='outgroup',lineage=lineage,contamination_status='not_assessed',ecology='unknown',structure_coverage='not_assessed',inclusion_reason=reason,status='provisional_selected_for_QC');return m
for name,lineage,reason in config['ncbi']:
 eligible=[r for r in rows if (r['organism_name']==name or r['organism_name'].startswith(name+' ')) and r['annotation_name'] not in ('','na') and r['version_status']=='latest']
 if not eligible:raise ValueError('No annotated genome: '+name)
 def score(r):return (r['refseq_category']=='reference genome',r['catalog_source']=='refseq',r['seq_rel_date'],r['assembly_accession'])
 r=max(eligible,key=score);base=r['ftp_path'].replace('ftp://','https://').rstrip('/');stem=base.rsplit('/',1)[-1]
 m=blank(name,lineage,reason);m.update(taxon_id='O'+r['species_taxid'],species_taxid=r['species_taxid'],assembly_accession=r['assembly_accession'],annotation_source_version=r['annotation_name']+' | '+r['annotation_date'],proteome_url=f'{base}/{stem}_protein.faa.gz',genome_url=f'{base}/{stem}_genomic.fna.gz',cds_url=f'{base}/{stem}_cds_from_genomic.fna.gz',assembly_level=r['assembly_level']);result.append(m)
external=json.loads((ROOT/'metadata/external_genome_receipts.json').read_text())
for aid,name,lineage in config['external']:
 files=[r for r in external if r['article_id']==aid]
 pep=next(r for r in files if 'proteins' in r)
 genome=next(r for r in files if 'gDNA' in r['name'] or '.genome.' in r['name'])
 species_ids={r['species_taxid'] for r in rows if r['organism_name']==name}
 m=blank(name,lineage,'Genome-derived unicellular holozoan annotation; publisher genome/proteome bundle')
 m.update(taxon_id='OFS'+str(aid),species_taxid=next(iter(species_ids)) if len(species_ids)==1 else '',assembly_accession='figshare:'+str(aid)+'.v'+str(pep['article_version']),annotation_source_version='Published Figshare bundle; see external_genome_receipts.json',proteome_url=pep['url'],genome_url=genome['url'])
 cds=next((r for r in files if '.cds.' in r['name']),None)
 if cds:m['cds_url']=cds['url']
 result.append(m)
assert len(result)==25 and len({r['species_name'] for r in result})==25
for filename,subset in [('outgroup_sampling_draft.tsv',result),('outgroup_ncbi_download_manifest.tsv',result[:20])]:
 with (ROOT/'metadata'/filename).open('w') as out:
  w=csv.DictWriter(out,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(subset)
(ROOT/'metadata/source_receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')
print('25 distinct outgroup taxa: 20 NCBI assemblies plus 5 verified publisher genome bundles')
