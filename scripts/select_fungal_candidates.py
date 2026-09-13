#!/usr/bin/env python3
"""Select 500 provisional genomes for QC using explicit taxonomic coverage rules."""
import csv,json
from collections import Counter,defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader((ROOT/'metadata/species_assembly_candidates.tsv').open(),delimiter='\t'))
rows=[r for r in rows if r['phylum'] and r['kingdom']=='Fungi']
config=json.loads((ROOT/'config/sampling.json').read_text())
level={'Complete Genome':4,'Chromosome':3,'Scaffold':2,'Contig':1}
def quality(r):return (r['refseq_category']=='reference genome',r['catalog_source']=='refseq_fungi',level.get(r['assembly_level'],0),r['seq_rel_date'],r['assembly_accession'])
chosen={};why={}
def add(r,reason):chosen[r['species_taxid']]=r;why[r['species_taxid']]=reason
# More closely sampled genera provide candidate comparisons; ecological states remain unassigned.
genera=config['focal_genera']
for genus in genera:
 pool=sorted([r for r in rows if r['genus']==genus],key=quality,reverse=True)
 for r in pool[:config['focal_per_genus']]:add(r,'Focal genus coverage; ecological contrasts require independent trait curation')
# Prefer unrepresented orders, then families, then genera within each phylum.
def diverse(pool,n,reason):
 for _ in range(n):
  available=[r for r in pool if r['species_taxid'] not in chosen]
  if not available:break
  counts={k:Counter(r[k] for r in chosen.values() if r[k]) for k in ('order','family','genus')}
  def score(r):return tuple(-counts[k].get(r[k],0) if r[k] else -100000 for k in ('order','family','genus'))+quality(r)
  add(max(available,key=score),reason)
phyla=sorted({r['phylum'] for r in rows})
for phylum in phyla:
 if phylum not in ('Ascomycota','Basidiomycota'):
  pool=[r for r in rows if r['phylum']==phylum]
  diverse(pool,min(config['minor_phylum_cap'],len(pool)),'Broad coverage of sparsely sampled phylum; quality exceptions subject to QC')
while len(chosen)<config['fungal_target']:
 counts=Counter(r['phylum'] for r in chosen.values())
 phylum=min(('Ascomycota','Basidiomycota'),key=lambda p:counts[p]/config['dikarya_weights'][p])
 before=len(chosen)
 diverse([r for r in rows if r['phylum']==phylum],1,'Broad order/family/genus coverage within Dikarya')
 if len(chosen)==before:raise RuntimeError('Insufficient eligible species')
assert len(chosen)==config['fungal_target']
fields=next(csv.reader((ROOT/'metadata/sampling_manifest.tsv').open(),delimiter='\t'))
manifest=[]
for sid,r in sorted(chosen.items(),key=lambda x:(x[1]['phylum'],x[1]['order'],x[1]['taxon_name'])):
 base=r['ftp_path'].replace('ftp://','https://').rstrip('/');stem=base.rsplit('/',1)[-1]
 m={k:'' for k in fields};m.update(taxon_id='F'+sid,species_name=r['taxon_name'],species_taxid=sid,study_role='ingroup',lineage=';'.join(r[k] for k in ('phylum','class','order','family','genus') if r[k]),assembly_accession=r['assembly_accession'],annotation_source_version=r['annotation_provider']+' | '+r['annotation_name']+' | '+r['annotation_date'],proteome_url=r['provisional_proteome_url'],cds_url=f'{base}/{stem}_cds_from_genomic.fna.gz',genome_url=f'{base}/{stem}_genomic.fna.gz',assembly_level=r['assembly_level'],contamination_status='not_assessed',ecology='unknown',structure_coverage='not_assessed',inclusion_reason=why[sid],status='provisional_selected_for_QC')
 manifest.append(m)
with (ROOT/'metadata/fungal_sampling_draft.tsv').open('w') as out:
 w=csv.DictWriter(out,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(manifest)
summary={'selected_fungal_candidates':len(manifest),'counts_by_phylum':dict(Counter(r['phylum'] for r in chosen.values())),'distinct_orders':len({r['order'] for r in chosen.values() if r['order']}),'distinct_families':len({r['family'] for r in chosen.values() if r['family']}),'status':'Draft for full-scale QC, not final sample; no phylogenetic branch-length optimization or verified ecological transitions yet'}
(ROOT/'metadata/fungal_sampling_draft_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
