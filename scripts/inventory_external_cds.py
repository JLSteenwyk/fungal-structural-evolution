#!/usr/bin/env python3
"""Record verified available CDS sources and extraction gaps for external taxa."""
import csv
import json
from pathlib import Path
from Bio import SeqIO
from audit_busco_gene_copies import ROOT, sha, read_table


def main():
    manifest_path=ROOT/'metadata/analysis_manifest.tsv'
    external_path=ROOT/'metadata/external_genome_receipts.json'
    derived_path=ROOT/'metadata/sanchytrid_coordinate_audit.json'
    manifest=read_table(manifest_path)
    external=json.loads(external_path.read_text())
    derived={r['taxon_id']:r for r in json.loads(derived_path.read_text())['taxa']}
    rows=[]
    for taxon in manifest:
        if taxon['cds_url'].startswith('https://ftp.ncbi.nlm.nih.gov/'):continue
        name=taxon['taxon_id'];path=None
        row={'taxon_id':name,'species_name':taxon['species_name']}
        if name in derived:
            source=derived[name];path=ROOT/source['cds_path'];expected=source['cds_sha256']
            row.update(status='available_previously_translation_verified_subset',source_url='',
                translation_scope='Exact translation for retained derived CDS; source audit exceptions remain excluded and explicit',
                source_exception_counts=source['status_counts'])
        elif taxon['cds_url']:
            matches=[r for r in external if r.get('url')==taxon['cds_url']]
            if len(matches)!=1:raise ValueError('Cannot resolve unique external CDS provenance')
            source=matches[0];path=ROOT/source['path'];expected=source['sha256']
            row.update(status='available_published_cds_translation_pending',source_url=taxon['cds_url'],translation_scope='Translation and protein identifier correspondence still require validation')
        else:
            row.update(status='annotation_guided_extraction_pending',source_url='',translation_scope='No ready CDS FASTA recorded; use deposited genome and exact annotation with strand/phase validation')
        if path:
            if sha(path)!=expected:raise ValueError('Changed external CDS source')
            identifiers=set();bases=0
            for record in SeqIO.parse(path,'fasta'):
                sequence=str(record.seq).upper()
                if record.id in identifiers or not sequence or not set(sequence)<=set('ACGTRYSWKMBDHVN'):raise ValueError('Invalid external CDS identity or alphabet')
                identifiers.add(record.id);bases+=len(sequence)
            if not identifiers:raise ValueError('Empty external CDS source')
            row.update(path=str(path.relative_to(ROOT)),sha256=expected,cds_records=len(identifiers),nucleotide_bases=bases)
        rows.append(row)
    result={'status':'complete_external_cds_source_inventory_not_complete_extraction','taxa':rows,
        'manifest_sha256':sha(manifest_path),'external_source_receipt_sha256':sha(external_path),
        'derived_cds_receipt_sha256':sha(derived_path),'script_sha256':sha(Path(__file__)),
        'interpretation':'Available source hashes, unique FASTA identifiers and DNA alphabet checked. Published CDS availability does not establish exact translation, orthology or codon-model suitability. Pending extraction remains explicit; no taxon removed.'}
    target=ROOT/'metadata/external_cds_source_inventory.json'
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
