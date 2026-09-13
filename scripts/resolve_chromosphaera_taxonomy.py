#!/usr/bin/env python3
"""Verify Chromosphaera's deposited taxonomy identity using ENA and NCBI records."""
import csv
import hashlib
import io
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import xml.etree.ElementTree as ET
ROOT = Path(__file__).resolve().parents[1]

def main():
    output = ROOT / 'data/taxonomy/chromosphaera-resolution-v1'
    if output.exists():
        raise FileExistsError('Use existing pinned evidence or a new version')
    urls = {
        'study_runs.tsv': 'https://www.ebi.ac.uk/ena/portal/api/search?' + urlencode({'result':'read_run', 'query':'study_accession="PRJNA360047"', 'fields':'run_accession,sample_accession,scientific_name,tax_id,sample_alias,sample_title', 'format':'tsv'}),
        'taxonomy.xml': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=1932427&retmode=xml',
        'article.xml': 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5560861/fullTextXML'}
    bodies = {}
    for name, url in urls.items():
        with urlopen(url, timeout=60) as response:
            bodies[name] = response.read()
    taxonomy = ET.fromstring(bodies['taxonomy.xml']).find('Taxon')
    if taxonomy.findtext('TaxId') != '1932427' or taxonomy.findtext('ScientificName') != 'Ichthyosporea sp. XGB-2017a':
        raise ValueError('Unexpected taxonomy identity')
    aliases = {n.findtext('DispName'): n.findtext('ClassCDE') for n in taxonomy.findall('./OtherNames/Name')}
    if aliases.get('Chromosphaera perkinsii') != 'unpublished name':
        raise ValueError('Expected explicit NCBI name linkage missing')
    rows = list(csv.DictReader(io.StringIO(bodies['study_runs.tsv'].decode()), delimiter='\t'))
    selected = [r for r in rows if r['sample_accession'] == 'SAMN06200032']
    if len(selected) != 2 or any(r['tax_id'] != '1932427' or r['sample_alias'] != 'Chromosphaera perkinsii DNA-seq' for r in selected):
        raise ValueError('Genome sample linkage differs')
    article = ET.fromstring(bodies['article.xml'])
    if 'PRJNA360047' not in ''.join(article.itertext()) or 'Chromosphaera' not in ''.join(article.itertext()):
        raise ValueError('Genome publication linkage missing')
    output.mkdir(parents=True)
    for name, body in bodies.items():
        (output / name).write_bytes(body)
    receipt = {'status':'verified_deposition_taxonomy_link', 'taxon_id':'OFS5426494',
        'project_species_name':'Chromosphaera perkinsii', 'species_taxid':'1932427',
        'ncbi_scientific_name':'Ichthyosporea sp. XGB-2017a', 'ncbi_alias_class':'unpublished name',
        'biosample':'SAMN06200032', 'bioproject':'PRJNA360047', 'run_accessions':sorted(r['run_accession'] for r in selected),
        'paper_doi':'10.7554/eLife.26036', 'retrieved_utc':datetime.now(timezone.utc).isoformat(),
        'sources': {name:{'url':urls[name], 'path':str((output/name).relative_to(ROOT)), 'sha256':hashlib.sha256(body).hexdigest()} for name,body in bodies.items()},
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'interpretation':'NCBI explicitly links the project name to this species-rank record, and the original genome BioSample independently links the same taxid and name. Keep the publication name and database name distinct. NCBI alias-class wording is recorded, not interpreted as a ruling on nomenclatural validity. Historical sampling manifests remain unchanged; versioned query override supplies the verified identifier.'}
    path = ROOT/'metadata/chromosphaera_taxonomy_resolution.json'
    path.write_text(json.dumps(receipt,indent=2)+'\n')
    config = ROOT/'config/taxonomy_overrides.json'
    if config.exists():
        raise FileExistsError('Review existing overrides instead of overwriting')
    config.write_text(json.dumps({'OFS5426494':{'species_name':'Chromosphaera perkinsii','species_taxid':'1932427',
        'evidence_path':str(path.relative_to(ROOT)), 'evidence_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}},indent=2)+'\n')
    print(json.dumps(receipt,indent=2))

if __name__ == '__main__':
    main()
