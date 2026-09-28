"""Add assembly-linked BioSample ecology while preserving all45 existing rows."""
import csv
import json
from pathlib import Path
from bs4 import BeautifulSoup
from ancestral_chain_attempt import sha,write_json


def read(path):
    with open(path) as h:return list(csv.DictReader(h,delimiter='\t'))


def main():
    config_path=Path('config/aphelid_ecology_evidence_20260928.json')
    config=json.loads(config_path.read_text())
    for path,h in config['pins'].items():assert sha(path)==h
    for path,source in config['sources'].items():assert sha(path)==source['sha256']
    original=read('metadata/species_ecology_evidence_reviewed_20260927.tsv');assert len(original)==45
    manifest={r['taxon_id']:r for r in read('metadata/analysis_manifest.tsv')}
    fields=list(original[0]);added=[];links=[]
    for spec in config['records']:
        selected=manifest[spec['taxon_id']]
        assert selected['study_role']=='ingroup'
        assert selected['species_name']==spec['species_name'] and selected['assembly_accession']==spec['assembly_accession']
        report={line[2:].split(':',1)[0].strip():line.split(':',1)[1].strip() for line in Path(spec['paths']['assembly.txt']).read_text().splitlines() if line.startswith('# ') and ':' in line}
        assert report['GenBank assembly accession']==spec['assembly_accession']
        assert report['BioSample']==spec['biosample'] and report['BioProject']==spec['bioproject']
        assert report['Infraspecific name']=='strain='+spec['strain']
        assert report['Organism name'].startswith(spec['species_name']+' (')
        sample=BeautifulSoup(Path(spec['paths']['biosample.html']).read_text(),'html.parser').get_text(' ',strip=True)
        project=BeautifulSoup(Path(spec['paths']['bioproject.html']).read_text(),'html.parser').get_text(' ',strip=True)
        for text in [sample,project]:
            assert spec['species_name']+' strain '+spec['strain'] in text
            assert spec['host'] in text and 'parasitoid' in text and spec['bioproject'] in text
        assert spec['biosample'] in sample and 'axenic cultures of '+spec['host'] in sample
        row=dict.fromkeys(fields,'');row.update(taxon_id=spec['taxon_id'],assembly_accession=spec['assembly_accession'],
            species_name=spec['species_name'],source_species_name=spec['species_name'],trait='trophic_role_as_reported_by_sample_submitter',
            state='algal_parasitoid',source_url=config['sources'][spec['paths']['biosample.html']]['url'],
            source_locator='BioSample Description and Attributes; assembly report BioSample, BioProject and strain fields',
            evidence='Submitter reports this strain parasitizing Scenedesmus dimorphus and cultivation with axenic host cultures.',
            provisional_transition_group='unassigned_'+spec['taxon_id'],
            scope_note='Selected assembly linked to named BioSample and strain. Submitter annotation; no independent experimental replication or inferred transition.',
            evidence_level='assembly_linked_submitter_sample_annotation',selected_isolate_experimentally_verified='False',
            confirmatory_test_status='pending_phylogeny_and_replicated_transition_review')
        assert spec['taxon_id'] not in {r['taxon_id'] for r in original+added};added.append(row)
        links.append({k:spec[k] for k in ['taxon_id','species_name','assembly_accession','strain','biosample','bioproject','host']})
    output=Path('metadata/species_ecology_evidence_aphelids_20260928.tsv')
    with output.open('x') as h:
        w=csv.DictWriter(h,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(original+added)
    sample_path=Path('metadata/aphelid_ecology_sample_links_20260928.tsv')
    with sample_path.open('x') as h:
        w=csv.DictWriter(h,list(links[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(links)
    reread=read(output);assert reread[:45]==original and reread[45:]==added and len(reread)==47
    write_json(Path('metadata/aphelid_ecology_evidence_receipt_20260928.json'),dict(status='complete_two_assembly_linked_aphelid_ecology_additions',preserved_rows=45,added_rows=2,total_rows=47,config_sha256=sha(config_path),script_sha256=sha(__file__),artifacts={str(p):sha(p) for p in [output,sample_path]},scope=config['scope']))
    print('Preserved45 species rows; added2 assembly-linked submitter ecological statements')

if __name__=='__main__':main()
