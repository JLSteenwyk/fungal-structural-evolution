#!/usr/bin/env python3
"""Verify cached primary identity sources and trace all sampled amphioxus markers.

Source URLs and retrieval hashes reside beside the cached XML/assembly reports.
Writes exclusive new evidence and a display/provenance overlay; no frozen
manifests, query identifiers, native inputs or inference outputs are edited.
"""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
from urllib.parse import unquote
import xml.etree.ElementTree as ET

from audit_selected_taxon_identity_snapshot_v2 import rows, sha, table


def attributes(text):
    parts = [s.split('=', 1) for s in text.split(';') if s]
    assert all(len(p) == 2 for p in parts)
    assert len({p[0] for p in parts}) == len(parts)
    return {key: unquote(value) for key,value in parts}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sources', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    p.add_argument('--overlay', type=Path, required=True)
    a = p.parse_args()
    assert not any(path.exists() for path in [a.output,a.receipt,a.overlay])
    manifest = Path('metadata/analysis_manifest.tsv')
    inventory = Path('results/taxonomy/selected-identity-snapshot-20261003-v2/taxon_identity.tsv')
    inventory_receipt = Path('metadata/selected_taxon_identity_snapshot_completed_20261003_v2.json')
    mapping = Path('results/phylogeny/markers-full-v1/protein_mapping.tsv')
    gff = Path('data/annotations/GCF_019207075.1_genomic.gff.gz')
    proteins = Path('data/qc_proteomes/O2700040.faa')
    stats_receipt = Path('data/assembly_statistics/GCF_019207075.1.receipt.json')
    taxon = 'O2700040'
    selected = [r for r in rows(manifest) if r['taxon_id'] == taxon]
    assert len(selected) == 1
    selected = selected[0]
    assert selected['assembly_accession'] == 'GCF_019207075.1'
    assert selected['species_name'] == 'Branchiostoma floridae'
    assert selected['species_taxid'] == '2700040' and selected['study_role'] == 'outgroup'
    audit = json.loads(inventory_receipt.read_text())
    assert audit['artifacts'][str(inventory)] == sha(inventory)
    identity = next(r for r in rows(inventory) if r['taxon_id'] == taxon)
    assert identity['snapshot_scientific_name'] == 'Branchiostoma floridae x Branchiostoma belcheri'
    assert identity['catalogue_paired_accession'] == 'GCA_019207075.1'
    assert identity['catalogue_paired_comparison'] == 'identical'
    assert identity['catalogue_biosample'] == 'SAMN13907882'
    assert identity['catalogue_pubmed_ids'] == '36867684'
    bindings = {str(path):sha(path) for path in [Path(__file__),
        Path('scripts/audit_selected_taxon_identity_snapshot_v2.py'),manifest,inventory,
        inventory_receipt,mapping,gff,proteins,stats_receipt]}
    retrievals = []
    for name in ['retrieval_attempts.json','assembly_report_retrieval_attempts.json']:
        path = a.sources/name
        bindings[str(path)] = sha(path)
        for record in json.loads(path.read_text()):
            assert record['status'] == 'downloaded_not_yet_content_validated'
            assert record['http_status'] == 200 and sha(record['path']) == record['sha256']
            bindings[record['path']] = record['sha256']
            retrievals.append(record)
    project = ET.parse(a.sources/'bioproject.xml').getroot()
    assert project.find('.//ArchiveID').attrib['accession'] == 'PRJNA706157'
    description = project.findtext('.//ProjectDescr/Description')
    assert description == 'This is the haploid genome assembly of Branchiostoma belcheri from the Bb/Bf hybrid'
    assert project.find('.//Organism').attrib['taxID'] == '2700040'
    # RefSeq catalogue has its own annotation project; verify explicit peer link.
    assert any(n.attrib.get('accession') == identity['catalogue_bioproject']
               for n in project.findall('.//ProjectLinks//ProjectIDRef'))
    sample = ET.parse(a.sources/'biosample.xml').getroot().find('BioSample')
    assert sample.attrib['accession'] == 'SAMN13907882'
    assert sample.find('./Description/Organism').attrib['taxonomy_id'] == '2700040'
    assert sample.findtext('./Description/Organism/OrganismName') == identity['snapshot_scientific_name']
    sample_name = next(n.text for n in sample.findall('./Ids/Id') if n.attrib.get('db_label')=='Sample name')
    assert sample_name == 'bbbf_hybrid_bb'
    sample_project = next(n.attrib['label'] for n in sample.findall('./Links/Link') if n.attrib.get('target')=='bioproject')
    article = ET.parse(a.sources/'article.xml').getroot()
    identifiers = {n.attrib['pub-id-type']:n.text for n in article.findall('.//article-id')}
    assert identifiers['pmid'] == '36867684' and identifiers['doi'] == '10.1073/pnas.2201504120'
    paragraphs = [''.join(n.itertext()) for n in article.findall('.//p')]
    method = next(s for s in paragraphs if 'We assigned a contig to either parental species' in s)
    assert 'partitioned the species-specific haploid reads' in method
    assert any(sample_project in s for s in paragraphs)
    report_paths = [a.sources/(acc+'_Bb-hap_assembly_report.txt') for acc in ['GCA_019207075.1','GCF_019207075.1']]
    assert report_paths[0].read_bytes() == report_paths[1].read_bytes()
    headers, sequences = {}, {}
    for line in report_paths[1].read_text().splitlines():
        if line.startswith('# '):
            if ': ' in line:
                key,value = line[2:].split(':',1)
                assert key not in headers
                headers[key] = value.strip()
            continue
        if not line or line.startswith('#'): continue
        fields = line.split('\t')
        assert len(fields) == 10
        name,role,molecule,location,genbank,relationship,refseq,unit,length,ucsc = fields
        assert refseq not in sequences and relationship == '='
        sequences[refseq] = dict(sequence_id=refseq,genbank_sequence_id=genbank,
            sequence_name=name,sequence_role=role,length=int(length))
    for key,value in [('GenBank assembly accession','GCA_019207075.1'),
                      ('RefSeq assembly accession','GCF_019207075.1'),
                      ('RefSeq assembly and GenBank assemblies identical','yes'),
                      ('BioSample','SAMN13907882'),('BioProject','PRJNA706157'),('Taxid','2700040')]:
        assert headers[key] == value
    old_stats = json.loads(stats_receipt.read_text())
    assert sha(old_stats['path']) == old_stats['sha256']
    bindings[old_stats['path']] = old_stats['sha256']
    for key,value in old_stats['headers'].items(): assert headers[key] == value, key
    assert len(sequences) == old_stats['all_assembly_metrics']['scaffold-count'] == 79
    assert sum(r['length'] for r in sequences.values()) == old_stats['all_assembly_metrics']['total-length']
    genome_sequences, segments, feature_counts = set(), defaultdict(list), Counter()
    metadata = []
    with gzip.open(gff,'rt') as stream:
        for line in stream:
            if line.startswith('#'):
                metadata.append(line.rstrip('\n'))
                continue
            fields = line.rstrip('\n').split('\t')
            assert len(fields) == 9
            sequence_id,source,feature,start,end,score,strand,phase,attrs = fields
            assert sequence_id in sequences
            assert 1 <= int(start) <= int(end) <= sequences[sequence_id]['length']
            genome_sequences.add(sequence_id); feature_counts[feature] += 1
            if feature != 'CDS': continue
            attrib = attributes(attrs)
            protein_id = attrib.get('protein_id')
            if protein_id:
                segments[protein_id].append(dict(sequence_id=sequence_id,start=int(start),
                    end=int(end),strand=strand,phase=phase,parent=attrib.get('Parent',''),
                    gene=attrib.get('gene','')))
    assert genome_sequences == set(sequences)
    assert '#!genome-build-accession NCBI_Assembly:GCF_019207075.1' in metadata
    assert '#!annotation-source NCBI RefSeq '+selected['annotation_source_version'].split(' | ')[0] in metadata
    from hashlib import sha256
    fasta, label, sequence = {}, None, []
    with proteins.open() as stream:
        for line in stream:
            if line.startswith('>'):
                if label is not None:
                    assert label not in fasta
                    fasta[label] = sha256(''.join(sequence).encode()).hexdigest()
                label, sequence = line[1:].split()[0], []
            else: sequence.append(line.strip())
    assert label not in fasta
    fasta[label] = sha256(''.join(sequence).encode()).hexdigest()
    all_mapping = rows(mapping)
    markers = sorted({r['marker'] for r in all_mapping})
    assert len(markers) == 125
    selected_markers = {r['marker']:r for r in all_mapping if r['taxon_id']==taxon}
    assert len(selected_markers) == 100
    output_rows = []
    for marker in markers:
        row = selected_markers.get(marker)
        coords = segments.get(row['protein_id'],[]) if row else []
        if row:
            assert fasta[row['protein_id']] == row['sequence_sha256']
            assert coords and len({c['sequence_id'] for c in coords}) == 1
            assert len({c['strand'] for c in coords}) == 1
        output_rows.append(dict(marker=marker,taxon_id=taxon,
            selection_status='selected' if row else 'not_selected',protein_id=row['protein_id'] if row else '',
            sequence_sha256=row['sequence_sha256'] if row else '',
            annotation_sequence_id=coords[0]['sequence_id'] if coords else '',
            annotation_genbank_sequence_id=sequences[coords[0]['sequence_id']]['genbank_sequence_id'] if coords else '',
            annotation_gene=';'.join(sorted({c['gene'] for c in coords})),
            cds_segment_count=len(coords),cds_segments_json=json.dumps(coords,separators=(',',':')),
            deposition_parental_assignment='Branchiostoma belcheri principal haplotype' if row else '',
            independent_locus_parental_purity='unverified'))
    a.output.mkdir(parents=True,exist_ok=False)
    table(a.output/'all_marker_provenance.tsv',output_rows)
    table(a.output/'assembly_sequences.tsv',list(sequences.values()))
    for path,h in bindings.items(): assert sha(path) == h, ('source changed',path)
    evidence = dict(status='verified_amphioxus_deposition_identity_and_all_marker_provenance',
        checked_utc=datetime.now(timezone.utc).isoformat(),taxon_id=taxon,
        frozen_project_name=selected['species_name'],selected_assembly=selected['assembly_accession'],
        genbank_assembly='GCA_019207075.1',deposition_taxid='2700040',
        deposition_organism_name=identity['snapshot_scientific_name'],biosample='SAMN13907882',
        biosample_name=sample_name,biosample_related_project=sample_project,
        assembly_project='PRJNA706157',refseq_annotation_project=identity['catalogue_bioproject'],
        principal_haplotype_parental_assignment='Branchiostoma belcheri',
        selected_genome_is_proven_mixed_parental_assembly=False,
        independent_genome_wide_parental_purity_verified=False,
        paper_doi=identifiers['doi'],paper_pmid=identifiers['pmid'],
        assembly_report_matches_original_statistics=True,assembly_sequences=len(sequences),
        gff_feature_counts=dict(feature_counts),annotated_protein_identifiers=len(segments),
        normalized_proteome_identifiers=len(fasta),full_marker_panel=len(markers),
        selected_markers=100,not_selected_markers=25,
        all_selected_marker_sequences_and_coordinates_checked=True,
        original_manifest_unchanged=True,source_hashes=bindings,
        artifacts={str(path):sha(path) for path in a.output.iterdir()},
        retrievals=retrievals,
        resource_plan=dict(cpus=1,address_space_gib=4,cpu_seconds=180,wall_seconds=300,
            output_allowance_mib=16,planning_runtime_seconds=[1,60],gpu=False,charges=False),
        scope='Exact RefSeq/GenBank/BioSample/BioProject links and identical assembly reports verified, including original cached statistics and complete annotation sequence coordinates. Primary study explains parental read partitioning; submitter assigns this principal haplotype to B.belcheri. All125marker slots retained, with100selected protein hashes and complete CDS coordinate provenance. This fixes an evidence/display identity, not species delimitation or independent locus parental purity. No retranslation, relabeling of frozen inputs, taxid replacement, taxon pruning or native inference.')
    with a.receipt.open('x') as stream: stream.write(json.dumps(evidence,indent=2)+'\n')
    overlay = dict(scope='Versioned evidence/display review; adoption requires explicit consumer use. Frozen IDs, sequences and manifests retain their historical identities.',
        records={taxon:dict(assembly_accession=selected['assembly_accession'],
            frozen_species_name=selected['species_name'],
            display_name='Branchiostoma belcheri (hybrid-derived principal haplotype)',
            ncbi_deposition_name=identity['snapshot_scientific_name'],ncbi_deposition_taxid='2700040',
            biological_parent_assignment_source='submitter deposition; independently verified locus purity remains open',
            review_class='species_label_disagrees_with_deposited_principal_haplotype',
            evidence_path=str(a.receipt),evidence_sha256=sha(a.receipt),
            analysis_policy='Retain this outgroup and its stable identifier with explicit assembly provenance; evaluate exclusion/reference sensitivity where relevant. Do not treat the hybrid source specimen as proof that this separated haplotype mixes parental genomes.')})
    with a.overlay.open('x') as stream: stream.write(json.dumps(overlay,indent=2)+'\n')
    print(json.dumps({k:v for k,v in evidence.items() if k not in ['source_hashes','artifacts','retrievals']},indent=2))


if __name__ == '__main__':
    main()
