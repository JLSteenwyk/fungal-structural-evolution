#!/usr/bin/env python3
"""Add six reviewed published classifications with explicit sample linkage."""
import csv
import json
from pathlib import Path
from import_miyauchi_species_ecology import ROOT, sha, read

ADDITIONS = {
    'Amanita rubescens': ('ectomycorrhizal', 'Amanita_single_origin'),
    'Gautieria morchelliformis': ('ectomycorrhizal', 'Phallomycetidae_transition_review_pending'),
    'Hysterangium stoloniferum': ('ectomycorrhizal', 'Phallomycetidae_transition_review_pending'),
    'Ramaria rubella': ('saprotrophic', 'Phallomycetidae_transition_review_pending'),
    'Thelephora ganbajun': ('ectomycorrhizal', 'Thelephora_transition_review_pending'),
    'Thelephora terrestris': ('ectomycorrhizal', 'Thelephora_transition_review_pending'),
}


def main():
    previous=ROOT/'metadata/species_ecology_evidence_expanded_20260922.tsv'
    matches=ROOT/'metadata/miyauchi_species_ecology_matches.tsv'
    identity=ROOT/'metadata/miyauchi_ecology_sample_identity.tsv'
    previous_receipt=ROOT/'metadata/ecology_evidence_expansion_20260922_receipt.json'
    match_receipt=ROOT/'metadata/miyauchi_species_ecology_import_receipt.json'
    identity_receipt=ROOT/'metadata/miyauchi_ecology_sample_identity_receipt.json'
    assert sha(previous)==json.loads(previous_receipt.read_text())['artifacts'][previous.name]
    assert sha(matches)==json.loads(match_receipt.read_text())['artifacts'][matches.name]
    assert sha(identity)==json.loads(identity_receipt.read_text())['output_sha256']
    rows=read(previous);fields=list(rows[0]);names={r['species_name'] for r in rows}
    source={r['species_name']:r for r in read(matches)}
    linked={r['species_name']:r for r in read(identity)}
    for name,(state,group) in ADDITIONS.items():
        assert name not in names
        s=source[name];link=linked[name]
        assert link['status']=='same_biosample_and_wgs_project'
        assert s['selected_assembly_accession']==link['selected_assembly_accession']
        assert s['source_ecology']=={'ectomycorrhizal':'Ectomycorrhizae','saprotrophic':'Saprotroph'}[state]
        note='Published classification linked by identical BioSample and WGS project to the selected assembly; assembly/annotation identity and experimental phenotype validation are not established.'
        if name=='Ramaria rubella':
            note+=' Species source overrides the genus-level ECM candidate for this coding; source also names R. acris. Taxonomic-concept uncertainty remains explicit and omission sensitivity is required.'
        row=dict.fromkeys(fields,'')
        row.update(taxon_id=s['taxon_id'],assembly_accession=s['selected_assembly_accession'],
                   trait=rows[0]['trait'],species_name=name,state=state,source_species_name=s['source_species_name'],
                   source_doi=s['source_doi'],source_url=s['source_url'],
                   source_locator=f"Supplementary Dataset 2, {s['source_sheet']}, Excel row {s['source_excel_row']}; Dataset 1 sample linkage in metadata/miyauchi_ecology_sample_identity.tsv",
                   evidence='Species-level classification explicitly reported in the supplementary table.',
                   provisional_transition_group=group,scope_note=note,evidence_level='published_species_classification',
                   genus_candidate_primary_lifestyle=s['genus_candidate_primary_lifestyle'],
                   conflicts_with_genus_ECM_candidate=s['genus_ECM_disagreement'],
                   selected_isolate_experimentally_verified='False',
                   confirmatory_test_status='pending_taxonomy_phylogeny_and_replicated_transition_review')
        rows.append(row)
    path=ROOT/'metadata/species_ecology_evidence_sample_linked_20260922.tsv'
    with path.open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)
    result=dict(status='complete_sample_linked_species_classification_extension',statements=len(rows),
                additions=len(ADDITIONS),inputs={str(p.relative_to(ROOT)):sha(p) for p in
                (previous,matches,identity,previous_receipt,match_receipt,identity_receipt)},
                script_sha256=sha(Path(__file__)),table_sha256=sha(path),
                scope='Source-specific species classifications; source conflict retained, no experimental phenotype or independent-transition claim. Original 26 rows preserved.')
    (ROOT/'metadata/ecology_sample_linked_curation_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
