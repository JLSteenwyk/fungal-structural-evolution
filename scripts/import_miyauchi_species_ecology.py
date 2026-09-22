#!/usr/bin/env python3
"""Import exact-name species classifications, retaining conflicts and identity limits."""
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
import openpyxl

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    paths = {k:ROOT/v for k,v in {
        'downloads':'metadata/miyauchi_ecology_source_downloads.json',
        'manifest':'metadata/analysis_manifest.tsv',
        'genus':'metadata/ecology_candidates.tsv',
        'curated':'metadata/species_ecology_evidence_expanded_20260922.tsv',
    }.items()}
    downloads = json.loads(paths['downloads'].read_text())
    for record in downloads:
        assert sha(ROOT/record['path']) == record['sha256']
    record = next(r for r in downloads if 'MOESM6' in r['path'])
    source = ROOT/record['path']
    workbook = openpyxl.load_workbook(source, read_only=True, data_only=True)
    sheet = workbook['S_Table 2']
    values = list(sheet.values)
    assert values[4][1]=='Species' and values[4][7]=='Ecology'
    source_rows = [(i,r) for i,r in enumerate(values[5:],6) if r[1] is not None]
    assert len(source_rows)==135
    taxa = defaultdict(list)
    for row in read(paths['manifest']):
        if row['study_role']=='ingroup':
            taxa[row['species_name']].append(row)
    genus = {r['taxon_id']:r for r in read(paths['genus'])}
    curated = {r['taxon_id']:r for r in read(paths['curated'])}
    counts = Counter(str(r[1]).strip() for _,r in source_rows)
    output = []
    dispositions = []
    for number, row in source_rows:
        raw_name = str(row[1])
        name = raw_name.strip()
        matches = taxa.get(name, [])
        status = 'exact_unique_species_name' if len(matches)==1 and counts[name]==1 else (
            'unmatched_species_name' if not matches else 'ambiguous_species_name')
        dispositions.append(dict(source_excel_row=number, source_species_name=raw_name,
                                 source_ecology=row[7], match_status=status))
        if status!='exact_unique_species_name':
            continue
        taxon = matches[0]
        existing = curated.get(taxon['taxon_id'], {})
        candidate = genus[taxon['taxon_id']]
        # Only flag the explicit ECM/non-ECM classification discrepancy here.
        # Other category strings remain available for review without an ontology guess.
        ecm_disagreement = candidate['candidate_primary_lifestyle']=='ectomycorrhizal' and row[7] in ('Saprotroph','Wood decayer')
        output.append(dict(taxon_id=taxon['taxon_id'],species_name=name,
                           selected_assembly_accession=taxon['assembly_accession'],
                           source_species_name=raw_name,source_ecology=row[7],
                           source_jgi_id=row[0],source_strain_description=row[8],
                           source_sheet=sheet.title,source_excel_row=number,
                           source_doi='https://doi.org/10.1038/s41467-020-18795-w',
                           source_url=record['url'],source_sha256=record['sha256'],
                           previous_curated_state=existing.get('state',''),
                           genus_candidate_primary_lifestyle=candidate['candidate_primary_lifestyle'],
                           genus_ECM_disagreement=ecm_disagreement,
                           taxon_identity_flags=candidate['taxon_identity_flags'],
                           selected_isolate_equivalence='not_verified',
                           ecological_effect_eligibility='pending_source_conflict_identity_and_phylogenetic_review'))
    artifacts = {}
    for name, rows in [('miyauchi_species_ecology_matches.tsv', output),
                       ('miyauchi_species_ecology_source_disposition.tsv', dispositions)]:
        path = ROOT/'metadata'/name
        with path.open('w') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
            writer.writeheader();writer.writerows(rows)
        artifacts[name]=sha(path)
    result=dict(status='complete_exact_species_name_source_import_pending_review',
                source_records=len(source_rows),matched_records=len(output),
                previously_curated=sum(bool(r['previous_curated_state']) for r in output),
                newly_linked_species=sum(not r['previous_curated_state'] for r in output),
                genus_ECM_disagreements=sum(r['genus_ECM_disagreement'] for r in output),
                source_state_counts=dict(Counter(r['source_ecology'] for r in output)),
                inputs={k:sha(v) for k,v in paths.items()},script_sha256=sha(Path(__file__)),
                openpyxl_version=openpyxl.__version__,artifacts=artifacts,
                scope='Published species-name classifications with original labels and source rows. Whitespace trimming only; no synonym, strain or assembly equivalence inferred. Not automatically merged into curated states or ecological tests; conflicts and noncommensurate classifications require review.')
    (ROOT/'metadata/miyauchi_species_ecology_import_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
