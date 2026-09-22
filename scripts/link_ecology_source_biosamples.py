#!/usr/bin/env python3
"""Link published sample identifiers to exact selected assembly records."""
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
import openpyxl

ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle,delimiter='\t'))


def main():
    files={key:ROOT/path for key,path in {
        'source':'data/traits/miyauchi-2020/41467_2020_18795_MOESM5_ESM.xlsx',
        'downloads':'metadata/miyauchi_ecology_source_downloads.json',
        'matches':'metadata/miyauchi_species_ecology_matches.tsv',
        'match_receipt':'metadata/miyauchi_species_ecology_import_receipt.json',
        'assemblies':'metadata/assembly_candidates.tsv',
        'quality':'metadata/assembly_quality_metrics.tsv',
    }.items()}
    download=next(r for r in json.loads(files['downloads'].read_text()) if 'MOESM5' in r['path'])
    assert sha(files['source'])==download['sha256']
    imported=json.loads(files['match_receipt'].read_text())
    assert sha(files['matches'])==imported['artifacts'][files['matches'].name]
    matches=read(files['matches'])
    accessions={r['selected_assembly_accession'] for r in matches}
    assemblies=defaultdict(list)
    for row in read(files['assemblies']):
        if row['assembly_accession'] in accessions:
            assemblies[row['assembly_accession']].append(row)
    quality={r['assembly_accession']:r for r in read(files['quality'])}
    values=list(openpyxl.load_workbook(files['source'],read_only=True,data_only=True)['S_Table 1'].values)
    assert values[4][0]=='Species' and values[4][10].strip()=='BioSample'
    paper=defaultdict(list)
    for number,row in enumerate(values[5:34],6):
        assert row[0] and str(row[10]).startswith('SAMN')
        paper[str(row[0]).strip()].append((number,row))
    results=[]
    for match in matches:
        accession=match['selected_assembly_accession']
        candidates=assemblies[accession]
        assert len(candidates)==1,accession
        selected=candidates[0]
        q=quality[accession]
        assert q['statistics_status']=='verified'
        # The assembly catalogue uses "na" where the parsed report uses blank.
        catalog_sample = '' if selected['biosample']=='na' else selected['biosample']
        assert catalog_sample==q['biosample']
        entries=paper.get(match['species_name'],[])
        row=dict(taxon_id=match['taxon_id'],species_name=match['species_name'],
                 selected_assembly_accession=accession,selected_biosample=selected['biosample'],
                 selected_bioproject=selected['bioproject'],selected_wgs_master=selected['wgs_master'],
                 selected_infraspecific_name=selected['infraspecific_name'],
                 source_excel_row='',source_strain='',source_biosample='',source_bioproject='',source_wgs_master='',
                 biosample_equal='',wgs_project_equal_ignoring_version='',bioproject_equal='',
                 status='no_unique_published_sample_identifier',scope='Same sample/project identifiers do not prove identical assembly sequence, annotation, taxonomic concept or experimental phenotype validation.')
        if len(entries)==1:
            number,record=entries[0]
            row.update(source_excel_row=number,source_strain=str(record[1]).strip(),
                       source_biosample=str(record[10]).strip(),source_bioproject=str(record[9]).strip(),
                       source_wgs_master=str(record[11]).strip())
            row['biosample_equal']=row['source_biosample']==row['selected_biosample']
            row['wgs_project_equal_ignoring_version']=row['source_wgs_master'].split('.')[0]==row['selected_wgs_master'].split('.')[0]
            row['bioproject_equal']=row['source_bioproject']==row['selected_bioproject']
            row['status']='same_biosample_and_wgs_project' if row['biosample_equal'] and row['wgs_project_equal_ignoring_version'] else 'identifier_disagreement_requires_review'
        results.append(row)
    path=ROOT/'metadata/miyauchi_ecology_sample_identity.tsv'
    with path.open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(results[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(results)
    receipt=dict(status='complete_published_sample_identifier_linkage',records=len(results),
                 dispositions=dict(Counter(r['status'] for r in results)),
                 inputs={k:sha(v) for k,v in files.items()},script_sha256=sha(Path(__file__)),
                 output_sha256=sha(path),scope='Sample and WGS project linkage only. Selected assembly accession matched exactly; WGS version suffix ignored explicitly. Unmatched entries remain unresolved, not negative matches.')
    (ROOT/'metadata/miyauchi_ecology_sample_identity_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':
    main()
