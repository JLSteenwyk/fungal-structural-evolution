#!/usr/bin/env python3
"""Import source-ranked ecological candidates without converting genera to species evidence."""
import csv
import datetime
import hashlib
import json
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path
import openpyxl

ROOT = Path(__file__).resolve().parents[1]
URL = 'https://static-content.springer.com/esm/art%3A10.1007%2Fs13225-020-00466-2/MediaObjects/13225_2020_466_MOESM4_ESM.xlsx'
EXPECTED = '8f8f13b77fb50b5324ea099c876feb749529cb0f1a6ae86de0700f4430de5e44'
DOI = 'https://doi.org/10.1007/s13225-020-00466-2'
FIELDS = {'primary_lifestyle': 'primary_lifestyle', 'secondary_lifestyle': 'Secondary_lifestyle',
          'endophytic_capability': 'Endophytic_interaction_capability_template',
          'plant_pathogenic_capacity': 'Plant_pathogenic_capacity_template',
          'decay_substrate': 'Decay_substrate_template', 'decay_type': 'Decay_type_template',
          'aquatic_habitat': 'Aquatic_habitat_template', 'animal_biotrophic_capacity': 'Animal_biotrophic_capacity_template',
          'growth_form': 'Growth_form_template', 'fruitbody_type': 'Fruitbody_type_template'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def match_genus(row, genera):
    if row['study_role'] != 'ingroup':
        return 'outgroup_not_assigned', []
    genus = row['species_name'].split()[0]
    matches = genera.get(genus, [])
    return ('no_exact_genus_match' if not matches else
            'ambiguous_genus_rows' if len(matches) > 1 else 'genus_candidate_requires_species_verification'), matches


def main():
    path = ROOT / 'data/traits/fungaltraits_genera_publisher.xlsx'
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        with urllib.request.urlopen(URL, timeout=90) as response:
            data = response.read()
        if hashlib.sha256(data).hexdigest() != EXPECTED:
            raise ValueError('Publisher supplement differs from the pinned snapshot')
        temp = path.with_suffix('.partial')
        temp.write_bytes(data)
        temp.replace(path)
    if sha(path) != EXPECTED:
        raise ValueError('Changed trait source workbook')
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    values = iter(workbook['data'].iter_rows(values_only=True))
    header = next(values)
    if not {'GENUS', 'jrk_template', *FIELDS.values()} <= set(header):
        raise ValueError('Unexpected trait schema')
    genera, source_rows = defaultdict(list), []
    for excel_row, values in enumerate(values, 2):
        row = {key: ('' if value is None else str(value).strip()) for key, value in zip(header, values) if key}
        if not row.get('GENUS'):
            continue
        row['excel_row'] = excel_row
        genera[row['GENUS']].append(row)
        source_rows.append(row)
    workbook.close()
    manifest = ROOT / 'metadata/analysis_manifest.tsv'
    review_path = ROOT / 'metadata/taxon_label_review.tsv'
    with manifest.open() as handle:
        taxa = list(csv.DictReader(handle, delimiter='\t'))
    with review_path.open() as handle:
        identity_flags = {r['taxon_id']: r['flags'] for r in csv.DictReader(handle, delimiter='\t')}
    rows = []
    for taxon in taxa:
        status, matches = match_genus(taxon, genera)
        candidate = matches[0] if len(matches) == 1 else {}
        row = {k: taxon[k] for k in ['taxon_id', 'species_name', 'study_role', 'assembly_accession']}
        row.update(match_status=status, queried_genus=taxon['species_name'].split()[0],
                   source_genus=candidate.get('GENUS', ''), source_rank='genus' if matches else '',
                   source_record_id=candidate.get('jrk_template', ''),
                   source_excel_rows=';'.join(str(r['excel_row']) for r in matches),
                   source_doi=DOI if matches else '', source_phylum=candidate.get('Phylum', ''),
                   source_class=candidate.get('Class', ''), source_family=candidate.get('Family', ''),
                   taxon_identity_flags=identity_flags.get(taxon['taxon_id'], ''),
                   source_genus_comment_present=bool(candidate.get('COMMENT on genus', '')),
                   source_lifestyle_comment_present=bool(candidate.get('Comment_on_lifestyle_template', '')),
                   species_verified=False, eligible_for_confirmatory_tests=False)
        for field, source_field in FIELDS.items():
            row['candidate_' + field] = candidate.get(source_field, '')
        rows.append(row)
    table = ROOT / 'metadata/ecology_candidates.tsv'
    with table.open('w') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    candidates = [r for r in rows if r['match_status'] == 'genus_candidate_requires_species_verification']
    receipt = {'source_url': URL, 'source_doi': DOI, 'source_path': str(path.relative_to(ROOT)),
        'source_sha256': EXPECTED, 'source_bytes': path.stat().st_size,
        'source_sheet': 'data', 'source_genus_rows': len(source_rows), 'source_distinct_genus_labels': len(genera),
        'duplicate_source_genera': {k: len(v) for k, v in genera.items() if len(v) > 1},
        'snapshot_recorded_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'manifest_sha256': sha(manifest), 'taxon_review_sha256': sha(review_path), 'taxa': len(rows),
        'match_status_counts': dict(Counter(r['match_status'] for r in rows)),
        'candidate_primary_lifestyle_counts': dict(Counter(r['candidate_primary_lifestyle'] or 'not_recorded' for r in candidates)),
        'species_verified': 0, 'eligible_for_confirmatory_tests': 0,
        'openpyxl_version': openpyxl.__version__, 'script_sha256': sha(Path(__file__)), 'table_sha256': sha(table),
        'interpretation': 'Genus-level candidate evidence only. No synonym resolution, species-level inheritance, negative-state imputation or replicated-transition inference. Supplement contents are pinned as retrieved, not asserted to reproduce historical paper counts.'}
    (ROOT / 'metadata/ecology_candidates_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
