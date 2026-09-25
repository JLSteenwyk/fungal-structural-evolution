#!/usr/bin/env python3
"""Union disjoint AFDB mappings and unchanged confidence-qualified encodings."""
import argparse
import csv
import gzip
import json
import shutil
from collections import Counter
from pathlib import Path
from compare_marker_structures import sha
from assess_pae_sensitivity import checked_receipt


def read_table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_table(path, rows):
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    a = p.parse_args(); plan = json.loads(a.plan.read_text())
    pins = {str(a.plan): sha(a.plan), **plan['pins']}
    def verify():
        for path, digest in pins.items():
            if sha(Path(path)) != digest:
                raise ValueError('Changed source: ' + path)
    verify()
    output = Path(plan['output'])
    if output.exists():
        raise FileExistsError('Use a new immutable union')
    if shutil.disk_usage(output.parent).free < plan['minimum_free_disk_gib'] * 2**30:
        raise ValueError('Insufficient disk')
    links = {}; models = {}; sequences = set(); encodings = {}; sources = []; totals = Counter()
    matrix_hash = None; policy = None; coverage = {}; cohort_links = []
    for cohort in plan['cohorts']:
        mapping, encoded = Path(cohort['mapping']), Path(cohort['encodings'])
        mr, er = checked_receipt(mapping), checked_receipt(encoded)
        for root, receipt in [(mapping, mr), (encoded, er)]:
            pins[str(root / 'receipt.json')] = sha(root / 'receipt.json')
            pins.update({str(root / name): digest for name, digest in receipt['artifacts'].items()})
        if mr['source_policy']['provider'] != 'GDM' or mr['source_policy']['tool'] != 'AlphaFold Monomer v2.0 pipeline':
            raise ValueError('Expected AFDB GDM monomers')
        if matrix_hash is None:
            matrix_hash, policy = mr['matrix_receipt_sha256'], mr['source_policy']
        if matrix_hash != mr['matrix_receipt_sha256'] or policy != mr['source_policy']:
            raise ValueError('Different matrix or source-selection policy')
        if er['mapping_receipt_sha256'] != sha(mapping / 'receipt.json') or er['confidence_stage'] != 'plddt_and_pae':
            raise ValueError('Unbound or unqualified encodings')
        local_models = {}
        for model in json.loads((mapping / 'model_provenance.json').read_text()):
            name = Path(model['path']).stem
            if name in models or model['sequence_sha256'] in sequences:
                raise ValueError('Overlapping model or sequence')
            models[name] = model; local_models[name] = model; sequences.add(model['sequence_sha256'])
        local_encodings = {}; local_totals = Counter()
        for row in read_table(encoded / 'model_summary.tsv'):
            name = row['model_name']
            if name in encodings or name in local_encodings or name not in local_models:
                raise ValueError('Unexpected encoding identity')
            if row['sequence_sha256'] != local_models[name]['sequence_sha256'] or sha(Path(row['encoding_path'])) != row['encoding_sha256']:
                raise ValueError('Changed encoding')
            local_encodings[name] = row; encodings[name] = row
            for key in er['totals']:
                local_totals[key] += int(row[key])
        if set(local_encodings) != set(local_models) or len(local_models) != er['models'] or dict(local_totals) != er['totals']:
            raise ValueError('Encoding scope or totals differ')
        totals.update(local_totals)
        local_links = {}
        for row in read_table(mapping / 'marker_structure_links.tsv'):
            key = row['marker'], row['taxon_id']
            if key in links or key in local_links:
                raise ValueError('Overlapping marker/taxon')
            model = local_models[Path(row['model_path']).stem]
            if row['sequence_sha256'] != model['sequence_sha256'] or row['model_sha256'] != model['sha256']:
                raise ValueError('Link/model mismatch')
            links[key] = row; local_links[key] = row
        if len(local_links) != mr['marker_proteins_linked'] or len(local_models) != mr['distinct_models']:
            raise ValueError('Mapping scope differs')
        cohort_links.append((mapping, local_links, mr['matrix_residue_links']))
        for row in read_table(mapping / 'marker_coverage.tsv'):
            marker = row['marker']; count = int(row['single_copy_taxa'])
            if marker in coverage and coverage[marker] != count:
                raise ValueError('Different source marker universe')
            coverage[marker] = count
        sources.append(dict(cohort, mapping_receipt_sha256=sha(mapping / 'receipt.json'), encoding_receipt_sha256=sha(encoded / 'receipt.json')))
    if len(models) != plan['models'] or len(links) != plan['marker_links']:
        raise ValueError('Planned union dimensions differ')
    output.mkdir(); mapping_out = output / 'mapping'; mapping_out.mkdir()
    encoding_out = output / 'encodings'; encoding_out.mkdir()
    residue_total = 0
    with gzip.open(mapping_out / 'matrix_to_structure_residues.tsv.gz', 'wt') as target:
        writer = None; header = None
        for mapping, local, expected_count in cohort_links:
            last = {}; counts = Counter(); count = 0
            with gzip.open(mapping / 'matrix_to_structure_residues.tsv.gz', 'rt') as handle:
                reader = csv.DictReader(handle, delimiter='\t')
                if writer is None:
                    header = reader.fieldnames
                    writer = csv.DictWriter(target, header, delimiter='\t', lineterminator='\n'); writer.writeheader()
                if reader.fieldnames != header:
                    raise ValueError('Residue schema differs')
                for row in reader:
                    key = row['marker'], row['taxon_id']; link = local[key]
                    column = int(row['matrix_column_1based'])
                    if column <= last.get(key, 0) or any(row[k] != link[k] for k in ['protein_id', 'model_id', 'model_version']):
                        raise ValueError('Duplicate, unordered, or mismatched residue link')
                    last[key] = column; counts[key] += 1; count += 1; writer.writerow(row)
            if count != expected_count or any(counts[k] != int(v['retained_marker_residues']) for k, v in local.items()):
                raise ValueError('Incomplete source residue grid')
            residue_total += count
            print('Copied', mapping, count, 'residue links', flush=True)
    if residue_total != plan['residue_links']:
        raise ValueError('Union residue total differs')
    write_table(mapping_out / 'marker_structure_links.tsv', [links[k] for k in sorted(links)])
    counts = Counter(k[0] for k in links)
    write_table(mapping_out / 'marker_coverage.tsv', [dict(marker=k, single_copy_taxa=coverage[k], taxa_with_exact_sequence_structure=counts[k]) for k in sorted(coverage)])
    (mapping_out / 'model_provenance.json').write_text(json.dumps([models[k] for k in sorted(models)], indent=2) + '\n')
    write_table(encoding_out / 'model_summary.tsv', [encodings[k] for k in sorted(encodings)])
    verify()
    shared = dict(source_cohorts=sources, plan_sha256=sha(a.plan), script_sha256=sha(Path(__file__)))
    mr = dict(shared, status='complete_disjoint_afdb_mapping_union', source_policy=policy,
              matrix_receipt_sha256=matrix_hash, distinct_models=len(models), marker_proteins_linked=len(links),
              distinct_taxa_linked=len({k[1] for k in links}), matrix_residue_links=residue_total,
              artifacts={f.name: sha(f) for f in mapping_out.iterdir()})
    (mapping_out / 'receipt.json').write_text(json.dumps(mr, indent=2) + '\n')
    er = dict(shared, status='complete_disjoint_afdb_qualified_encoding_union', models=len(encodings), totals=dict(totals),
              confidence_stage='plddt_and_pae', mapping_receipt_sha256=sha(mapping_out / 'receipt.json'),
              artifacts={'model_summary.tsv': sha(encoding_out / 'model_summary.tsv')},
              scope='Unchanged source arrays, rehashed at union; source qualification receipts retained. No new structure inference, confidence recalibration, or fitted evolutionary results.')
    (encoding_out / 'receipt.json').write_text(json.dumps(er, indent=2) + '\n')
    print(json.dumps(dict(models=len(models), marker_links=len(links), residue_links=residue_total)))


if __name__ == '__main__':
    main()
