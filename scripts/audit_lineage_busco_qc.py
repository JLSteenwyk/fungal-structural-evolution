#!/usr/bin/env python3
"""Validate every lineage BUSCO marker call and join panel-specific quality summaries."""
import argparse
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT, read_table, sha
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('run', 'datasets', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    receipt = checked_receipt(a.run)
    config = json.loads((a.run / 'config.json').read_text())
    if receipt['status'] != 'complete_full_ingroup_lineage_qc' or receipt['config_sha256'] != sha(a.run / 'config.json'):
        raise ValueError('Incomplete or changed run')
    for field, path in [('manifest_sha256', ROOT / 'metadata/analysis_manifest.tsv'),
                        ('input_receipts_sha256', ROOT / 'metadata/qc_input_receipts.json'),
                        ('dataset_receipt_sha256', a.datasets / 'receipt.json')]:
        if config[field] != sha(path):
            raise ValueError('Changed source')
    manifest = {r['taxon_id']: r for r in read_table(ROOT / 'metadata/analysis_manifest.tsv')}
    inputs = {r['taxon_id']: r for r in json.loads((ROOT / 'metadata/qc_input_receipts.json').read_text())}
    datasets = {r['dataset']: r for r in json.loads((a.datasets / 'receipt.json').read_text())['datasets']}
    universes = {}
    for name, ds in datasets.items():
        folder = ROOT / ds['path']
        for path, digest in ds['files_sha256'].items():
            if sha(folder / path) != digest:
                raise ValueError('Changed dataset artifact')
        universes[name] = {line.split()[0] for line in (folder / 'scores_cutoff').read_text().splitlines() if line.strip()}
        if len(universes[name]) != int(ds['config']['number_of_BUSCOs']):
            raise ValueError('Dataset marker universe differs')
    assignments = read_table(a.run / 'taxon_dataset_assignments.tsv')
    broad = {r['taxon_id']: r for r in read_table(ROOT / 'metadata/busco_eukaryota_qc.tsv')}
    if len(assignments) != len(manifest) or {r['taxon_id'] for r in assignments} != manifest.keys() or broad.keys() != manifest.keys() or inputs.keys() != manifest.keys():
        raise ValueError('Taxon universe differs')
    expected_jobs = {f"jobs/{r['taxon_id']}.json" for r in assignments if r['dataset']}
    if expected_jobs != {k for k in receipt['artifacts'] if k.startswith('jobs/')}:
        raise ValueError('Job grid differs')
    rows, pins, dataset_values = [], {}, defaultdict(list)
    total_calls = total_hits = 0
    categories = {'Complete': 'Single copy BUSCOs', 'Duplicated': 'Multi copy BUSCOs',
                  'Fragmented': 'Fragmented BUSCOs', 'Missing': 'Missing BUSCOs'}
    for assignment in assignments:
        taxon = assignment['taxon_id']; source = inputs[taxon]; m = manifest[taxon]
        expected_dataset = config['selection'].get(m['lineage'].split(';')[0], config['fallback']) if m['study_role'] == 'ingroup' else ''
        if assignment['dataset'] != expected_dataset or assignment['species_name'] != m['species_name'] or assignment['study_role'] != m['study_role']:
            raise ValueError('Dataset assignment differs')
        b = broad[taxon]; bp = ROOT / b['summary_path']
        if sha(bp) != b['summary_sha256'] or b['input_sha256'] != source['sha256']:
            raise ValueError('Broad QC source differs')
        br = json.loads(bp.read_text())
        if br['lineage_dataset']['name'] != 'eukaryota_odb12.2' or int(br['results']['n_markers']) != int(b['markers']):
            raise ValueError('Broad denominator differs')
        for table_field, raw_field in [('complete_percent', 'Complete percentage'), ('single_copy_percent', 'Single copy percentage'),
                                       ('duplicated_percent', 'Multi copy percentage'), ('fragmented_percent', 'Fragmented percentage'), ('missing_percent', 'Missing percentage')]:
            if float(b[table_field]) != float(br['results'][raw_field]):
                raise ValueError('Broad percentage differs')
        row = dict(assignment)
        row.update({'broad_' + k: b[k] for k in ('markers', 'complete_percent', 'single_copy_percent', 'duplicated_percent', 'fragmented_percent', 'missing_percent')})
        row['lineage_qc_status'] = 'not_applied_outgroup'
        for field in ('markers', 'complete', 'single_copy', 'duplicated', 'fragmented', 'missing'):
            row['lineage_' + field] = ''
        for field in ('complete', 'single_copy', 'duplicated', 'fragmented', 'missing'):
            row['lineage_' + field + '_percent'] = ''
        if expected_dataset:
            job = json.loads((a.run / 'jobs' / (taxon + '.json')).read_text())
            if job['status'] != 'complete_lineage_protein_qc' or job['returncode'] != 0 or job['config_sha256'] != receipt['config_sha256'] or job['input_sha256'] != source['sha256'] or job['dataset'] != expected_dataset:
                raise ValueError('Invalid taxon job')
            sp = ROOT / job['summary_path']; summary = json.loads(sp.read_text())
            if sha(sp) != job['summary_sha256'] or summary['lineage_dataset']['name'] != expected_dataset or summary['lineage_dataset']['creation_date'] != datasets[expected_dataset]['config']['creation_date']:
                raise ValueError('Changed summary or dataset identity')
            fasta = ROOT / source['input_path']
            if sha(fasta) != source['sha256']:
                raise ValueError('Changed proteome')
            with fasta.open() as handle:
                ids = [title.split()[0] for title, seq in SimpleFastaParser(handle)]
            if len(ids) != len(set(ids)):
                raise ValueError('Duplicate protein identifier')
            ids = set(ids)
            full = a.run / taxon / ('run_' + expected_dataset) / 'full_table.tsv'
            calls = defaultdict(list)
            for line in full.read_text().splitlines():
                if not line or line.startswith('#'):
                    continue
                values = line.split('\t'); marker, status = values[:2]
                if status not in categories:
                    raise ValueError('Unexpected BUSCO status')
                if status != 'Missing':
                    if len(values) < 5 or values[2] not in ids or not math.isfinite(float(values[3])) or int(values[4]) <= 0:
                        raise ValueError('Invalid BUSCO protein hit')
                    total_hits += 1
                calls[marker].append(values)
            if calls.keys() != universes[expected_dataset]:
                raise ValueError('Full marker grid differs')
            counts = Counter()
            for marker, hits in calls.items():
                status = hits[0][1]
                if any(h[1] != status for h in hits) or (status != 'Duplicated' and len(hits) != 1) or (status == 'Duplicated' and (len(hits) < 2 or len({h[2] for h in hits}) != len(hits))):
                    raise ValueError('Inconsistent marker multiplicity')
                counts[categories[status]] += 1
            n = len(calls); total_calls += n
            expected_counts = {'n_markers': n, 'Complete BUSCOs': counts['Single copy BUSCOs'] + counts['Multi copy BUSCOs'],
                               **{field: counts[field] for field in categories.values()}}
            if job['counts'] != expected_counts or any(int(summary['results'][k]) != v for k, v in expected_counts.items()):
                raise ValueError('Raw marker counts differ from summaries')
            row['lineage_qc_status'] = 'audited_lineage_panel'; row['lineage_markers'] = n
            for short, field in [('complete', 'Complete BUSCOs'), ('single_copy', 'Single copy BUSCOs'), ('duplicated', 'Multi copy BUSCOs'), ('fragmented', 'Fragmented BUSCOs'), ('missing', 'Missing BUSCOs')]:
                row['lineage_' + short] = expected_counts[field]
                row['lineage_' + short + '_percent'] = 100 * expected_counts[field] / n
            dataset_values[expected_dataset].append(row)
            pins[str(full)] = sha(full); pins[str(sp)] = sha(sp)
        rows.append(row)
    if sum(len(v) for v in dataset_values.values()) != receipt['completed'] or len(rows) - receipt['completed'] != receipt['outgroups_retaining_broad_qc']:
        raise ValueError('Completion accounting differs')
    summaries = []
    for dataset, values in sorted(dataset_values.items()):
        summaries.append({'dataset': dataset, 'taxa': len(values), 'markers': values[0]['lineage_markers'],
                          **{field + '_percent_median': statistics.median(r['lineage_' + field + '_percent'] for r in values)
                             for field in ('complete', 'single_copy', 'duplicated', 'fragmented', 'missing')}})
    a.output.mkdir(parents=True)
    write_table(a.output / 'taxon_quality.tsv', rows); write_table(a.output / 'dataset_summary.tsv', summaries)
    (a.output / 'raw_output_hashes.json').write_text(json.dumps(pins, indent=2) + '\n')
    result = {'status': 'passed_full_lineage_busco_readback', 'taxa': len(rows), 'lineage_audited_taxa': receipt['completed'],
              'outgroups_broad_only': receipt['outgroups_retaining_broad_qc'], 'lineage_marker_calls': total_calls, 'lineage_protein_hit_rows': total_hits,
              'source_receipt_sha256': sha(a.run / 'receipt.json'), 'script_sha256': sha(Path(__file__)),
              'interpretation': 'Every lineage full-table marker identity/status/count and hit protein identity checked, with broad QC retained for all taxa. Dataset-specific denominators are not interchangeable. Protein-mode completeness cannot establish assembly completeness, contamination, ploidy or biological loss; no taxon filtering.',
              'artifacts': {path.name: sha(path) for path in a.output.iterdir()}}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
