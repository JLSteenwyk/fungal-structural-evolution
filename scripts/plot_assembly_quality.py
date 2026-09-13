#!/usr/bin/env python3
"""Compare contiguity and marker recovery while preserving assembly metric scopes."""
import argparse
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from retrieve_assembly_statistics import ROOT, sha


def read_table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new immutable quality snapshot')
    ncbi_path = ROOT / 'metadata/assembly_quality_metrics.tsv'
    ext_path = ROOT / 'metadata/external_assembly_quality.tsv'
    nr = json.loads((ROOT / 'metadata/assembly_quality_receipt.json').read_text())
    er = json.loads((ROOT / 'metadata/external_assembly_quality_receipt.json').read_text())
    if sha(ncbi_path) != nr['metrics_sha256'] or sha(ext_path) != er['table_sha256']:
        raise ValueError('Changed assembly-quality metrics')
    if nr['ncbi_failed'] or nr['ncbi_verified'] != nr['ncbi_requested']:
        raise ValueError('Complete NCBI statistics collection required for this figure')
    ncbi = read_table(ncbi_path)
    external = {r['taxon_id']: r for r in read_table(ext_path)}
    busco_path = ROOT / 'metadata/busco_eukaryota_qc.tsv'
    busco = {r['taxon_id']: r for r in read_table(busco_path)}
    for r in busco.values():
        path = ROOT / r['summary_path']
        if sha(path) != r['summary_sha256']:
            raise ValueError('Changed BUSCO summary')
        result = json.loads(path.read_text())['results']
        if float(r['complete_percent']) != result['Complete percentage'] or int(r['markers']) != result['n_markers']:
            raise ValueError('BUSCO table differs from source summary')
    joined = []
    for row in ncbi:
        ext = external.get(row['taxon_id'])
        complete = float(busco[row['taxon_id']]['complete_percent'])
        joined.append({k: row[k] for k in ['taxon_id', 'species_name', 'study_role', 'assembly_accession']} | {
            'assembly_source': 'external_FASTA' if ext else 'NCBI_report',
            'ncbi_all_contig_N50': row['all_contig-N50'],
            'external_record_N50': ext['record_N50'] if ext else '',
            'busco_complete_percent': complete, 'busco_marker_panel': 'eukaryota_odb12.2_125_markers'})
    if set(busco) != {r['taxon_id'] for r in joined} or len(external) != er['taxa']:
        raise ValueError('Quality and taxon sets differ')
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.7))
    counts = {}
    for role, color, name in [('ingroup', '#39877a', 'Fungal entries'), ('outgroup', '#a05d3d', 'Outgroups')]:
        rows = [r for r in joined if r['assembly_source'] == 'NCBI_report' and r['study_role'] == role
                and r['ncbi_all_contig_N50'] != '' and float(r['ncbi_all_contig_N50']) > 0]
        counts[role] = len(rows)
        axes[0].scatter([float(r['ncbi_all_contig_N50']) for r in rows], [r['busco_complete_percent'] for r in rows],
                        color=color, s=15, alpha=.6, linewidths=0, label=f'{name} (n={len(rows)})')
    axes[0].legend(frameon=False, fontsize=8, loc='lower right')
    ext_rows = [r for r in joined if r['assembly_source'] == 'external_FASTA']
    for r in ext_rows:
        axes[1].scatter(float(r['external_record_N50']), r['busco_complete_percent'],
                        color='#39877a' if r['study_role'] == 'ingroup' else '#a05d3d', s=35)
        if r['species_name'] in ['Pirum gemmata', 'Abeoforma whisleri']:
            offset = (8, 10) if r['species_name'].startswith('Abeoforma') else (8, -15)
            axes[1].annotate(r['species_name'], (float(r['external_record_N50']), r['busco_complete_percent']),
                             xytext=offset, textcoords='offset points', fontsize=8)
    for ax in axes:
        ax.set_xscale('log')
        ax.set_ylim(-3, 103)
        ax.set_ylabel('Complete broad BUSCO markers (%)')
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].set(xlabel='NCBI whole-assembly contig N50 (bp)', title='Deposited assembly statistics')
    axes[1].set(xlabel='Deposited FASTA-record N50 (bp)', title=f'External genomes (n={len(ext_rows)})')
    fig.suptitle('Assembly contiguity and broad marker recovery')
    fig.text(.5, .02, 'Metric scopes differ between panels. Neither contiguity nor marker recovery establishes contamination status.',
             ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .06, 1, .95))
    args.output.mkdir(parents=True)
    for extension in ['svg', 'png', 'pdf']:
        fig.savefig(args.output / f'assembly_quality.{extension}', dpi=200)
    plt.close(fig)
    table = args.output / 'assembly_and_marker_quality.tsv'
    with table.open('w') as handle:
        writer = csv.DictWriter(handle, list(joined[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(joined)
    receipt = {'ncbi_quality_receipt_sha256': sha(ROOT / 'metadata/assembly_quality_receipt.json'),
        'external_quality_receipt_sha256': sha(ROOT / 'metadata/external_assembly_quality_receipt.json'),
        'busco_table_sha256': sha(busco_path), 'script_sha256': sha(Path(__file__)), 'taxa': len(joined),
        'ncbi_taxa_plotted': counts, 'external_taxa_plotted': len(ext_rows),
        'artifacts': {p.name: sha(p) for p in args.output.iterdir()}}
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
