#!/usr/bin/env python3
"""Prepare full-matrix mixture sensitivity resources from current guide log estimates."""
import argparse
import hashlib
import json
from pathlib import Path
import re
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use a new immutable resource snapshot')
    matrices = []
    for label, matrix_name, guide_name in [
        ('profile', 'profile-matrix-50-v1', 'initial-guide-v1'),
        ('mafft', 'mafft-matrix-50-v1', 'initial-mafft-guide-v1')]:
        folder = ROOT / 'results/phylogeny' / matrix_name
        receipt_path = folder / 'receipt.json'; receipt = json.loads(receipt_path.read_text())
        for name, digest in receipt['artifacts'].items():
            if sha(folder / name) != digest:
                raise ValueError('Changed matrix artifact: ' + name)
        with (folder / 'matrix.faa').open() as handle:
            records = list(SeqIO.parse(handle, 'fasta'))
        if len(records) != 526 or len({r.id for r in records}) != 526:
            raise ValueError('Full taxon identity grid differs')
        if any(len(r.seq) != receipt['columns'] for r in records):
            raise ValueError('Matrix length differs')
        log_path = ROOT / 'results/phylogeny' / guide_name / 'guide.log'
        # The producer is live. Freeze one read, not a hash of a later mutable log.
        raw = log_path.read_bytes(); text = raw.decode()
        counts = re.search(r'Alignment has (\d+) sequences with (\d+) columns, (\d+) distinct patterns', text)
        mem = re.search(r'NOTE: (\d+) MB RAM .* is required!', text)
        if not counts or not mem or int(counts[1]) != 526 or int(counts[2]) != receipt['columns']:
            raise ValueError('Guide startup estimates do not match the matrix')
        base_mb = int(mem[1])
        matrices.append({'label': label, 'matrix': str(folder.relative_to(ROOT)),
            'matrix_receipt_sha256': sha(receipt_path), 'matrix_sha256': sha(folder / 'matrix.faa'),
            'taxa': 526, 'columns': receipt['columns'], 'patterns': int(counts[3]),
            'guide_log_path': str(log_path.relative_to(ROOT)),
            'guide_log_snapshot_sha256': hashlib.sha256(raw).hexdigest(),
            'guide_memory_report_line': mem[0], 'homogeneous_reported_mb': base_mb,
            'mixture_memory_scenarios_reported_mb_units': [
                {'components_with_empirical_profile': k+1, 'model': f'LG+C{k}+F+G4',
                 'linear_component_scenario_mb': base_mb*(k+1),
                 'with_25_percent_headroom_mb': base_mb*(k+1)*1.25} for k in [20,60]],
            'guide_log_snapshot': raw})
    expected = None
    for row in matrices:
        with (ROOT / row['matrix'] / 'matrix.faa').open() as handle:
            ids = {r.id for r in SeqIO.parse(handle, 'fasta')}
        if expected is not None and ids != expected:
            raise ValueError('Guide sensitivity requires identical taxon universes')
        expected = ids
    a.output.mkdir(parents=True)
    for row in matrices:
        name = row['label'] + '_guide.log'
        (a.output / name).write_bytes(row.pop('guide_log_snapshot'))
    result = {'status': 'prepared_full_matrix_mixture_resource_scenarios',
        'script_sha256': sha(Path(__file__)), 'matrices': matrices,
        'source': 'https://iqtree.github.io/doc/Complex-Models#site-specific-frequency-models',
        'source_method_doi': '10.1093/sysbio/syx068',
        'interpretation': 'Linear component scaling of IQ-TREE startup RAM estimates, with explicit 25% scenario headroom; neither measured mixture memory nor an upper bound. PMSF profile estimation retains mixture memory costs. C20 is the initial feasible full-matrix sensitivity; C60 needs a separately verified memory-saving configuration. Neither model addresses all across-lineage compositional variation.',
        'proposed_design': {'matrix_labels': ['profile','mafft'], 'guide_labels': ['profile','mafft'],
            'crossed_runs': 4, 'model': 'LG+C20+F+G4', 'method': 'PMSF',
            'concurrent_runs': 1, 'threads': 16, 'iqtree_memory_limit': '600G',
            'minimum_available_memory_gib': 750, 'output_allowance_gb': 100,
            'planning_wall_hours_per_run': [24,336],
            'runtime_basis': 'Conservative planning envelope for full-matrix mixture profile estimation and supported tree search, not a measured prediction.',
            'prerequisites': 'Both completed guide receipts, exact tip grids and finite branch lengths; freeze each guide before use. All four guide/alignment combinations preserve 526 taxa. Support, model adequacy, gene discordance and taxon/marker sensitivity remain required.'},
        'artifacts': {p.name: sha(p) for p in a.output.iterdir()}}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({r['label']:r['mixture_memory_scenarios_reported_mb_units'] for r in matrices},indent=2))


if __name__ == '__main__':
    main()
