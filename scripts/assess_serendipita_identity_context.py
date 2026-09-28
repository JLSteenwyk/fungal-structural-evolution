#!/usr/bin/env python3
"""Describe all selected Serendipita pairs in both frozen marker matrices."""
import csv
import hashlib
import itertools
import json
from pathlib import Path
from Bio import SeqIO


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    manifest = Path('metadata/analysis_manifest.tsv')
    with manifest.open() as stream:
        taxa = {r['taxon_id']: r for r in csv.DictReader(stream, delimiter='\t')
                if r['species_name'].startswith('Serendipita ')}
    assert len(taxa) == 5
    pins = {str(manifest): sha(manifest)}
    rows = []
    canonical = set('ARNDCQEGHILKMFPSTWYV')
    for method in ['profile', 'mafft']:
        root = Path('results/phylogeny') / (method + '-matrix-50-v1')
        rp = root / 'receipt.json'
        receipt = json.loads(rp.read_text())
        p = root / 'matrix.faa'
        assert sha(p) == receipt['artifacts']['matrix.faa']
        assert receipt['manifest_sha256'] == sha(manifest)
        pins.update({str(x): sha(x) for x in [rp, p]})
        seqs = {r.id: str(r.seq) for r in SeqIO.parse(p, 'fasta') if r.id in taxa}
        assert set(seqs) == set(taxa)
        for a, b in itertools.combinations(sorted(taxa), 2):
            assert len(seqs[a]) == len(seqs[b]) == receipt['columns']
            pairs = [(x, y) for x, y in zip(seqs[a], seqs[b]) if x in canonical and y in canonical]
            matches = sum(x == y for x, y in pairs)
            # Independent position-set and mismatch accounting.
            overlap = {i for i, aa in enumerate(seqs[a]) if aa in canonical} & {i for i, aa in enumerate(seqs[b]) if aa in canonical}
            differences = sum(seqs[a][i] != seqs[b][i] for i in overlap)
            assert len(overlap) == len(pairs) and differences + matches == len(pairs)
            rows.append(dict(method=method, taxon_a=a, name_a=taxa[a]['species_name'],
                taxon_b=b, name_b=taxa[b]['species_name'], matrix_columns=receipt['columns'],
                paired_canonical_columns=len(pairs), identical_columns=matches,
                differing_columns=differences, excluded_columns=receipt['columns']-len(pairs),
                observed_difference_fraction=differences/len(pairs) if pairs else 'NA'))
    out = Path('results/phylogeny/serendipita-identity-context-20260927-v1')
    out.mkdir(exist_ok=False)
    p = out / 'all_pair_observed_differences.tsv'
    with p.open('w') as stream:
        writer = csv.DictWriter(stream, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    result = dict(status='complete_selected_serendipita_pairwise_context', taxa=5, comparisons=20,
        pins=pins, script_sha256=sha(__file__), artifacts={p.name:sha(p)},
        scope='All pairs among the five selected Serendipita entries, both frozen marker matrices. Canonical residues only; gaps and unknowns excluded explicitly. Raw observed amino-acid differences, not corrected evolutionary distance, independent sites, species delimitation, contamination assessment or evidence that three unnamed isolates are distinct species.')
    rp = out / 'receipt.json'
    rp.write_text(json.dumps(result, indent=2) + '\n')
    result.update(completed_receipt_path=str(rp), completed_receipt_sha256=sha(rp))
    Path('metadata/serendipita_identity_context_completed_20260927.json').write_text(json.dumps(result, indent=2)+'\n')
    for row in rows:
        print(row['method'], row['taxon_a'], row['taxon_b'], row['paired_canonical_columns'], round(row['observed_difference_fraction'],6))


if __name__ == '__main__':
    main()
