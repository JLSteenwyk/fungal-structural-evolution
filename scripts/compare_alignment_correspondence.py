#!/usr/bin/env python3
"""Compare residue homology edges between profile and de novo alignments."""
import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from Bio import AlignIO, SeqIO
from build_species_matrix import AA, sha
from map_marker_structures import residue_positions


def choose2(n):
    return n * (n - 1) // 2


def agreement(profile, mafft):
    """Count identical cross-taxon residue-pair edges without expanding all pairs."""
    common = profile.keys() & mafft.keys()
    pcounts, mcounts = Counter(), Counter()
    cells = Counter()
    for residue in common:
        p, m = profile[residue], mafft[residue]
        pcounts[p] += 1
        mcounts[m] += 1
        cells[(p, m)] += 1
    pedges = sum(choose2(n) for n in pcounts.values())
    medges = sum(choose2(n) for n in mcounts.values())
    same = sum(choose2(n) for n in cells.values())
    union = pedges + medges - same
    result = {'profile_retained_residues': len(profile), 'mafft_retained_residues': len(mafft),
              'common_retained_residues': len(common), 'profile_edges_on_common_residues': pedges,
              'mafft_edges_on_common_residues': medges, 'shared_edges': same,
              'edge_jaccard': same / union if union else '',
              'fraction_profile_edges_shared': same / pedges if pedges else '',
              'fraction_mafft_edges_shared': same / medges if medges else ''}
    all_profile = Counter(profile.values())
    by_column = defaultdict(int)
    for (p, _), n in cells.items():
        by_column[p] += choose2(n)
    columns = [{'profile_column_1based': p, 'profile_retained_residues': total,
                'common_retained_residues': pcounts[p],
                'fraction_residues_retained_by_both': pcounts[p] / total,
                'profile_edges_on_common_residues': choose2(pcounts[p]), 'shared_edges': by_column[p],
                'fraction_profile_edges_shared': by_column[p] / choose2(pcounts[p]) if pcounts[p] >= 2 else ''}
               for p, total in sorted(all_profile.items())]
    return result, columns


def read_audit(path):
    receipt = json.loads((path / 'receipt.json').read_text())
    if receipt['status'] != 'complete_audit':
        raise ValueError('Complete alignment audit required')
    return receipt, json.loads((path / 'column_masks.json').read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', type=Path, required=True)
    parser.add_argument('--profile-audit', type=Path, required=True)
    parser.add_argument('--mafft', type=Path, required=True)
    parser.add_argument('--mafft-audit', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    preceipt, pmasks = read_audit(args.profile_audit)
    mreceipt, mmasks = read_audit(args.mafft_audit)
    expected = {r['marker']: r['alignment_sha256'] for r in preceipt['alignment_receipts']}
    mexpected = {r['marker']: r['alignment_sha256'] for r in mreceipt['alignment_receipts']}
    if set(expected) != set(mexpected) or set(pmasks) != set(expected) or set(mmasks) != set(expected):
        raise ValueError('Marker sets differ')
    if args.output.exists():
        raise FileExistsError('Choose a new immutable snapshot')
    rows, column_rows, provenance = [], [], []
    for marker in sorted(expected):
        pp, mp = args.profile / f'{marker}.faa', args.mafft / f'{marker}.faa'
        if sha(pp) != expected[marker] or sha(mp) != mexpected[marker]:
            raise ValueError('Changed audited alignment')
        pr_path = args.profile / f'{marker}.receipt.json'
        mr_path = args.mafft / f'{marker}.receipt.json'
        pr, mr = json.loads(pr_path.read_text()), json.loads(mr_path.read_text())
        if pr['input_sha256'] != mr['input_sha256']:
            raise ValueError('Methods aligned different source proteins')
        sto = args.profile / f'{marker}.sto'
        if sha(sto) != pr['stockholm_sha256']:
            raise ValueError('Changed profile residue coordinates')
        palign = {r.id: str(r.seq).upper() for r in SeqIO.parse(pp, 'fasta')}
        malign = {r.id: str(r.seq).upper() for r in SeqIO.parse(mp, 'fasta')}
        full = {r.id: str(r.seq).upper() for r in AlignIO.read(sto, 'stockholm')}
        if not palign.keys() == malign.keys() == full.keys():
            raise ValueError('Different taxon sets')
        # Recompute occupancy rules instead of trusting serialized masks.
        for sequences, mask in [(palign, pmasks[marker]['0.5']), (malign, mmasks[marker]['0.5'])]:
            actual = [i for i, column in enumerate(zip(*sequences.values()), 1)
                      if sum(aa in AA for aa in column) / len(sequences) >= .5]
            if mask != actual:
                raise ValueError('Mask differs from declared occupancy rule')
        profile, mafft = {}, {}
        for taxon in sorted(full):
            if full[taxon].replace('-', '').replace('.', '') != malign[taxon].replace('-', ''):
                raise ValueError('Source residues differ between alignments')
            retained = pr['retained_stockholm_columns_1based']
            if ''.join(full[taxon][i - 1] for i in retained).replace('.', '-') != palign[taxon]:
                raise ValueError('Profile match columns differ from Stockholm')
            positions = residue_positions(full[taxon])
            for column in pmasks[marker]['0.5']:
                index = retained[column - 1] - 1
                if full[taxon][index] in AA:
                    profile[(taxon, positions[index])] = column
            positions = residue_positions(malign[taxon])
            for column in mmasks[marker]['0.5']:
                if malign[taxon][column - 1] in AA:
                    mafft[(taxon, positions[column - 1])] = column
        result, columns = agreement(profile, mafft)
        rows.append(dict(marker=marker, taxa=len(full), profile_columns=len(pmasks[marker]['0.5']),
                         mafft_columns=len(mmasks[marker]['0.5']), **result))
        column_rows.extend(dict(marker=marker, **r) for r in columns)
        provenance.append({'marker': marker, 'profile_receipt_sha256': sha(pr_path), 'mafft_receipt_sha256': sha(mr_path)})
        print(marker, result['edge_jaccard'], flush=True)
    args.output.mkdir(parents=True)
    for filename, records in [('marker_correspondence.tsv', rows), ('profile_column_correspondence.tsv', column_rows)]:
        with (args.output / filename).open('w') as handle:
            writer = csv.DictWriter(handle, list(records[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(records)
    totals = {key: sum(r[key] for r in rows) for key in ['profile_retained_residues', 'mafft_retained_residues',
              'common_retained_residues', 'profile_edges_on_common_residues', 'mafft_edges_on_common_residues', 'shared_edges']}
    union = totals['profile_edges_on_common_residues'] + totals['mafft_edges_on_common_residues'] - totals['shared_edges']
    receipt = {'markers': len(rows), 'occupancy': .5, 'totals': totals,
        'pooled_edge_jaccard': totals['shared_edges'] / union if union else None,
        'profile_audit_sha256': sha(args.profile_audit / 'receipt.json'),
        'mafft_audit_sha256': sha(args.mafft_audit / 'receipt.json'),
        'script_sha256': sha(Path(__file__)), 'inputs': provenance,
        'interpretation': 'Residue-pair correspondence among residues retained by both methods. Agreement is not alignment accuracy; methods can share errors. Unequal retained residue sets are reported separately.',
        'artifacts': {p.name: sha(p) for p in args.output.iterdir()}}
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k:v for k,v in receipt.items() if k != 'inputs'}, indent=2))


if __name__ == '__main__':
    main()
