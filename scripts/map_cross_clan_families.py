#!/usr/bin/env python3
"""Link all cross-clan candidate members to source proteins and both family guides."""
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path
import sqlite3


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8388608), b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--candidates', type=Path, required=True)
    ap.add_argument('--composition', type=Path, required=True)
    ap.add_argument('--composition-audit', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    table = a.candidates / 'candidate_members.tsv'
    database = a.composition / 'domain_cluster_composition.sqlite'
    cr, dr = a.candidates / 'receipt.json', a.composition / 'receipt.json'
    sources = {str(p): sha(p) for p in [table, database, cr, dr, a.composition_audit, Path(__file__)]}
    c, d, audit = [json.loads(p.read_text()) for p in [cr, dr, a.composition_audit]]
    assert c['status'] == 'complete_cross_clan_candidate_inventory'
    assert audit['status'] == 'passed_full_domain_cluster_composition_source_readback'
    assert audit['producer_receipt_sha256'] == sources[str(dr)]
    assert c['artifacts'][table.name] == sources[str(table)]
    assert d['artifacts'][database.name] == sources[str(database)]
    db = sqlite3.connect('file:' + str(database) + '?mode=ro', uri=True)
    with table.open() as f:
        reader = csv.DictReader(f, delimiter='\t')
        header = reader.fieldnames
        members = list(reader)
    mapped = []
    for r in members:
        identity = db.execute('SELECT i.model,m.representative FROM intervals i JOIN members m USING(interval_id) WHERE interval_id=?', (r['interval_id'],)).fetchall()
        assert identity == [(r['model'], r['representative'])]
        proteins = db.execute('SELECT native_gene_id,taxon_id,protein_id FROM proteins WHERE model=?', (r['model'],)).fetchall()
        assert proteins
        for gene, taxon, protein in sorted(proteins):
            families = dict(db.execute('SELECT guide,family FROM family_links WHERE native_gene_id=? AND guide IN (?,?)', (gene, 'profile', 'mafft')))
            assert set(families) == {'profile', 'mafft'}
            mapped.append({**r, 'native_gene_id': gene, 'taxon_id': taxon, 'protein_id': protein,
                           'profile_family': families['profile'], 'mafft_family': families['mafft']})
    db.close()
    keys = ['boundary', 'representative', 'pfam_left', 'pfam_right']
    groups = defaultdict(list)
    for r in mapped:
        groups[tuple(r[k] for k in keys)].append(r)
    summaries = []
    attributes = ['model', 'native_gene_id', 'taxon_id', 'profile_family', 'mafft_family']
    for key, rows in sorted(groups.items()):
        for scope in ['all_members', 'exclusive_models']:
            selected = [r for r in rows if scope == 'all_members' or r['model_exclusive_within_pair_cluster'] == '1']
            s = dict(zip(keys, key)); s['membership_scope'] = scope
            for attr in attributes:
                left, right = [{r[attr] for r in selected if r['side'] == side} for side in ['left', 'right']]
                s[attr + '_left'] = len(left)
                s[attr + '_right'] = len(right)
                s[attr + '_shared'] = len(left & right)
            summaries.append(s)
    assert len(groups) == sum(c['candidates_by_boundary'].values())
    a.output.mkdir(parents=True, exist_ok=False)
    for filename, rows, fields in [('candidate_proteins.tsv', mapped, header + ['native_gene_id', 'taxon_id', 'protein_id', 'profile_family', 'mafft_family']),
                                    ('candidate_family_summary.tsv', summaries, list(summaries[0]))]:
        with (a.output / filename).open('w') as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
            w.writeheader(); w.writerows(rows)
    for p, digest in sources.items():
        assert sha(p) == digest, p
    receipt = dict(status='complete_cross_clan_source_protein_family_mapping', source_hashes=sources,
                   candidate_entries=len(groups), input_members=len(members), mapped_rows=len(mapped),
                   proteins=len({r['native_gene_id'] for r in mapped}), taxa=len({r['taxon_id'] for r in mapped}),
                   artifacts={p.name: sha(p) for p in a.output.glob('*.tsv')},
                   scope='All candidate members retained, including shared models. Exclusive-model summaries use pair/cluster-local exclusivity. Counts precede confidence and direct-alignment screening. Shared family assignments are not proof of orthology or homology; shared taxa are not evolutionary transitions. Boundary views overlap. One structure can map to multiple source proteins.')
    (a.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
