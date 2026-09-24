#!/usr/bin/env python3
"""Join every cross-clan interval to audited archive confidence and locations."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(8388608), b''): h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ['candidates', 'extraction', 'audit', 'output']:
        ap.add_argument('--' + name, type=Path, required=True)
    args = ap.parse_args()
    crp = args.candidates / 'receipt.json'; erp = args.extraction / 'receipt.json'
    member_path = args.candidates / 'candidate_members.tsv'
    pins = {str(p): sha(p) for p in [crp, erp, args.audit, member_path, Path(__file__)]}
    cr = json.loads(crp.read_text()); er = json.loads(erp.read_text()); audit = json.loads(args.audit.read_text())
    if cr['status'] != 'complete_cross_clan_candidate_inventory': raise ValueError('Incomplete candidates')
    if audit['status'] != 'passed_all_exported_domain_atoms_and_full_disposition_scope': raise ValueError('Full coordinate audit required')
    if audit['producer_receipt_sha256'] != pins[str(erp)]: raise ValueError('Coordinate audit binding differs')
    if cr['artifacts'][member_path.name] != pins[str(member_path)]: raise ValueError('Candidate artifact differs')
    with member_path.open() as f: members = list(csv.DictReader(f, delimiter='\t'))
    expected = {}
    for r in members:
        if r['interval_id'] in expected and expected[r['interval_id']] != r['model']: raise ValueError('Inconsistent interval/model')
        expected[r['interval_id']] = r['model']
    found = {}; scanned = 0; manifest_hashes = {}
    for proof in er['proofs']:
        label = Path(proof['job']).stem
        p = args.extraction / 'shards' / (label + '.jsonl')
        digest = sha(p)
        if digest != proof['receipt']['artifacts'][p.name]: raise ValueError('Manifest changed')
        manifest_hashes[str(p)] = digest
        with p.open() as f:
            for line in f:
                row = json.loads(line); scanned += 1; iid = row['interval_id']
                if iid not in expected: continue
                if iid in found or row['model_key'] != expected[iid] or row['status'] != 'exported': raise ValueError('Candidate identity/disposition differs')
                mean = float(row['mean_ca_plddt']); fraction = float(row['fraction_ca_plddt_ge70'])
                if not math.isfinite(mean) or not 0 <= mean <= 100 or not 0 <= fraction <= 1: raise ValueError('Invalid confidence')
                if row['residues'] != row['end'] - row['start'] + 1: raise ValueError('Invalid interval length')
                row['archive_path'] = str(p.with_suffix('.tar'))
                row['archive_sha256'] = proof['receipt']['artifacts'][label + '.tar']
                found[iid] = row
    if set(found) != set(expected) or scanned != er['counts']['intervals']: raise ValueError('Incomplete interval scope')
    args.output.mkdir(parents=True, exist_ok=False)
    fields = ['interval_id', 'model_key', 'start', 'end', 'residues', 'mean_ca_plddt',
              'fraction_ca_plddt_ge70', 'archive_path', 'archive_sha256', 'member',
              'pdb_sha256', 'fragment_sequence_sha256', 'source_sha256']
    table = args.output / 'interval_confidence.tsv'
    with table.open('w') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        w.writeheader(); w.writerows(found[k] for k in sorted(found))
    groups = defaultdict(lambda: [set(), set()])
    for r in members:
        if r['model_exclusive_within_pair_cluster'] == '1':
            key = tuple(r[k] for k in ['boundary', 'representative', 'pfam_left', 'pfam_right'])
            groups[key][int(r['side'] == 'right')].add(r['interval_id'])
    workloads = {}
    for threshold in [0, .5, .8, .9]:
        counts = Counter(); unique_pairs = set()
        for key, sides in groups.items():
            selected = [{i for i in side if found[i]['fraction_ca_plddt_ge70'] >= threshold} for side in sides]
            counts[key[0] + ':candidate_entries'] += bool(selected[0] and selected[1])
            for a in selected[0]:
                for b in selected[1]:
                    if found[a]['model_key'] == found[b]['model_key']: raise ValueError('Nonexclusive model pair')
                    counts[key[0] + ':interval_pairs'] += 1
                    unique_pairs.add(tuple(sorted((a, b))))
        workloads[str(threshold)] = dict(counts, unique_unordered_interval_pairs=len(unique_pairs))
    for p, digest in {**pins, **manifest_hashes}.items():
        if sha(p) != digest: raise ValueError('Input changed')
    result = dict(status='complete_cross_clan_confidence_manifest_join', source_hashes=pins,
                  manifest_hashes=manifest_hashes, scanned_intervals=scanned,
                  candidate_intervals=len(found), workloads_by_fraction_ge70=workloads,
                  artifacts={table.name: sha(table)},
                  scope='Full candidate confidence/location join to audited manifests. '
                  'Fraction thresholds are descriptive sensitivity choices, not validation of folds. '
                  'Workloads retain only models exclusive within each accession pair and cluster. '
                  'No coordinates re-extracted or archives rehashed here; member bytes must be checked before alignment. '
                  'No PAE filtering, homology, or evolutionary conclusions.')
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'manifest_hashes'}, indent=2))


if __name__ == '__main__': main()
