#!/usr/bin/env python3
"""Compare audited PMSF results sharing an alignment but using different guides."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from Bio import SeqIO
from audit_busco_gene_copies import sha


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['run-a', 'run-b', 'audit-a', 'audit-b', 'output']:
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    configs, supports, pins = [], [], {}
    for run, audit in [(a.run_a, a.audit_a), (a.run_b, a.audit_b)]:
        receipt = json.loads((audit / 'receipt.json').read_text())
        if (receipt['status'] != 'passed_pmsf_profile_tree_and_bootstrap_readback'
                or receipt['source_receipt_sha256'] != sha(run / 'receipt.json')):
            raise ValueError('Unbound source audit')
        source = json.loads((run / 'receipt.json').read_text())
        for folder, artifacts in [(run, source['artifacts']), (audit, receipt['artifacts'])]:
            for name, expected in artifacts.items():
                if sha(folder / name) != expected:
                    raise ValueError('Changed source artifact')
        configs.append(json.loads((run / 'config.json').read_text()))
        with (audit / 'branch_support.tsv').open() as handle:
            rows = [r for r in csv.DictReader(handle, delimiter='\t') if r['tree'] == 'ml']
        lookup = {frozenset(json.loads(r['split_taxa_json'])): r for r in rows}
        if len(lookup) != len(rows):
            raise ValueError('Duplicated split')
        supports.append(lookup)
        pins[str(audit / 'receipt.json')] = sha(audit / 'receipt.json')
        pins[str(run / 'receipt.json')] = sha(run / 'receipt.json')
    if configs[0]['matrix_receipt_sha256'] != configs[1]['matrix_receipt_sha256']:
        raise ValueError('Alignment must be identical for this guide comparison')
    if configs[0]['guide_audit_receipt_sha256'] == configs[1]['guide_audit_receipt_sha256']:
        raise ValueError('Distinct guide inputs required')
    cmd = configs[0]['command']
    alignment = Path(cmd[cmd.index('-s') + 1])
    if sha(alignment) != configs[0]['pinned_files'][str(alignment)]:
        raise ValueError('Alignment changed')
    taxa = {r.id for r in SeqIO.parse(alignment, 'fasta')}
    left, right = supports
    common = set(left) & set(right)
    high = [{s for s, row in table.items() if float(row['sh_alrt_percent']) >= 80
             and float(row['empirical_ufboot_percent']) >= 95} for table in supports]
    conflicts = [(x, y) for x in high[0] for y in high[1]
                 if x & y and x - y and y - x and taxa - (x | y)]
    profiles = [np.loadtxt(run / 'pmsf.sitefreq') for run in [a.run_a, a.run_b]]
    if profiles[0].shape != profiles[1].shape or not np.array_equal(profiles[0][:, 0], profiles[1][:, 0]):
        raise ValueError('Profile site grids differ')
    l1 = np.abs(profiles[0][:, 1:] - profiles[1][:, 1:]).sum(axis=1)
    a.output.mkdir(parents=True)
    with (a.output / 'split_comparison.tsv').open('w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['split_taxa_json', 'present_a', 'present_b', 'sh_alrt_a', 'sh_alrt_b',
                         'ufboot_a', 'ufboot_b', 'branch_length_a', 'branch_length_b'])
        for key in sorted(set(left) | set(right), key=lambda s: (len(s), sorted(s))):
            x, y = left.get(key, {}), right.get(key, {})
            writer.writerow([json.dumps(sorted(key)), bool(x), bool(y)] +
                            [r.get(field, '') for field in ['sh_alrt_percent', 'empirical_ufboot_percent', 'branch_length'] for r in [x, y]])
    with (a.output / 'supported_conflicts.tsv').open('w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['split_a_taxa_json', 'split_b_taxa_json'])
        writer.writerows((json.dumps(sorted(x)), json.dumps(sorted(y))) for x, y in
                         sorted(conflicts, key=lambda pair: (sorted(pair[0]), sorted(pair[1]))))
    result = dict(status='complete_same_alignment_two_guide_descriptive_sensitivity',
                  taxa=len(taxa), sites=len(l1), internal_splits_a=len(left), internal_splits_b=len(right),
                  shared_internal_splits=len(common), rf_distance=len(set(left) ^ set(right)),
                  high_support_splits_a=len(high[0]), high_support_splits_b=len(high[1]),
                  incompatible_high_support_split_pairs=len(conflicts),
                  site_profile_l1_mean=float(l1.mean()), site_profile_l1_max=float(l1.max()),
                  inputs=pins, script_sha256=sha(Path(__file__)),
                  artifacts={p.name: sha(p) for p in a.output.iterdir()},
                  scope='Descriptive guide sensitivity on one fixed alignment. High support means SH-aLRT >=80 and empirical UFBoot >=95. Different guide-estimated profiles prevent interpreting cross-run likelihood differences as model preference. Two other crossed runs, rooting, adequacy and gene discordance remain required.')
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
