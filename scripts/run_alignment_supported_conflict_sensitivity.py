#!/usr/bin/env python3
"""Complete MAFFT guide-conflict diagnostics and exhaustive readback for both alignments."""
import argparse
import csv
import fcntl
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

from screen_duplication_alignment_reuse import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    bindings = {str(args.plan): sha(args.plan), **plan['pins']}

    def verify():
        for path, digest in bindings.items():
            if sha(path) != digest:
                raise ValueError('Changed conflict source: ' + path)

    verify()
    root = Path(plan['output'])
    root.mkdir(parents=True, exist_ok=True)
    lock = (root / '.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    config = root / 'plan_binding.json'
    binding = dict(plan_sha256=sha(args.plan), script_sha256=sha(__file__))
    if config.exists() and json.loads(config.read_text()) != binding:
        raise ValueError('Changed restart configuration')
    config.write_text(json.dumps(binding, indent=2) + '\n')
    stages = {}
    for stage in plan['stages']:
        verify()
        path = Path(stage['receipt'])
        if not path.exists():
            with (root / (stage['id'] + '.log')).open('a') as handle:
                print('Starting full conflict stage', stage['id'], flush=True)
                subprocess.run([sys.executable, *stage['arguments']], stdout=handle,
                               stderr=subprocess.STDOUT, check=True)
        result = json.loads(path.read_text())
        for field, expected in stage['expected'].items():
            if result[field] != expected:
                raise ValueError('Incomplete conflict stage: ' + stage['id'] + ':' + field)
        for name, digest in result.get('artifacts', {}).items():
            if sha(path.parent / name) != digest:
                raise ValueError('Changed conflict stage artifact: ' + name)
        stages[stage['id']] = dict(receipt=str(path), receipt_sha256=sha(path))
    for label, source in plan['sources'].items():
        proof = json.loads(Path(source['proof']).read_text())
        if proof['diagnostic_receipt_sha256'] != sha(Path(source['diagnostic']) / 'receipt.json'):
            raise ValueError('Conflict proof lineage differs: ' + label)
    key_fields = ['marker', 'guide', 'full_split_taxa', 'sh_alrt_cutoff']

    def load(source):
        path = Path(source['diagnostic']) / 'marker_guide_conflict.tsv'
        records = {}
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                key = tuple(row[k] for k in key_fields)
                if key in records:
                    raise ValueError('Repeated conflict cell')
                records[key] = row
        if len(records) != plan['expected_rows']:
            raise ValueError('Incomplete paired conflict scope')
        return records

    profile, mafft = [load(plan['sources'][label]) for label in ['profile', 'mafft']]
    if profile.keys() != mafft.keys():
        raise ValueError('Full-guide comparison grids differ')
    path = root / 'paired_guide_conflict.tsv'
    fields = key_fields + ['profile_marker_taxa', 'mafft_marker_taxa',
                           'profile_status', 'mafft_status', 'classification_changed']
    totals = Counter()
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for key, first in profile.items():
            second = mafft[key]
            row = dict(zip(key_fields, key))
            row.update(profile_marker_taxa=first['marker_taxa'], mafft_marker_taxa=second['marker_taxa'],
                       profile_status=first['status'], mafft_status=second['status'],
                       classification_changed=int(first['status'] != second['status']))
            writer.writerow(row)
            totals[(key[1], key[3], first['status'], second['status'])] += 1
    # Read the complete serialized join back against both independently audited tables.
    checked = set()
    observed = Counter()
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = tuple(row[k] for k in key_fields)
            if key in checked:
                raise ValueError('Repeated serialized comparison')
            checked.add(key)
            first, second = profile[key], mafft[key]
            assert row['profile_status'] == first['status'] and row['mafft_status'] == second['status']
            assert row['profile_marker_taxa'] == first['marker_taxa'] and row['mafft_marker_taxa'] == second['marker_taxa']
            assert int(row['classification_changed']) == int(first['status'] != second['status'])
            observed[(key[1], key[3], row['profile_status'], row['mafft_status'])] += 1
    assert checked == profile.keys() and observed == totals
    summary = root / 'classification_transition_counts.tsv'
    with summary.open('w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['guide', 'sh_alrt_cutoff', 'profile_status', 'mafft_status', 'marker_edge_cells'])
        writer.writerows([*key, value] for key, value in sorted(totals.items()))
    verify()
    result = dict(status='complete_full_alignment_supported_conflict_sensitivity_with_exhaustive_readbacks',
                  plan_sha256=sha(args.plan), stages=stages, source_hashes=bindings,
                  paired_rows_checked=len(checked), transition_rows=len(totals),
                  classification_changed_rows=sum(v for key, v in totals.items() if key[2] != key[3]),
                  artifacts={p.name: sha(p) for p in [path, summary]}, scientific_eligibility=False,
                  scope='All 125 markers x two full homogeneous guide trees x every internal edge x two SH-aLRT cutoffs. '
                        'Both methods have exhaustive independent set-based conflict maxima, exact supports, witnesses and grid/aggregate checks. '
                        'Serialized status join read back against those source tables, not an independent inference rerun. '
                        'Thirty markers have different original taxon coverage; restricted splits may collapse guide paths. '
                        'No independent-event assumption, preferred alignment, concordance factor, bootstrap probability, species-tree acceptance or biological discordance-cause claim.')
    (root / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True)


if __name__ == '__main__':
    main()
