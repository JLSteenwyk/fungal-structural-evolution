#!/usr/bin/env python3
"""Compare refreshed native candidate exports with the independently audited catalog join."""
import argparse
import csv
import gzip
import json
from pathlib import Path
from estimate_expanded_duplication_pairs import sha


def identity(row, prefix=''):
    assignments = tuple(sorted((row['gene_' + side], row[prefix + 'model_' + side]) for side in ['a', 'b']))
    return row['taxon_id'], assignments, row[prefix + 'same_model']


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['coverage', 'expanded-readback', 'output']:
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    audit = json.loads(a.expanded_readback.read_text())
    assert audit['status'] == 'passed_full_expanded_duplication_event_and_sql_join_readback'
    rp = a.coverage / 'receipt.json'
    receipt = json.loads(rp.read_text())
    assert receipt['status'] == 'complete_duplication_frozen_structure_coverage'
    sources = dict(audit['sources'])
    sources.update({str(a.expanded_readback): sha(a.expanded_readback), str(rp): sha(rp)})
    sources.update({str(a.coverage / name): h for name, h in receipt['artifacts'].items()})
    for path, h in sources.items():
        assert sha(path) == h, path
    candidates = {}
    for guide in ['profile', 'mafft']:
        with (a.coverage / (guide + '_terminal_two_model_candidates.tsv')).open() as f:
            for r in csv.DictReader(f, delimiter='\t'):
                key = (guide, r['family'], r['gene_node'])
                assert key not in candidates
                candidates[key] = identity(r)
    expected_count = len(candidates)
    event_path, = [Path(s) for s in audit['sources'] if s.endswith('/event_coverage.tsv.gz')]
    n = 0
    checked = dict(profile=0, mafft=0)
    with gzip.open(event_path, 'rt') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            n += 1
            if r['new_coverage_class'] != 'both_models':
                continue
            key = (r['guide'], r['family'], r['gene_node'])
            expected = identity(r, 'new_')
            assert candidates.pop(key) == expected, key
            checked[r['guide']] += 1
    assert n == audit['event_rows'] and not candidates and sum(checked.values()) == expected_count
    for guide, count in checked.items():
        assert count == audit['guides'][guide]['new_both_models']
    for path, h in sources.items():
        assert sha(path) == h, path
    result = dict(status='passed_full_expanded_candidate_export_join', events_scanned=n,
                  candidates_checked=checked, source_sha256=sources, script_sha256=sha(__file__),
                  scope='Every candidate identity, taxon, gene/model assignment and same-model flag agrees with the independently audited expanded catalog join. Native support, sequence flags, full noncandidate coverage rows and tree correspondence are outside this check.')
    a.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(checked))


if __name__ == '__main__':
    main()
