#!/usr/bin/env python3
"""Publish unit-labeled full fit dispositions and separate/combined screen counts."""
import argparse
import csv
import json
from pathlib import Path
from full_triad_fit_sources import SUMMARY_FIELDS
from run_ortholog_pair_guide_comparison import sha


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--completion', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); c = json.loads(args.completion.read_text())
    assert c['status'] == 'complete_verified_full_triad_same_residue_geometry' and c['exact_process_journals_checked'] == 2
    assert sha(c['producer_receipt']) == c['producer_receipt_sha256'] and sha(c['independent_readback']) == c['independent_readback_sha256']
    assert sha(c['full_hash_archive']) == c['full_hash_archive_sha256']
    r = json.loads(Path(c['producer_receipt']).read_text()); a = json.loads(Path(c['independent_readback']).read_text())
    assert a['status'] == 'passed_full_triad_same_residue_quaternion_readback' and a['producer_receipt_sha256'] == c['producer_receipt_sha256']
    assert all(c[key] == r[key] == a[key] for key in SUMMARY_FIELDS)
    for name, digest in r['artifacts'].items(): assert sha(Path(c['producer_receipt']).parent / name) == digest
    rows = []
    for key, value in sorted(c['counts'].items()):
        mask, definition, status = key.split(':'); rows.append(dict(section='fit_status', mask=mask, mapping_definition=definition, screen='', disposition=status, unit='definition_fit_order_state', count=value))
    assert sum(c['counts'].values()) == c['fit_rows'] == 865792
    for key in sorted(c['screen_pass_counts']):
        mask, definition, screen = key.split(':'); core = c['core_screen_pass_counts'][key]; combined = c['screen_pass_counts'][key]
        assert 0 <= combined <= core <= c['mapping_states']
        for kind, value in [('shared_core_pass', core), ('shared_core_and_inherited_three_edge_pass', combined)]:
            rows.append(dict(section='screen_pass', mask=mask, mapping_definition=definition, screen=screen, disposition=kind, unit='definition_fit_order_state', count=value))
    table = Path('docs/tables/full_triad_geometry_dispositions_20261001.tsv')
    with table.open('x') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n'); writer.writeheader(); writer.writerows(rows)
    proof = dict(status='published_bound_full_triad_same_residue_fit_counts', completion_sha256=sha(args.completion), producer_receipt_sha256=c['producer_receipt_sha256'],
                 independent_readback_sha256=c['independent_readback_sha256'], table=str(table), table_sha256=sha(table), rows=len(rows),
                 script_sha256=sha(__file__), scientific_eligibility=False,
                 scope='All completed865792definition fit states and separate shared-core/combined3edge six-screen counts checked against complete producer/independent reader/journal-closed summary. Original large fitted table hash rechecked. Alternative nativeorders/masks/definitions are repeated states, not independent events or accepted evolutionary contrasts.')
    with args.output.open('x') as f: f.write(json.dumps(proof, indent=2) + '\n')
    print(json.dumps(proof, indent=2))


if __name__ == '__main__': main()
