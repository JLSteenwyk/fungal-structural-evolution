#!/usr/bin/env python3
"""Close full catalog change replay and release small complete taxon tables."""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    parser.add_argument('--coverage-copy', type=Path, required=True)
    parser.add_argument('--change-copy', type=Path, required=True)
    args = parser.parse_args()
    assert not any(p.exists() for p in [args.receipt, args.coverage_copy, args.change_copy])
    plan = json.loads(args.plan.read_text())
    raw_path = Path(plan['output']) / 'receipt.json'
    raw = json.loads(raw_path.read_text())
    atlas_path = Path(plan['new_catalog_closure'])
    atlas = json.loads(atlas_path.read_text())
    assert atlas['status'] == 'complete_verified_full_completed_retrieval_catalog_refresh'
    assert raw['status'] == 'complete_validated_catalog_comparison'
    assert raw['plan_sha256'] == sha(args.plan)
    assert raw['new_links'] == atlas['proteins_linked']
    reader_path = Path('metadata/completed_afdb_catalog_comparison_readback_20261005_v1.json')
    reader = json.loads(reader_path.read_text())
    assert reader['status'] == 'passed_full_catalog_comparison_source_replay'
    assert reader['source_receipt_sha256'] == sha(raw_path)
    assert reader['links_replayed'] == dict(old=raw['old_links'], new=raw['new_links'])
    assert reader['taxa'] == atlas['taxa'] == 526
    assert reader['dispositions'] == raw['dispositions']
    assert reader['unlinked_in_both'] == raw['unlinked_in_both']
    pins = {}; transports = []
    for prefix, validation in [('completed_afdb_catalog_comparison', Path(plan['comparison_adapter_receipt'])),
                               ('completed_afdb_catalog_comparison_readback', reader_path)]:
        ep = Path(f'metadata/{prefix}_execution_20261005_v1.json')
        tp = Path(f'metadata/{prefix}_transport_20261005_v1.json')
        pp = Path(f'metadata/{prefix}_original_tool_payloads_20261005_v1.json')
        e, t, tool = [json.loads(p.read_text()) for p in [ep, tp, pp]]
        assert e['exit_code'] == 0 and not e['timed_out']
        assert e['receipt_sha256'] == t['validation_sha256'] == sha(validation)
        assert tool['initial']['session_id'] == tool['original_tool_session_id'] == t['original_tool_session_id']
        assert tool['terminal']['exit_code'] == t['original_tool_terminal_exit_code'] == 0
        assert t['whole_wrapper_initial_and_terminal_payloads_matched']
        assert t['manager_start_records'] == t['manager_completion_records'] == 1
        rows = [json.loads(line) for line in subprocess.check_output(
            ['journalctl', '--user', '-u', t['unit'], '--all', '-o', 'json', '--no-pager'], text=True).splitlines()]
        inv = e['invocation_id']
        assert inv == t['invocation_id'] and inv in tool['initial']['output']
        rows = [r for r in rows if inv in [r.get('_SYSTEMD_INVOCATION_ID'), r.get('USER_INVOCATION_ID')]]
        exact = [r for r in rows if r.get('_PID') == str(e['wrapper']['pid'])
                 and r.get('_CMDLINE') == ' '.join(e['wrapper']['cmdline'])]
        assert len(exact) == 2
        assert json.loads(exact[0]['MESSAGE']) == dict(original_wrapper=e['wrapper'], invocation_id=inv)
        assert json.loads(exact[1]['MESSAGE']) == {k: v for k, v in e.items()
            if k not in ['source_hashes', 'artifacts', 'command', 'wrapper', 'child', 'scope']}
        assert not any('Failed with result' in r.get('MESSAGE', '') for r in rows)
        jp = Path(f'metadata/{prefix}_original_whole_journal_20261005_v1.jsonl')
        with jp.open('x') as handle:
            for row in rows: handle.write(json.dumps(row, sort_keys=True) + '\n')
        for mapping in [e['source_hashes'], e['artifacts'], t['source_hashes']]:
            for path, digest in mapping.items(): bind(pins, path, digest)
        for path in [ep, tp, pp, validation, jp]: bind(pins, path)
        transports.append(dict(unit=t['unit'], invocation_id=inv,
            original_tool_session_id=t['original_tool_session_id'], actual_terminal_exit_code=0,
            full_original_payloads_verified=True))
    for mapping in [reader['source_hashes'], atlas['source_hashes'], plan['pins']]:
        for path, digest in mapping.items(): bind(pins, path, digest)
    for path, digest in raw['sources'].items(): bind(pins, path, digest)
    for name, digest in raw['artifacts'].items(): bind(pins, Path(plan['output']) / name, digest)
    for path in [args.plan, raw_path, atlas_path, reader_path, Path(__file__)]: bind(pins, path)
    source = Path(plan['output']) / 'taxon_coverage_change.tsv'
    with source.open() as handle: rows = list(csv.DictReader(handle, delimiter='\t'))
    assert len(rows) == len({r['taxon_id'] for r in rows}) == 526
    summary = defaultdict(Counter); totals = Counter()
    dispositions = ('new_catalog_link', 'lost_catalog_link', 'changed_selected_model', 'unchanged_model')
    for row in rows:
        n = int(row['representative_proteins']); counts = {k: int(row[k]) for k in dispositions}
        assert n > 0 and all(v >= 0 for v in counts.values())
        assert int(row['proteins_with_model_old']) == counts['lost_catalog_link'] + counts['changed_selected_model'] + counts['unchanged_model']
        assert int(row['proteins_with_model_new']) == counts['new_catalog_link'] + counts['changed_selected_model'] + counts['unchanged_model']
        assert int(row['unlinked_in_both']) + sum(counts.values()) == n
        totals.update(counts)
        summary[row['study_role']].update(taxa=1, representative_proteins=n,
            proteins_with_model_old=int(row['proteins_with_model_old']),
            proteins_with_model_new=int(row['proteins_with_model_new']), **counts)
    assert {k: totals[k] for k in dispositions} == {k: raw['dispositions'].get(k, 0) for k in dispositions}
    assert sum(s['representative_proteins'] for s in summary.values()) == atlas['proteins_screened']
    assert summary['ingroup']['taxa'] == 501 and summary['outgroup']['taxa'] == 25
    verify(pins)
    coverage_source = Path(plan['new']) / 'taxon_coverage.tsv'
    for original, copy in [(coverage_source, args.coverage_copy), (source, args.change_copy)]:
        with copy.open('xb') as handle: handle.write(original.read_bytes())
        assert sha(copy) == sha(original)
        bind(pins, copy)
    result = dict(status='complete_verified_full_afdb_catalog_comparison',
        checked_utc=datetime.now(timezone.utc).isoformat(), taxa=526,
        representative_proteins=atlas['proteins_screened'], old_links=raw['old_links'], new_links=raw['new_links'],
        net_link_change=raw['net_link_change'], dispositions={k: totals[k] for k in dispositions},
        unlinked_in_both=raw['unlinked_in_both'], study_role_change={k: dict(v) for k, v in summary.items()},
        public_coverage_table=str(args.coverage_copy), public_change_table=str(args.change_copy),
        original_transports=transports, complete_bound_files=len(pins), source_hashes=pins,
        scientific_eligibility=False, confidence_qualified_atlas_complete=False,
        all_eight_aims_incomplete=True, new_predictions=0, gpu=False,
        scope='Both full old/new link tables independently replayed across the fixed 526taxon universe. '
              'Actual original producer/reader tool terminals, complete wrapper journals and source/output '
              'hashes close. Both public 526row tables are byte-identical source copies. Coverage gains '
              'are newly linked existing models, not new predictions or evolutionary discoveries; confidence, '
              'ESMFold integration, domain/structural families and all eight aims remain unfinished.')
    with args.receipt.open('x') as handle: json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__': main()
