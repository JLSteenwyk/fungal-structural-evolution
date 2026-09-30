#!/usr/bin/env python3
"""Publish exactly verified full-context work counts with explicit counting units."""
import argparse
import csv
import json
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['plan', 'completion', 'readback', 'output', 'index']: p.add_argument('--' + name, type=Path, required=True)
    args = p.parse_args(); assert not args.output.exists() and not args.index.exists()
    config = json.loads(args.plan.read_text()); rp = Path(config['output']) / 'receipt.json'
    r, a, c = [json.loads(path.read_text()) for path in [rp, args.readback, args.completion]]
    assert c['status'] == 'complete_verified_full_reference_context_measurement_design' and len(c['services']) == 2
    assert r['status'] == 'complete_full_reference_context_measurement_design_pending_independent_readback'
    assert a['status'] == 'passed_full_reference_context_measurement_design_sql_readback'
    assert r['plan_sha256'] == a['plan_sha256'] == sha(args.plan)
    assert a['producer_receipt_sha256'] == sha(rp) == c['source_hashes'][str(rp)]
    assert c['source_hashes'][str(args.readback)] == sha(args.readback)
    bindings = {str(path): sha(path) for path in [args.plan, rp, args.readback, args.completion]}
    summary = c['summary']; assert all(r[k] == a[k] == v for k, v in summary.items())
    rows = []
    for cell, count in sorted(summary['measurement_disposition_counts'].items()):
        guide, design, parent, status = cell.split('|')
        assert guide in ['profile', 'mafft'] and design in ['availability', 'sequence_first'] and parent in ['parent=0', 'parent=1']
        assert status in ['no_reference_gene', 'in_full_reference_measurement_design', 'outside_full_reference_measurement_design', 'identical_model_not_independent', 'reference_model_missing']
        assert isinstance(count, int) and count > 0
        rows.append(dict(guide=guide, design=design, parent_context_eligible=parent[-1], measurement_disposition=status,
                         counting_unit='source_contexts' if status == 'no_reference_gene' else 'duplicate_reference_side_links', count=str(count)))
    assert sum(int(x['count']) for x in rows if x['counting_unit'] == 'duplicate_reference_side_links') == summary['duplicate_reference_links']
    with args.output.open('x') as handle:
        w = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
    with args.output.open() as handle: assert list(csv.DictReader(handle, delimiter='\t')) == rows
    for path, h in bindings.items(): assert sha(path) == h, path
    result = dict(status='published_exact_verified_full_reference_context_measurement_work_counts', table=str(args.output), table_sha256=sha(args.output),
                  rows=len(rows), target_contexts=summary['target_contexts'], duplicate_reference_links=summary['duplicate_reference_links'],
                  script_sha256=sha(__file__), source_hashes=bindings, scientific_eligibility=False,
                  scope='Exact independently checked source/journal-closed work-disposition counts. Empty-reference cells count contexts; other cells count duplicate/reference side links, with both designs/guides/ties and parent exclusions explicit. Units must not be pooled; overlapping logical side links are not unique predictions, independent events, coverage-qualified cohorts or effects. Being in a physical pair design does not override parent ineligibility.')
    with args.index.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__': main()
