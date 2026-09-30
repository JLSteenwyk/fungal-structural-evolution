#!/usr/bin/env python3
"""Publish unit-labeled triad source-work counts only after full verified closure."""
import argparse
import csv
import json
from pathlib import Path
from reference_measurement_union_sources import verify
from reference_triad_design_sources import SUMMARY_FIELDS
from run_ortholog_pair_guide_comparison import sha


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--completion', type=Path, required=True)
    p.add_argument('--table', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); assert not args.table.exists() and not args.output.exists()
    c = json.loads(args.completion.read_text())
    assert c['status'] == 'complete_verified_full_reference_triad_work_design' and len(c['services']) == 2
    verify(c['source_hashes']); summary = c['summary']; assert set(summary) == set(SUMMARY_FIELDS)
    rows = []
    for key, counts in sorted(summary['reference_link_counts'].items()):
        guide, design, parent = key.split('|')
        for field, count in sorted(counts.items()): rows.append([guide, design, parent.split('=')[1], field, 'logical_reference_records', count])
    for key, count in sorted(summary['empty_reference_context_counts'].items()):
        guide, design, parent = key.split('|'); rows.append([guide, design, parent.split('=')[1], 'no_reference_gene', 'source_contexts', count])
    for key, count in sorted(summary['duplicate_design_counts'].items()):
        guide, status, work = key.split('|'); rows.append([guide, '', '', status + '|' + work, 'source_contexts', count])
    for key, count in sorted(summary['triad_model_identity_counts'].items()): rows.append(['', '', '', key, 'unique_ordered_versioned_model_triads', count])
    rows.append(['', '', '', 'correspondence_work_triads', 'unique_ordered_versioned_model_triads', summary['correspondence_work_triads']])
    rows.append(['', '', '', 'potential_correspondence_states_two_masks_eight_orders', 'triad_mask_order_states', summary['potential_correspondence_states']])
    with args.table.open('x') as handle:
        w = csv.writer(handle, delimiter='\t', lineterminator='\n')
        w.writerow(['guide', 'reference_design', 'parent_context_eligible', 'quantity', 'counting_unit', 'count']); w.writerows(rows)
    result = dict(status='published_verified_full_reference_triad_work_counts', completion=str(args.completion), completion_sha256=sha(args.completion),
                  publisher_sha256=sha(__file__), table=str(args.table), table_sha256=sha(args.table), rows=len(rows), checked_source_hashes=len(c['source_hashes']),
                  exact_original_journals=2, scientific_eligibility=False,
                  scope='Exact unit-labeled source-work diagnostics from full independent SQL readback and producer/reader hash/journal closure. Logical references, contexts, ordered physical model triples and future mask/order grid states have distinct units. All source parents/missing structures/lexical choices retained. Source-only work readiness is not measured coverage, residue correspondence, biological orthology, prediction accuracy or calibrated asymmetry; do not pool overlapping guides/designs/ties.')
    with args.output.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__': main()
