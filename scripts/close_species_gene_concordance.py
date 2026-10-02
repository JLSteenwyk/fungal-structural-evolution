#!/usr/bin/env python3
"""Close full native gene concordance only after full readback and two journals."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha
from species_gcf_sources import load


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    config = json.loads(Path(plan['source_plan']).read_text())
    _, _, _, bindings = load(config, plan['source_plan'])
    bind(bindings, args.plan)
    for path, digest in plan['pins'].items(): bind(bindings, path, digest)
    root = Path(config['output'])
    rp, ap = root / 'receipt.json', Path(plan['readback'])
    producer, reader = [json.loads(path.read_text()) for path in [rp, ap]]
    assert producer['status'] == 'complete_full_native_species_gene_concordance_batch_pending_independent_readback'
    assert reader['status'] == 'passed_full_native_species_gene_concordance_independent_readback'
    assert producer['plan_sha256'] == reader['plan_sha256'] == sha(plan['source_plan'])
    assert reader['producer_receipt_sha256'] == sha(rp)
    summary = dict(references=8, marker_alignments=2, markers_per_alignment=125,
                   native_runs=16, branch_marker_cells=1046000)
    assert all(producer[k] == reader[k] == v for k, v in summary.items())
    assert reader['branch_summary_rows'] == 8368
    assert sum(reader['state_counts'].values()) == 1046000 and len(reader['run_summaries']) == 16
    summary.update(branch_summary_rows=8368, state_counts=reader['state_counts'],
                   run_summaries=reader['run_summaries'])
    for path, record in [(rp, producer), (ap, reader)]:
        bind(bindings, path)
        assert record['scientific_eligibility'] is False
        for source, digest in record['source_hashes'].items(): bind(bindings, source, digest)
    for name, digest in reader['artifacts'].items(): bind(bindings, root / name, digest)
    verify(bindings)
    archive = root / 'completion_archive.json'
    inner = root / 'completion_closure_plan.json'
    closure = dict(output=str(archive), completed_status='complete_verified_species_gene_concordance_archive',
                   evidence=dict(producer=dict(path=str(rp), expected=dict(status=producer['status'], **{k: producer[k] for k in ['references', 'native_runs', 'branch_marker_cells']})),
                                 reader=dict(path=str(ap), expected=dict(status=reader['status'], **summary))),
                   links=[dict(**{'from': 'reader'}, field='producer_receipt_sha256', to=str(rp))],
                   launches=plan['launches'], pins=bindings, summary=summary, scope=plan['scope'])
    with inner.open('x') as handle: handle.write(json.dumps(closure, indent=2) + '\n')
    subprocess.run([sys.executable, 'scripts/record_completed_order_marker_handoffs.py', '--plan', str(inner)], check=True)
    closed = json.loads(archive.read_text())
    assert len(closed['services']) == 2
    result = dict(status='complete_verified_full_native_species_gene_concordance', **summary,
                  full_hash_archive=str(archive), full_hash_archive_sha256=sha(archive),
                  bound_source_hashes=len(closed['source_hashes']), exact_process_journals_checked=2,
                  producer_receipt=str(rp), producer_receipt_sha256=sha(rp),
                  independent_readback=str(ap), independent_readback_sha256=sha(ap),
                  completion_plan_sha256=sha(args.plan), scientific_eligibility=False, scope=plan['scope'])
    with Path(plan['output']).open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['status', 'native_runs', 'branch_marker_cells', 'branch_summary_rows', 'bound_source_hashes', 'state_counts']}), flush=True)


if __name__ == '__main__':
    main()
