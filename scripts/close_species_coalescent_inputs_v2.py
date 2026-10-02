#!/usr/bin/env python3
"""Close all 30 actual coalescent inputs with full readback and two journals."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

from species_coalescent_sources import load
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', required=True, type=Path)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    config = json.loads(Path(plan['source_plan']).read_text())
    _, _, _, _, bindings = load(config, plan['source_plan'])
    bind(bindings, args.plan)
    for path, digest in plan['pins'].items(): bind(bindings, path, digest)
    root = Path(config['output'])
    rp, ap = root / 'receipt.json', Path(plan['readback'])
    producer, reader = [json.loads(path.read_text()) for path in [rp, ap]]
    assert producer['status'] == 'complete_full_coalescent_input_preparation_pending_independent_readback'
    assert reader['status'] == 'passed_full_coalescent_input_independent_split_projection_readback'
    assert producer['plan_sha256'] == reader['plan_sha256'] == sha(plan['source_plan'])
    assert reader['producer_receipt_sha256'] == sha(rp)
    summary = dict(cases=30, alignments=2, cohorts=5, support_policies=3,
                   marker_states=3750, original_gene_split_support_decisions=1781685)
    assert all(producer[key] == reader[key] == value for key, value in summary.items())
    summary.update(retained_internal_split_states=reader['retained_internal_split_states'], case_summaries=reader['case_summaries'])
    for path, record in [(rp, producer), (ap, reader)]:
        bind(bindings, path)
        assert record['scientific_eligibility'] is False
        for source, digest in record['source_hashes'].items(): bind(bindings, source, digest)
    for name, digest in producer['artifacts'].items(): bind(bindings, root / name, digest)
    verify(bindings)
    archive, inner = root / 'completion_archive.json', root / 'completion_closure_plan.json'
    closure = dict(output=str(archive), completed_status='complete_verified_species_coalescent_input_archive',
                   evidence=dict(producer=dict(path=str(rp), expected=dict(status=producer['status'], **{key: producer[key] for key in ['cases', 'marker_states', 'original_gene_split_support_decisions']})),
                                 reader=dict(path=str(ap), expected=dict(status=reader['status'], **summary))),
                   links=[dict(**{'from': 'reader'}, field='producer_receipt_sha256', to=str(rp))],
                   launches=plan['launches'], pins=bindings, summary=summary, scope=plan['scope'])
    with inner.open('x') as handle: handle.write(json.dumps(closure, indent=2) + '\n')
    subprocess.run([sys.executable, 'scripts/record_completed_process_handoffs_v2.py', '--plan', str(inner)], check=True)
    closed = json.loads(archive.read_text())
    assert len(closed['services']) == 2
    result = dict(status='complete_verified_full_species_coalescent_inputs', **summary,
                  full_hash_archive=str(archive), full_hash_archive_sha256=sha(archive),
                  bound_source_hashes=len(closed['source_hashes']), exact_process_journals_checked=2,
                  producer_receipt=str(rp), producer_receipt_sha256=sha(rp),
                  independent_readback=str(ap), independent_readback_sha256=sha(ap),
                  completion_plan_sha256=sha(args.plan), scientific_eligibility=False, scope=plan['scope'])
    with Path(plan['output']).open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: result[key] for key in ['status', 'cases', 'marker_states', 'original_gene_split_support_decisions', 'retained_internal_split_states', 'bound_source_hashes']}), flush=True)


if __name__ == '__main__':
    main()
