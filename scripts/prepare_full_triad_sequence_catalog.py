#!/usr/bin/env python3
"""Inventory every original triad and full sequence for correspondence controls."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path

from full_triad_sequence_sources import load_original
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def run(plan_path):
    plan = json.loads(Path(plan_path).read_text())
    assert shutil.disk_usage('.').free >= plan['resources']['minimum_free_disk_gib'] * 2**30
    catalog, ready, inputs, mapping, bindings = load_original(plan, plan_path)
    sets, links = {}, []
    for triad in catalog:
        link = dict(triad, sequence_control_disposition='no_source_ready_logical_link',
                    sequence_set_id=None, role_to_sequence_indices=None)
        if triad['source_design_ready_links']:
            models = sorted(triad['models'])
            identity = json.dumps(models, separators=(',', ':'))
            sid = hashlib.sha256(identity.encode()).hexdigest()
            sequences = [inputs[(*m, 'full')]['sequence'] for m in models]
            record = dict(sequence_set_id=sid, models=models, sequences=sequences,
                          sequence_sha256=[hashlib.sha256(s.encode()).hexdigest() for s in sequences],
                          original_lengths=list(map(len, sequences)))
            assert sid not in sets or sets[sid] == record
            sets[sid] = record
            link.update(sequence_control_disposition='scheduled_full_sequence_alignment',
                        sequence_set_id=sid,
                        role_to_sequence_indices=[models.index(m) for m in triad['models']])
        links.append(link)
    out = Path(plan['output']); out.mkdir(parents=True, exist_ok=False)
    for name, records in [('sequence_sets.jsonl', [sets[k] for k in sorted(sets)]),
                          ('triad_sequence_links.jsonl', links)]:
        with (out / name).open('x') as f:
            for r in records:
                f.write(json.dumps(r, separators=(',', ':')) + '\n')
    verify(bindings)
    summary = dict(ordered_model_triads=len(catalog), source_ready_triads=len(ready),
                   unscheduled_original_triads=len(catalog)-len(ready),
                   unique_sequence_model_sets=len(sets), native_alignment_states=12*len(sets),
                   target_contexts=mapping['target_contexts'],
                   reference_tie_records=mapping['reference_tie_records'],
                   duplicate_reference_links=mapping['duplicate_reference_links'],
                   maximum_sequence_length=max(len(s) for r in sets.values() for s in r['sequences']),
                   maximum_three_sequence_length=sum(sorted([len(s) for r in sets.values() for s in r['sequences']], reverse=True)[:3]),
                   sequence_letters=''.join(sorted({c for r in sets.values() for s in r['sequences'] for c in s})))
    # Exact maximum across real sets, not three unrelated long proteins.
    summary['maximum_three_sequence_length'] = max(sum(r['original_lengths']) for r in sets.values())
    result = dict(status='complete_full_triad_sequence_catalog_pending_independent_readback',
                  plan_sha256=sha(plan_path), **summary, source_hashes=bindings,
                  artifacts={name:sha(out/name) for name in ['sequence_sets.jsonl','triad_sequence_links.jsonl']},
                  scientific_eligibility=False, scope=plan['scope'])
    with (out/'receipt.json').open('x') as f: f.write(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']}), flush=True)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True)
    run(p.parse_args().plan)
