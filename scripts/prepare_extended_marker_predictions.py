#!/usr/bin/env python3
"""Freeze 769–1024-residue markers using measured long-control resource observations."""
import argparse
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT, read_table, sha
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs', type=Path, nargs='+', required=True)
    p.add_argument('--predictions', type=Path, nargs='+', required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--minimum-length', type=int, default=769)
    p.add_argument('--maximum-length', type=int, default=1024)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use a new immutable queue')
    if not 769 <= a.minimum_length <= a.maximum_length <= 1024:
        raise ValueError('This resource extrapolation is restricted to lengths 769–1024')
    sequences, origin, links, inventory = {}, {}, [], []
    for folder in a.inputs:
        checked_receipt(folder)
        for record in SeqIO.parse(folder / 'candidates.faa', 'fasta'):
            sid, seq = record.id, str(record.seq)
            if sid in sequences or sid != 'S' + hashlib.sha256(seq.encode()).hexdigest():
                raise ValueError('Overlapping queues or incorrect sequence identity')
            sequences[sid] = seq; origin[sid] = str(folder)
            state = ('noncanonical' if not set(seq) <= set('ACDEFGHIKLMNPQRSTVWY') else
                     'short_queue' if len(seq) < a.minimum_length else
                     'longer_deferred' if len(seq) > a.maximum_length else 'selected_length_band')
            inventory.append({'sequence_id': sid, 'source_queue': str(folder), 'length': len(seq), 'length_disposition': state})
    selected = {r['sequence_id'] for r in inventory if r['length_disposition'] == 'selected_length_band'}
    for folder in a.inputs:
        for row in read_table(folder / 'all_marker_links.tsv'):
            if row['sequence_id'] in selected and origin[row['sequence_id']] == str(folder):
                links.append(row)
    if {r['sequence_id'] for r in links} != selected:
        raise ValueError('Incomplete selected marker links')
    raw = (ROOT / 'data/raw/afdb_models.jsonl').read_bytes().splitlines(keepends=True)
    if raw and not raw[-1].endswith(b'\n'):
        raw.pop()
    frozen = b''.join(raw); latest = {}; reusable = defaultdict(list)
    for line in frozen.decode().splitlines():
        row = json.loads(line); latest[row['uniprot_accession']] = row
    for row in latest.values():
        if row['status'] != 'verified':
            continue
        for model in row['models']:
            sid = 'S' + model['sequence_sha256']
            if sid in selected:
                if sha(ROOT / model['path']) != model['sha256']:
                    raise ValueError('Changed reusable structure')
                reusable[sid].append({'path': model['path'], 'sha256': model['sha256'],
                                      'provider': model.get('provider'), 'tool': model.get('tool')})
    observations = []
    for folder in a.predictions:
        config_hash = sha(folder / 'config.json')
        # Freeze the current list of atomically written per-model receipts.
        for path in sorted(folder.glob('S*.json')):
            if path.name.endswith('.oom.json'):
                continue
            row = json.loads(path.read_text()); sid = row['sequence_id']
            if row['status'] != 'verified_prediction' or row['config_sha256'] != config_hash or row['sequence_sha256'] != sid[1:]:
                raise ValueError('Mixed or invalid prediction receipt')
            if sid in selected:
                for name, digest in row['artifacts'].items():
                    if sha(folder / name) != digest:
                        raise ValueError('Changed local reusable structure')
                reusable[sid].append({'path': str(path), 'sha256': sha(path), 'provider': 'local', 'tool': row['source']})
            if 788 <= row['length'] <= 1063:
                seconds, memory = float(row['inference_seconds']), int(row['peak_gpu_allocated_bytes'])
                if not math.isfinite(seconds) or seconds <= 0 or memory <= 0:
                    raise ValueError('Invalid resource observation')
                observations.append({'receipt_path': str(path), 'receipt_sha256': sha(path),
                                     'length': row['length'], 'seconds': seconds, 'peak_allocated_bytes': memory})
    if len(observations) < 10:
        raise ValueError('Insufficient long-control timing observations')
    dispositions = [{'sequence_id': sid, 'length': len(sequences[sid]),
                     'queue_disposition': 'verified_model_available' if reusable[sid] else 'long_prediction_candidate',
                     'existing_models_json': json.dumps(reusable[sid], sort_keys=True)} for sid in sorted(selected)]
    pending = [r for r in dispositions if r['queue_disposition'] == 'long_prediction_candidate']
    scenarios = []
    for exponent in (2, 3):
        normalized = [r['seconds'] / r['length'] ** exponent for r in observations]
        for percentile in (50, 90):
            coefficient = float(np.percentile(normalized, percentile, method='linear'))
            hours = coefficient * sum(r['length'] ** exponent for r in pending) / 3600
            scenarios.append({'length_exponent': exponent, 'observed_percentile': percentile,
                              'seconds_per_length_power': coefficient, 'projected_inference_gpu_hours': hours})
    memory_bound = max(r['peak_allocated_bytes'] * (a.maximum_length / r['length']) ** 2 for r in observations)
    plan = {'status': 'extrapolated_planning_scenarios_based_on_long_control_timings',
            'observed_receipts': len(observations), 'observed_length_range': [788, 1063],
            'target_length_range': [a.minimum_length, a.maximum_length], 'pending_sequences': len(pending),
            'scenarios': scenarios, 'planning_hours_with_50percent_overhead': 1.5 * max(r['projected_inference_gpu_hours'] for r in scenarios),
            'quadratic_total_allocated_memory_extrapolation_bytes': memory_bound,
            'planned_gpu_count': 1, 'planned_cpu_threads': 4, 'planned_output_gb': 100,
            'execution': 'Not launched or queued. Use a separate output configuration and wait for authorized GPU availability. Preserve the existing stop-on-OOM behavior. Recheck exact-sequence reuse before launch.',
            'limitations': 'Quadratic/cubic timing and quadratic total-memory scaling are explicit scenarios, not validated complexity bounds or runtime/confidence intervals. Runtime, allocator overhead and fragmentation can differ between controls and marker proteins. Timing receipts are observations, not an independent structure-quality audit. No extrapolation above 1024 residues or to the full proteome atlas.'}
    a.output.mkdir(parents=True)
    (a.output / 'afdb_snapshot.jsonl').write_bytes(frozen)
    with (a.output / 'candidates.faa').open('w') as handle:
        for row in pending:
            handle.write('>' + row['sequence_id'] + '\n' + sequences[row['sequence_id']] + '\n')
    for name, rows in [('source_length_inventory.tsv', inventory), ('queue_disposition.tsv', dispositions),
                       ('all_marker_links.tsv', links), ('resource_observations.tsv', observations)]:
        write_table(a.output / name, rows)
    (a.output / 'resource_plan.json').write_text(json.dumps(plan, indent=2) + '\n')
    result = {'status': 'complete_long_marker_input_preparation', 'selected_length_band_sequences': len(selected),
              'prediction_candidates': len(pending), 'reusable_sequences': len(selected) - len(pending),
              'selected_taxon_marker_links': len(links), 'selected_taxa': len({r['taxon_id'] for r in links}),
              'selected_markers': len({r['marker'] for r in links}),
              'source_dispositions': dict(Counter(r['length_disposition'] for r in inventory)),
              'source_receipts': {str(folder): sha(folder / 'receipt.json') for folder in a.inputs},
              'prediction_configs': {str(folder): sha(folder / 'config.json') for folder in a.predictions},
              'script_sha256': sha(Path(__file__)), 'resource_plan': plan,
              'interpretation': 'Full selected length band from disjoint original and follow-on marker queues. Complete proteins and all source marker links retained. Exact sequence reuse checked at this frozen snapshot. No prediction execution or structural coverage claim.',
              'artifacts': {path.name: sha(path) for path in a.output.iterdir()}}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
