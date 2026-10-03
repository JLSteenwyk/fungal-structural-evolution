#!/usr/bin/env python3
"""Full selected-chain native-source census; native sample decoding remains pending."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path

from independent_baliphy_scalar_sources import load
from independent_native_ancestral_alignment import ALPHABET, fasta_records, tree_labels
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def run(plan_path):
    plan = json.loads(plan_path.read_text()); source, bindings = load(plan, plan_path)
    producer = source['doc'](source['complete']['producer_receipt'], source['complete']['producer_receipt_sha256'])
    states = producer['states']; assert set(states) == set(source['native'])
    groups = {}; totals = Counter(); largest = Counter(); inputs = {}
    for cid, row in sorted(source['native'].items()):
        native = row['chain']; selected = row['selected_disposition']
        entry = dict(chain_id=cid, model_input_identity=row['model_input_identity'], seed=native['seed'],
            native_disposition=selected['status'], scientific_eligibility=False)
        if selected['status'] != 'all_saved_alignments_and_candidate_nodes_checked':
            entry['status'] = 'unresolved_failed_native_chain_retained'; groups[cid] = entry
            totals['failed_chains'] += 1; continue
        audit = source['doc'](selected['sample_audit'], selected['sample_audit_sha256'])
        receipt = source['doc'](audit['attempt_receipt'], audit['attempt_receipt_sha256'])
        assert receipt['exit_code'] == 0 and audit['iterations'] == 1000
        paths = {suffix:[p for p in receipt['artifacts'] if Path(p).name == suffix]
            for suffix in ['C1.P1.fastas', 'runtime-tree.nwk']}
        assert all(len(p) == 1 for p in paths.values())
        files = {}
        for suffix, keys in paths.items():
            p = Path(audit['attempt_receipt']).parent / keys[0]; digest = receipt['artifacts'][keys[0]]
            bind(bindings, p, digest)
            files[suffix] = dict(path=str(p), sha256=digest, bytes_at_observation=p.stat().st_size)
        ap = native['alignment']; bind(bindings, ap, native['alignment_sha256'])
        cache_key = (str(Path(ap).resolve()), native['alignment_sha256'])
        if cache_key not in inputs:
            with Path(ap).open() as handle: aligned = fasta_records(handle)
            inputs[cache_key] = {tip:seq.replace('-', '') for tip,seq in aligned.items()}
        observed = inputs[cache_key]
        assert len(observed) == native['proteins'] and all(observed.values())
        runtime_labels, runtime_tips = tree_labels(Path(files['runtime-tree.nwk']['path']).read_text())
        assert runtime_tips == set(observed)
        candidates = {}
        lengths = {}
        for sample in audit['candidate_samples']:
            iteration = sample['iteration']; node = sample['source_node']; label = sample['runtime_node']
            assert node not in candidates.setdefault(iteration, {})
            candidates[iteration][node] = label; lengths[(iteration,node)] = sample['ungapped_length']
        assert sorted(candidates) == list(range(0,1001,10))
        first = candidates[0]
        assert len(first) == len(set(first.values())) == 4
        assert all(c == first for c in candidates.values())
        assert set(first.values()) <= runtime_labels - runtime_tips
        state = states[cid]; report = source['doc'](state['receipt'], state['receipt_sha256'])
        assert len(report['summaries']) == 1
        summary = report['summaries'][0]
        assert summary['source_audit_sha256'] == selected['sample_audit_sha256']
        assert summary['samples'] == 101
        arrays = [p for p in summary['artifacts'] if Path(p).name == 'states.npz']
        coords = [p for p in summary['artifacts'] if Path(p).name == 'coordinates.json']
        assert len(arrays) == len(coords) == 1
        coordinates = source['doc'](coords[0], summary['artifacts'][coords[0]])
        assert coordinates['alphabet'] == ALPHABET and coordinates['nodes'] == sorted(first)
        assert coordinates['tips'] == [dict(tip=tip,length=len(observed[tip])) for tip in sorted(observed)]
        assert coordinates['input_alignment_sha256'] == native['alignment_sha256']
        coordinate_count = 4 * sum(len(seq) for seq in observed.values())
        assert summary['candidate_anchor_coordinates'] == coordinate_count
        assert summary['state_observations'] == 101 * coordinate_count
        for p,d in summary['artifacts'].items(): bind(bindings,p,d)
        entry.update(status='full_native_source_metadata_and_runtime_labels_checked_sample_decode_pending',
            sample_audit=selected['sample_audit'], sample_audit_sha256=selected['sample_audit_sha256'],
            selected_attempt_receipt=audit['attempt_receipt'], selected_attempt_receipt_sha256=audit['attempt_receipt_sha256'],
            input_alignment=ap, input_alignment_sha256=native['alignment_sha256'],
            native_files=files, source_to_runtime_candidates=first,
            original_state_array=dict(path=arrays[0],sha256=summary['artifacts'][arrays[0]]),
            coordinates=dict(path=coords[0],sha256=summary['artifacts'][coords[0]]),
            tips=len(observed), runtime_nodes=len(runtime_labels), expected_saved_alignments=101,
            candidate_anchor_coordinates=coordinate_count, state_observations=101*coordinate_count,
            inherited_unanchored_residue_observations=summary['unanchored_residue_observations'])
        groups[cid] = entry; totals['checked_chains'] += 1
        totals['saved_alignments_to_decode'] += 101; totals['candidate_frames_to_decode'] += 404
        totals['state_observations'] += entry['state_observations']
        totals['candidate_anchor_coordinates'] += coordinate_count
        totals['inherited_unanchored_residue_observations'] += summary['unanchored_residue_observations']
        totals['native_alignment_bytes'] += files['C1.P1.fastas']['bytes_at_observation']
        for key, value in [('native_alignment_bytes',files['C1.P1.fastas']['bytes_at_observation']),
                           ('state_array_uncompressed_bytes',entry['state_observations']),
                           ('candidate_anchor_coordinates',coordinate_count), ('runtime_nodes',len(runtime_labels))]:
            largest[key] = max(largest[key],value)
    assert len(groups) == 1620 and totals['checked_chains'] == 1618 and totals['failed_chains'] == 2
    closed = source['complete']
    assert totals['state_observations'] == closed['state_observations']
    assert totals['candidate_anchor_coordinates'] == closed['per_chain_anchor_occurrences']
    assert totals['inherited_unanchored_residue_observations'] == closed['unanchored_residue_observations']
    verify(bindings)
    details = Path(plan['details']); details.parent.mkdir(parents=True,exist_ok=False)
    with details.open('x') as handle: handle.write(json.dumps(groups,indent=2)+'\n')
    result = dict(status='complete_full_independent_native_alignment_inventory_sample_replay_pending',
        checked_utc=datetime.now(timezone.utc).isoformat(), plan_sha256=sha(plan_path),
        full_chains=1620, checked_chains=1618, failed_chains=2, full_groups=405, complete_groups=403,
        original_cutoffs=[250,500], alphabet=ALPHABET, totals=dict(totals), largest=dict(largest),
        full_chain_details=str(details), full_chain_details_sha256=sha(details),
        source_hashes_rechecked=len(bindings), full_source_hash_archive=plan['hash_archive'],
        proposed_replay_resources=dict(cpu_equivalents=2,memory_gib=32,swap_bytes=0,blas_threads=1,
            output_allowance_gib=4,minimum_free_disk_gib=104,uncalibrated_wall_hours_per_stage=[1,48],
            finish_eta=None,production_native_sample_replay_launched=False,gpu=False,new_cost_usd=0),
        scientific_eligibility=False,
        scope='Full1620selected-chain ledger/1618intact/twofailures, including six intact chains in '
              'unresolved quartets. All closed source hashes freshly rechecked, input FASTA and every '
              'runtime tree independently decoded; node/tip/seed/candidate metadata reconciled. '
              'Native saved alignment blocks, per-draw residue projections, unanchored counts and '
              'both cutoff count arrays are not replayed by this inventory. Source-to-runtime clade '
              'mapping inherits prior closed audit; no independent topology, posterior or model qualification.')
    archive = dict(status='full_native_alignment_inventory_source_hash_archive',source_hashes=bindings)
    with Path(plan['hash_archive']).open('x') as handle:handle.write(json.dumps(archive,indent=2)+'\n')
    result['full_source_hash_archive_sha256'] = sha(plan['hash_archive'])
    with Path(plan['output']).open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True)
    args=parser.parse_args();run(args.plan)
