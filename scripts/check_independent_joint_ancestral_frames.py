#!/usr/bin/env python3
"""Read-only native fixtures and full future-grid coordinate/metadata census."""
import argparse
import copy
from datetime import datetime, timezone
import json
from pathlib import Path
import resource
import time

import numpy as np

from ancestral_chain_attempt import sha
from independent_joint_ancestral_frames import GAP, decode, verify_arrays, write_arrays
from independent_native_ancestral_alignment import ALPHABET, fasta_records, project
from independent_native_ancestral_topology import match_trees, parse_tree
from independent_short_sampler_outputs_v2 import NATIVE_ALPHABET, strict_json, mapping_rows
from run_baliphy_reference_preflight import verify


def rejected(action):
    try:
        action()
    except (AssertionError, ValueError, KeyError, TypeError):
        return
    raise AssertionError('Altered joint frame or export accepted')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args(); started = time.monotonic()
    root = args.output.resolve(); root.mkdir(exist_ok=False)
    gate_path = Path('metadata/baliphy_joint_node_logger_software_validation_20261003_v4.json')
    gate = strict_json(gate_path.read_text()); verify({'pins':gate['source_hashes']})
    future_path = Path('metadata/baliphy_joint_node_logger_future_models_20261003_v4.json')
    future = strict_json(future_path.read_text()); verify(future)
    roles_path = Path(future['future_roles']); roles = strict_json(roles_path.read_text())
    assert len(roles) == len({r['chain']['chain_id'] for r in roles}) == 1620
    mapping = Path('results/ancestral/case-local-trees-20260927-v1/ancestral_node_mapping.tsv')
    bindings = {str(p):sha(p) for p in [Path(__file__),Path('scripts/independent_joint_ancestral_frames.py'),
        Path('scripts/independent_native_ancestral_alignment.py'),Path('scripts/independent_native_ancestral_topology.py'),
        Path('scripts/independent_short_sampler_outputs_v2.py'),gate_path,future_path,roles_path,mapping]}
    cached = {}; census = []; anchored_bytes = 0
    for role in roles:
        chain = role['chain']; key = chain['alignment'],chain['tree']
        if key not in cached:
            for name in ['alignment','tree']:
                assert sha(chain[name]) == chain[name+'_sha256']; bindings[chain[name]] = sha(chain[name])
            seqs = fasta_records(Path(chain['alignment']).read_text().splitlines())
            observed = {tip:seq.replace('-','') for tip,seq in seqs.items()}
            tree = parse_tree(Path(chain['tree']).read_text())
            assert sorted(observed) == tree['tips']
            assert all(seq and set(seq) <= set(NATIVE_ALPHABET+'X') for seq in observed.values())
            cached[key] = tree,observed
        tree,observed = cached[key]
        assert len(observed) == chain['proteins']
        candidates = mapping_rows(mapping,chain)
        labels = {node['label']:node for node in tree['nodes']}
        assert len(candidates) == 4 and all(r['source_node'] in labels and labels[r['source_node']]['children'] for r in candidates)
        anchors = sum(map(len,observed.values())); nbytes = 3*4*anchors*2
        anchored_bytes += nbytes
        census.append(dict(chain_id=chain['chain_id'],seed=chain['seed'],chain_role=chain['chain'],
            prior_label=chain['prior_label'],effective_input_group=chain['effective_input_group'],
            original_configuration_ids=chain['original_configuration_ids'],tips=len(observed),
            source_nodes=len(tree['nodes']),candidate_source_nodes=sorted(r['source_node'] for r in candidates),
            anchors_per_frame=anchors,three_frame_candidate_state_category_bytes=nbytes,
            native_execution_checked=False,posterior_qualified=False))
    assert len({r['effective_input_group'] for r in census}) == 135
    assert len({a for r in census for a in r['original_configuration_ids']}) == 324
    results = []; negatives = []; serial_cases = []; fixture = Path('data/software_audits/baliphy-joint-node-logger-20261003-v4')
    observed = {tip:seq.replace('-','') for tip,seq in fasta_records((fixture/'alignment.faa').read_text().splitlines()).items()}
    oracle_lookup = np.full(256,GAP,dtype=np.uint8)
    for index,letter in enumerate(ALPHABET):
        if letter in NATIVE_ALPHABET: oracle_lookup[index] = NATIVE_ALPHABET.index(letter)
    for prior in ['broad','centered','package']:
        directory = next((fixture/'native'/prior/'joint_node_logger').glob('attempt-0001/independent-chain-*'))
        tree_path = directory/'runtime-tree.nwk'; frames_path = directory/'C1.P1.site-property-samples.jsonl'
        matched = match_trees((fixture/'tree.nwk').read_text(),tree_path.read_text())
        candidates = {row['source_node']:row['runtime_node'] for row in matched['rows'] if row['source_node'] not in matched['tips']}
        assert len(candidates) == 4
        frames = [strict_json(line) for line in frames_path.read_text().splitlines()]
        assert [frame['iter'] for frame in frames] == [0,10,20]
        bindings[str(tree_path)] = sha(tree_path); bindings[str(frames_path)] = sha(frames_path)
        for frame in frames:
            summary,arrays = decode(frame,frame['iter'],observed,matched['runtime_labels'],candidates)
            sequences = fasta_records(frame['alignmentLines'])
            old_states,old_free = project(sequences,observed,candidates)
            assert np.array_equal(arrays['states'],oracle_lookup[old_states])
            assert np.array_equal(old_free,np.bincount(arrays['unanchored_candidate_indices'],minlength=4))
            # Independent per-cell oracle uses nongap rank within each native
            # node; it does not call this reader's column/category projection.
            for row,node in enumerate(sorted(candidates)):
                label = candidates[node]; sequence = sequences[label]
                for col_index,col in enumerate(arrays['anchor_columns']):
                    if sequence[col] == '-':
                        assert arrays['states'][row,col_index] == arrays['categories'][row,col_index] == GAP
                    else:
                        offset = len(sequence[:int(col)].replace('-',''))
                        assert arrays['categories'][row,col_index] == frame['catStates'][label]['categories'][offset]
                for index,candidate in enumerate(arrays['unanchored_candidate_indices']):
                    if candidate != row: continue
                    col = int(arrays['unanchored_columns'][index]); offset = len(sequence[:col].replace('-',''))
                    assert arrays['unanchored_states'][index] == frame['catStates'][label]['states'][offset]
                    assert arrays['unanchored_categories'][index] == frame['catStates'][label]['categories'][offset]
            output = root/(prior+'-'+str(frame['iter'])+'.npz')
            write_arrays(output,arrays); verify_arrays(output,arrays)
            bindings[str(output)] = sha(output)
            results.append(dict(prior=prior,**summary))
        frame = frames[0]; ancestor = sorted(candidates.values())[0]; tip = sorted(observed)[0]
        for case in ['missing_ancestor','extra_node','wrong_state','wrong_category','boolean_state','boolean_category',
                     'short_categories','wrong_iteration','wrong_width','missing_alignment','invalid_alignment',
                     'duplicate_alignment_label','unsafe_alignment','lowercase_alignment','wrong_concrete_observation',
                     'wrong_rate_shape','state_dependent_rate','nonfinite_rate','wrong_mean_rate','unexpected_condition']:
            bad = copy.deepcopy(frame); wrong_observed = dict(observed)
            if case == 'missing_ancestor': bad['catStates'].pop(ancestor)
            elif case == 'extra_node': bad['catStates']['invented'] = bad['catStates'][ancestor]
            elif case == 'wrong_state': bad['catStates'][ancestor]['states'][0] = (bad['catStates'][ancestor]['states'][0]+1)%20
            elif case == 'wrong_category': bad['catStates'][ancestor]['categories'][0] = 4
            elif case == 'boolean_state': bad['catStates'][ancestor]['states'][0] = True
            elif case == 'boolean_category': bad['catStates'][ancestor]['categories'][0] = False
            elif case == 'short_categories': bad['catStates'][ancestor]['categories'].pop()
            elif case == 'wrong_iteration': bad['iter'] = True
            elif case == 'wrong_width': bad['alignmentLines'][1] += '-'
            elif case == 'missing_alignment': bad.pop('alignmentLines')
            elif case == 'invalid_alignment': bad['alignmentLines'][1] = 'X'*len(bad['alignmentLines'][1])
            elif case == 'duplicate_alignment_label': bad['alignmentLines'].extend(bad['alignmentLines'][:2])
            elif case == 'unsafe_alignment': bad['alignmentLines'][0] += '\n'
            elif case == 'lowercase_alignment': bad['alignmentLines'][1] = bad['alignmentLines'][1].lower()
            elif case == 'wrong_concrete_observation': wrong_observed[tip] = 'V'+wrong_observed[tip][1:]
            elif case == 'wrong_rate_shape': bad['properties']['rate'].pop()
            elif case == 'state_dependent_rate': bad['properties']['rate'][0][1] += .01
            elif case == 'nonfinite_rate': bad['properties']['rate'][0] = [float('nan')]*20
            elif case == 'wrong_mean_rate': bad['properties']['rate'] = [[2]*20 for _ in range(4)]
            else: bad['conditions'] = {'extra':1}
            rejected(lambda:decode(bad,0,wrong_observed,matched['runtime_labels'],candidates))
            negatives.append(prior+':'+case)
        for case in ['changed_category','changed_state','changed_coordinate','changed_dtype','missing_array']:
            altered = {k:v.copy() for k,v in arrays.items()}
            if case == 'changed_category': altered['categories'][0,0] ^= 1
            elif case == 'changed_state': altered['states'][0,0] ^= 1
            elif case == 'changed_coordinate': altered['anchor_columns'][0] += 1
            elif case == 'changed_dtype': altered['categories'] = altered['categories'].astype(np.int16)
            else: altered.pop('category_rates')
            path = root/(prior+'-bad-'+case+'.npz')
            with path.open('xb') as handle: np.savez(handle,**altered)
            rejected(lambda:verify_arrays(path,arrays)); serial_cases.append(prior+':'+case)
    # Exhaustive joint state/category support with an unknown observed residue,
    # a duplicated anchor column, a deletion and a candidate-only insertion.
    support = 0
    for state,letter in enumerate(NATIVE_ALPHABET):
        for category in range(4):
            frame=dict(iter=0,alignmentLines=['>a','A-','>b',letter+'-','>n','-'+letter],
                catStates={'a':dict(states=[0],categories=[0]),'b':dict(states=[state],categories=[category]),
                           'n':dict(states=[state],categories=[category])},
                properties={'rate':[[1]*20 for _ in range(4)]},conditions={})
            summary,a = decode(frame,0,{'a':'A','b':'X'},{'a','b','n'},{'source':'n'})
            assert a['states'].tolist() == [[GAP,GAP]] and a['categories'].tolist() == [[GAP,GAP]]
            assert a['anchor_columns'].tolist() == [0,0]
            assert a['anchor_tip_indices'].tolist() == [0,1] and a['anchor_tip_offsets'].tolist() == [0,0]
            assert a['unanchored_columns'].tolist() == [1]
            assert a['unanchored_states'].tolist() == [state] and a['unanchored_categories'].tolist() == [category]
            assert summary['unknown_observed_tip_positions'] == 1
            support += 1
    rejected(lambda:strict_json('{"iter":0,"iter":0}'))
    rejected(lambda:strict_json('{"rate":NaN}'))
    census_path = root/'full_future_role_coordinate_census.json'
    results_path = root/'native_fixture_frame_readback.json'
    for p,value in [(census_path,census),(results_path,results)]:
        with p.open('x') as handle: handle.write(json.dumps(value,indent=2)+'\n')
        bindings[str(p)] = sha(p)
    usage = resource.getrusage(resource.RUSAGE_SELF)
    receipt = dict(status='passed_independent_joint_ancestral_frame_software_checks',
        checked_utc=datetime.now(timezone.utc).isoformat(),source_hashes=bindings,
        full_future_role_metadata_checked=1620,effective_inputs=135,original_configuration_aliases=324,
        input_alignment_tree_pairs=len(cached),full_three_frame_candidate_state_category_bytes=anchored_bytes,
        native_fixture_frames_checked=len(results),native_fixture_priors=['broad','centered','package'],
        native_ancestral_records_checked=sum(r['native_ancestors'] for r in results),
        native_ancestral_residue_category_pairs=sum(r['ancestral_pairs'] for r in results),
        exhaustive_state_category_support_cases=support,malformed_frames_rejected=negatives,
        altered_serializations_rejected=serial_cases,strict_json_negative_cases=2,
        elapsed_seconds=time.monotonic()-started,self_cpu_seconds=usage.ru_utime+usage.ru_stime,
        self_peak_rss_bytes=usage.ru_maxrss*1024,new_native_inference_runs=0,
        full_native_joint_sampling_output_checked=False,posterior_qualified=False,scientific_eligibility=False,
        scope='All1620futuremetadata/staticcoordinateinputs and9retainednativefixtureframes, not full native grid output. Independent per-cell category coordinate oracle, existing separate amino-acid projection oracle, all serialized values/dtypes/shapes/keys, all80state/categorysupportcombinations. No production runner, mixing, likelihood proof, posterior acceptance, native restart, GPU or charges. Anchored array storage estimate excludes per-frame coordinates, raw all-node JSON, unanchored residues and filesystem overhead; no ETA or full resource guarantee.')
    with args.receipt.open('x') as handle: handle.write(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['source_hashes','malformed_frames_rejected','altered_serializations_rejected']}))


if __name__ == '__main__': main()
