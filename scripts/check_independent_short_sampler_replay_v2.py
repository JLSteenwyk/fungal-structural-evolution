#!/usr/bin/env python3
"""Full-grid software serialization plus read-only actual native output checks."""
import argparse
import ast
import copy
from datetime import datetime, timezone
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np

from ancestral_chain_attempt import run_attempt, sha, write_json
from independent_native_ancestral_alignment import fasta_records, native_alignments, project
from independent_short_sampler_outputs_v2 import NATIVE_ALPHABET, SUCCESS, replay, site_frame, strict_json, summarize
from prepare_ancestral_state_traces import project_states as original_projection
import run_independent_short_sampler_replay_v2 as workflow


def rejected(action):
    try: action()
    except (AssertionError, KeyError, ValueError, FileNotFoundError): return
    raise AssertionError('Altered replay accepted')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True); parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args(); root = args.output.resolve(); root.mkdir(exist_ok=False)
    source = Path('metadata/baliphy_reference_sampler_qualification_plan_20261003.json')
    native = json.loads(source.read_text()); jobs = json.loads(Path(native['jobs']).read_text())
    assert len(jobs) == 1620 and len({j['chain']['seed'] for j in jobs}) == 1620
    by_id = {j['chain']['chain_id']: j for j in jobs}
    # Pure native API introspection: no inference, alignment sampling or MCMC.
    program = root / 'alphabet.hs'
    program.write_text('import Bio.Alphabet\nmain = print (getLetters aa)\n')
    binary = Path(jobs[0]['config']['command'][5]); prlimit = Path('/usr/bin/prlimit')
    config = dict(command=[str(prlimit), '--as=' + str(8 * 2**30), '--cpu=120', '--fsize=1048576',
                  '--', str(binary), 'run', str(program)], timeout_seconds=150,
                  pins={str(p): sha(p) for p in [prlimit, binary, program]})
    protocol_receipt = run_attempt(root / 'native-alphabet-probe', config)
    protocol = json.loads(protocol_receipt.read_text()); assert protocol['exit_code'] == 0
    assert ''.join(ast.literal_eval((protocol_receipt.parent / 'stdout.log').read_text().strip())) == NATIVE_ALPHABET
    bindings = {str(source): sha(source), str(protocol_receipt): sha(protocol_receipt),
                str(program): sha(program), str(Path(__file__)): sha(__file__)}
    for name, value in protocol['artifacts'].items(): bindings[str(protocol_receipt.parent / name)] = value
    available = [json.loads(p.read_text()) for p in (Path(native['output']) / 'chains').glob('*.json')]
    successful = [r for r in available if r['status'] == 'full_short_sampler_output_integrity_checked_not_posterior']
    selected = []; arrays = []; negatives = []
    for prior in ['broad', 'centered', 'package']:
        row = min((r for r in successful if r['prior_label'] == prior),
                  key=lambda r: Path(by_id[r['chain_id']]['chain']['alignment']).stat().st_size)
        job = by_id[row['chain_id']]; result, states, free = replay(job, row, native['mapping'])
        base = Path(row['native_receipt']).parent
        with Path(job['chain']['alignment']).open() as handle:
            observed = {k: v.replace('-', '') for k, v in fasta_records(handle).items()}
        with next(base.rglob('C1.P1.fastas')).open() as handle:
            original_states = []; original_free = []
            for _, sequences in native_alignments(handle):
                s, f = original_projection(sequences, observed, result['candidate_mapping'])
                original_states.append(s); original_free.append(f)
        assert np.array_equal(states, np.stack(original_states)) and np.array_equal(free, np.stack(original_free))
        selected.append(result); arrays.append((states, free))
        bindings[row['native_receipt']] = sha(row['native_receipt'])
        checkpoint = Path(native['output']) / 'chains' / (row['chain_id'] + '.json')
        bindings[str(checkpoint)] = sha(checkpoint)
        for name, value in json.loads(Path(row['native_receipt']).read_text())['artifacts'].items():
            bindings[str(base / name)] = value
        print('actual_existing_native_replay', prior, row['chain_id'], states.shape, flush=True)
    # Reuse retained native synthetic ambiguity draws; no inference rerun.
    from independent_short_sampler_outputs import site_frame as frozen_v1_site_frame
    from run_baliphy_reference_preflight import verify
    joint_gate_path=Path('metadata/baliphy_joint_node_logger_software_validation_20261003_v4.json')
    joint_gate=json.loads(joint_gate_path.read_text());verify({'pins':joint_gate['source_hashes']})
    bindings[str(joint_gate_path)]=sha(joint_gate_path)
    fixture_root=Path('data/software_audits/baliphy-joint-node-logger-20261003-v4')
    with (fixture_root/'alignment.faa').open() as h: fixture_observed={k:v.replace('-','') for k,v in fasta_records(h).items()}
    ambiguity_native=[];v1_rejected=[]
    for prior in ['broad','centered','package']:
        base=next((fixture_root/'native'/prior/'original_tip_logger').glob('attempt-*/independent-chain-*'))
        with (base/'C1.P1.fastas').open() as h: fixture_alignments=list(native_alignments(h))
        with (base/'C1.P1.site-property-samples.jsonl').open() as h: fixture_frames=[strict_json(line) for line in h]
        assert len(fixture_alignments)==len(fixture_frames)==3
        for frame,(iteration,sequences) in zip(fixture_frames,fixture_alignments):
            checked=site_frame(frame,iteration,sequences,fixture_observed)
            ambiguity_native.append(dict(prior=prior,iteration=iteration,**checked))
            if checked['tip_log_draw_disagreement_positions']:
                rejected(lambda:frozen_v1_site_frame(frame,iteration,sequences,fixture_observed))
                v1_rejected.append(prior+':'+str(iteration))
        for p in [base/'C1.P1.fastas',base/'C1.P1.site-property-samples.jsonl',fixture_root/'alignment.faa']:bindings[str(p)]=sha(p)
    assert len(ambiguity_native)==9 and sum(x['unknown_observed_tip_positions'] for x in ambiguity_native)==9
    assert v1_rejected and sum(x['tip_log_draw_disagreement_positions'] for x in ambiguity_native)==len(v1_rejected)
    # Complete support at one unknown position, including380 unequal draws.
    support_cases=0
    for alignment_letter in NATIVE_ALPHABET:
        for state in range(20):
            frame=dict(iter=0,catStates={'u':dict(categories=[0],states=[state])},properties={'rate':[[1]*20 for _ in range(4)]},conditions={})
            checked=site_frame(frame,0,{'u':alignment_letter},{'u':'X'})
            assert checked['unknown_observed_tip_positions']==1 and checked['concrete_observed_tip_positions']==0
            assert checked['tip_log_draw_disagreement_positions']==int(alignment_letter!=NATIVE_ALPHABET[state])
            assert checked['tip_logs_joint_trajectory_available'] is False
            support_cases+=1
    # Both logs agreeing on the wrong observed residue is still invalid.
    frame=dict(iter=0,catStates={'u':dict(categories=[0],states=[NATIVE_ALPHABET.index('C')])},properties={'rate':[[1]*20 for _ in range(4)]},conditions={})
    rejected(lambda:site_frame(frame,0,{'u':'C'},{'u':'A'}))
    rejected(lambda:site_frame(frame,0,{'u':'X'},{'u':'X'}))
    rejected(lambda:site_frame(frame,0,{'u':'C'},{'u':'?'}))
    unknown_negative_cases=['both_logs_change_concrete_observation','native_unresolved_letter','invalid_observed_letter']
    input_ambiguities={}
    for alignment_path in {j['chain']['alignment'] for j in jobs}:
        with Path(alignment_path).open() as h: observed_input=fasta_records(h)
        input_ambiguities[alignment_path]=[dict(tip=tip,ungapped_tip_offset=i,observed='X',alignment_logger_letter='A',category_logger_letter='C') for tip,seq in observed_input.items() for i,c in enumerate(seq.replace('-','')) if c=='X']
        bindings[alignment_path]=sha(alignment_path)
    assert sum(bool(input_ambiguities[j['chain']['alignment']]) for j in jobs)==24
    result = selected[0]; states, free = arrays[0]
    base = Path(result['native_receipt']).parent
    with next(base.rglob('C1.P1.fastas')).open() as handle: iteration, sequences = next(native_alignments(handle))
    frame = strict_json(next(next(base.rglob('*site-property-samples.jsonl')).open()))
    chain = by_id[result['chain_id']]['chain']
    with Path(chain['alignment']).open() as handle:
        observed = {k: v.replace('-', '') for k, v in fasta_records(handle).items()}
    tip = sorted(observed)[0]
    for case in ['changed_iteration', 'missing_tip', 'invented_ancestor', 'bad_category', 'bad_state',
                 'short_sequence', 'state_letter_disagreement', 'wrong_rate_dimensions',
                 'negative_rate', 'nonfinite_rate', 'wrong_mean_rate', 'state_dependent_rate', 'unexpected_condition']:
        altered = copy.deepcopy(frame)
        if case == 'changed_iteration': altered['iter'] = 1
        elif case == 'missing_tip': altered['catStates'].pop(tip)
        elif case == 'invented_ancestor': altered['catStates']['invented'] = altered['catStates'][tip]
        elif case == 'bad_category': altered['catStates'][tip]['categories'][0] = 4
        elif case == 'bad_state': altered['catStates'][tip]['states'][0] = 20
        elif case == 'short_sequence': altered['catStates'][tip]['categories'].pop()
        elif case == 'state_letter_disagreement': altered['catStates'][tip]['states'][0] = (altered['catStates'][tip]['states'][0]+1)%20
        elif case == 'wrong_rate_dimensions': altered['properties']['rate'].pop()
        elif case == 'negative_rate': altered['properties']['rate'][0] = [-1] * 20
        elif case == 'nonfinite_rate': altered['properties']['rate'][0] = [float('nan')] * 20
        elif case == 'wrong_mean_rate': altered['properties']['rate'] = [[2] * 20 for _ in range(4)]
        elif case == 'state_dependent_rate': altered['properties']['rate'][0][1] += .01
        else: altered['conditions'] = {'unexpected': 1}
        rejected(lambda: site_frame(altered, iteration, sequences, observed)); negatives.append(case)
    rejected(lambda: strict_json('{"iter":0,"iter":0}')); negatives.append('duplicate_json_key')
    rejected(lambda: strict_json('{"rate":NaN}')); negatives.append('nonfinite_json_constant')
    # A separate variable/gap/X fixture exercises all 22 projected states.
    toy = dict(a='A-CDEFGHIKLMNPQRSTVWY-', b='A-CDEFGHIKLMNPQRSTVWY-',
               n0='ACDEFGHIKLMNPQRSTVWYX-', n1='-XYWVTSRQPNMLKIHGFEDCA')
    observed_toy = {k: toy[k].replace('-', '') for k in ['a', 'b']}; observed_toy['b'] = 'X' + observed_toy['b'][1:]
    candidate_toy = {'c0': 'n0', 'c1': 'n1'}
    projected, unanchored = project(toy, observed_toy, candidate_toy)
    oracle, oracle_free = original_projection(toy, observed_toy, candidate_toy)
    assert np.array_equal(projected, oracle) and np.array_equal(unanchored, oracle_free)
    assert set(''.join(toy[k] for k in ['n0','n1'])) == set('ACDEFGHIKLMNPQRSTVWYX-')
    # Full real role metadata; replay results and source admission are explicitly
    # mocked here, not fabricated production measurements or completion evidence.
    template = selected[0]; states, free = arrays[0]; fake = {}
    for i, job in enumerate(jobs):
        c = job['chain']; row = copy.deepcopy(template)
        row.update(chain_id=c['chain_id'], seed=c['seed'], source_seed=job['source_seed'],
            effective_input_group=c['effective_input_group'], model_input_identity=c['effective_input_group']+'-'+c['prior_label'],
            prior_label=c['prior_label'], chain_role=c['chain'], family=c['family'], original_configuration_ids=c['original_configuration_ids'])
        # Explicitly mocked results exercise every real input's ambiguity metadata.
        unknown=input_ambiguities[c['alignment']]
        for saved in row['frames']:
            saved.update(unknown_observed_tip_positions=len(unknown),
                concrete_observed_tip_positions=saved['tip_category_state_pairs']-len(unknown),
                tip_log_draw_disagreement_positions=int(bool(unknown)),
                tip_log_draw_disagreements=unknown[:1],tip_logs_joint_trajectory_available=False)
        if i in [0, 4]:
            row.update(status='unresolved_original_short_sampler_disposition_retained', frames=[], node_rows=[], candidate_mapping={},
                       projected_state_observations=0, unanchored_candidate_residue_observations=0, negative_branch_pairs_require_review=0)
        fake[c['chain_id']] = row
    stage = root / 'artificial-full-stage'; plan = dict(output=str(stage), mapping=native['mapping'], pins={},
        resources=dict(minimum_free_disk_gib=1), scope='Artificial full serialization fixture only; no production measurements.')
    plan_path = root / 'artificial-plan.json'; write_json(plan_path, plan)
    original_dispositions = {cid: dict(status='full_short_sampler_output_integrity_checked_not_posterior'
        if row['status'] == SUCCESS else 'unsuccessful_sampler_qualification_attempt_retained')
        for cid, row in fake.items()}
    def fake_load(p, path): return jobs, original_dispositions, {}
    def fake_replay(job, disposition, mapping, array_path=None):
        row = copy.deepcopy(fake[job['chain']['chain_id']])
        if row['status'] != SUCCESS: return row
        if array_path is not None:
            with Path(array_path).open('xb') as handle:
                np.savez(handle, states=states, iterations=np.array([0,10,20], dtype=np.int32), unanchored_residue_counts=free)
        return row, states, free
    serialization_rejections = []
    with patch.object(workflow, 'load', side_effect=fake_load), patch.object(workflow, 'replay', side_effect=fake_replay):
        producer = workflow.run(plan_path, enforce_caps=False)
        assert producer['checked_chains'] == 1618 and producer['unresolved_chains'] == 2
        assert producer['complete_groups'] == 403 and producer['intact_chains_in_unresolved_groups'] == 6
        rejected(lambda: workflow.run(plan_path, enforce_caps=False)); serialization_rejections.append('completed_producer_restart')
        for case in ['missing_role', 'changed_seed', 'scientific_acceptance', 'false_ancestral_categories', 'false_joint_tip_logs', 'changed_metric', 'changed_ambiguity_count', 'hidden_draw_disagreement']:
            path = next(p for p in (stage / 'chains').glob('*.json') if json.loads(p.read_text())['status']==SUCCESS and json.loads(p.read_text())['frames'][0]['unknown_observed_tip_positions']); original = path.read_bytes()
            if case == 'missing_role': path.unlink()
            else:
                row = json.loads(original)
                if case == 'changed_seed': row['seed'] += 1
                elif case == 'scientific_acceptance': row['scientific_eligibility'] = True
                elif case == 'false_ancestral_categories': row['ancestral_categories_available'] = True
                elif case == 'false_joint_tip_logs': row['tip_logs_joint_trajectory_available']=True
                elif case == 'changed_ambiguity_count': row['frames'][0]['unknown_observed_tip_positions']+=1
                elif case == 'hidden_draw_disagreement': row['frames'][0]['tip_log_draw_disagreements']=[]
                else: row['projected_state_observations'] += 1
                write_json(path, row)
            rejected(lambda: workflow.run(plan_path, reader=True, enforce_caps=False))
            path.write_bytes(original); serialization_rejections.append(case)
        reader = workflow.run(plan_path, reader=True, enforce_caps=False)
        rejected(lambda: workflow.run(plan_path, reader=True, enforce_caps=False)); serialization_rejections.append('completed_reader_restart')
    rows = list(fake.values())
    for case in ['omitted_failure', 'duplicated_role', 'hidden_ancestral_categories','false_joint_tip_logs']:
        changed = copy.deepcopy(rows)
        if case == 'omitted_failure': changed = [r for r in changed if r['status'] == SUCCESS]
        elif case == 'duplicated_role': changed[-1] = changed[0]
        elif case == 'hidden_ancestral_categories': changed[0]['ancestral_categories_available'] = True
        else: changed[0]['tip_logs_joint_trajectory_available']=True
        rejected(lambda: summarize(changed, jobs)); serialization_rejections.append(case)
    for name in ['independent_short_sampler_outputs_v2', 'run_independent_short_sampler_replay_v2',
                 'independent_short_sampler_outputs',
                 'independent_native_ancestral_alignment', 'independent_native_ancestral_topology',
                 'prepare_ancestral_state_traces', 'ancestral_residue_anchors']:
        p = Path('scripts') / (name + '.py'); bindings[str(p)] = sha(p)
    receipt = dict(status='passed_full_independent_short_sampler_replay_v2_software_contracts',
        checked_utc=datetime.now(timezone.utc).isoformat(), full_real_roles=1620, full_real_groups=405,
        actual_existing_native_roles=3, actual_existing_native_alignments=9, actual_existing_candidate_frames=36,
        native_alphabet_verified=NATIVE_ALPHABET, pure_native_alphabet_receipt=str(protocol_receipt),
        full_producer_and_reader_serialization_checked=True, artificial_failed_roles_retained=2,
        actual_saved_native_ambiguity_frames=ambiguity_native,actual_v1_ambiguity_rejections=v1_rejected,
        full_unknown_state_support_cases=support_cases,unknown_negative_cases_rejected=unknown_negative_cases,
        artificial_unknown_role_count=24,artificial_unknown_positions=producer['unknown_observed_tip_positions'],
        artificial_draw_disagreements=producer['tip_log_draw_disagreement_positions'],tip_logs_joint_trajectory_available=False,posterior_qualified=False,
        intact_artificial_roles_in_unresolved_groups=6, malformed_property_cases_rejected=negatives,
        serialization_cases_rejected=serialization_rejections, source_hashes=bindings, scientific_eligibility=False,
        scope='Full1620real role metadata with explicitly mocked native results/source admission for complete producer/reader serialization. Three existing successful native production roles redecoded read-only across all priors and independently compared against the original anchor projection. Nine retained native synthetic ambiguity frames across allpriors,400unknown20x20supportcases and3extraobservationrejects. Fullmockedgridincludesall24affectedroles andtheirambiguity/disagreementserialization; no joint identity inferred. One pure native alphabet API probe with no inference; 22-state gap/X software fixture. No biological pilot, production restart, inferred ancestral categories, posterior qualification or full production closure.')
    with args.receipt.open('x') as handle: handle.write(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k != 'source_hashes'}), flush=True)


if __name__ == '__main__': main()
