#!/usr/bin/env python3
"""Prepare all corrected future models and four-chain roles after native software proof."""
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha, write_json
from baliphy_joint_node_logger_v5 import transform, restore, fresh_seeds
from run_baliphy_reference_preflight import verify


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--validation', type=Path, required=True)
    p.add_argument('--execution', type=Path, required=True)
    p.add_argument('--terminal', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    a = p.parse_args(); assert not a.receipt.exists()
    gate = json.loads(a.validation.read_text()); verify({'pins': gate['source_hashes']})
    assert gate['status'] == 'passed_joint_node_logger_v5_cjson_full_source_and_native_qualification'
    assert (gate['full_programs_checked'], gate['full_roles_checked'], gate['native_runs']) == (405, 1620, 9)
    assert (gate['same_record_frames_checked'], gate['ancestral_records_checked']) == (18, 72)
    assert gate['independent_all_node_sequence_category_coordinate_checks']
    assert gate['serialized_candidate_array_checks'] == 18
    assert gate['strict_mean_one_assertion_unchanged'] and gate['native_exact_double_roundtrips'] == 720
    assert len(gate['altered_frames_rejected_by_both_decoders']) == 108
    assert gate['original_encoder_normalization_failures_reproduced']
    execution = json.loads(a.execution.read_text())
    assert execution['status'] == 'exited_zero_with_receipt' and execution['exit_code'] == 0
    assert execution['receipt_sha256'] == sha(a.validation)
    verify({'pins': execution['source_hashes']})
    for path, digest in execution['artifacts'].items(): assert sha(path) == digest
    terminal = json.loads(a.terminal.read_text())
    assert terminal['actual_tool_session_id'] == 64924 and terminal['actual_terminal_exit_code'] == 0
    assert terminal['invocation_id'] == execution['invocation_id']
    verify({'pins': terminal['source_hashes']})
    source = Path('metadata/baliphy_reference_sampler_qualification_plan_20261003.json')
    native = json.loads(source.read_text()); verify(native)
    jobs = json.loads(Path(native['jobs']).read_text()); assert len(jobs) == 1620
    previous = Path('metadata/baliphy_joint_node_logger_future_models_20261003_v4.json')
    old = json.loads(previous.read_text()); verify(old)
    old_roles = Path(old['future_roles'])
    used = {j['chain']['seed'] for j in jobs} | {j['source_seed'] for j in jobs}
    used.update(r['chain']['seed'] for r in json.loads(old_roles.read_text()))
    used.update({20267001, 20267002, 20267003, 20268001, 20268002, 20268003})
    assert len(used) == gate['forbidden_seed_count']
    namespace = gate['future_seed_namespace']
    seeds = fresh_seeds([j['chain']['chain_id'] for j in jobs], used, namespace)
    root = a.output.resolve(); root.mkdir(exist_ok=False); (root/'models').mkdir()
    models = {}; roles = []; groups = defaultdict(list); pins = dict(native['pins'])
    for job in jobs:
        old_chain = job['chain']; group = old_chain['effective_input_group']+'-'+old_chain['prior_label']
        source_path = Path(old_chain['program']); assert sha(source_path) == old_chain['program_sha256']
        text = source_path.read_text(); changed = transform(text); target = root/'models'/(group+'.hs')
        if group not in models:
            target.write_text(changed); assert restore(target.read_text()) == text
            assert 'logEncoderChecks' not in changed and 'encoder-checks.jsonl' not in changed
            models[group] = dict(model_input_identity=group, prior_label=old_chain['prior_label'],
                family=old_chain['family'], source_program=str(source_path), source_program_sha256=sha(source_path),
                program=str(target), program_sha256=sha(target), reversible_logger_edit_verified=True,
                probability_model_changed=False, native_cjson_property_encoding=True,
                software_fixture_alpha_or_oracle_included=False)
            pins[str(target)] = sha(target)
        else: assert target.read_text() == changed
        chain = dict(old_chain, chain_id=group+'-joint-cjson-v5-chain'+str(old_chain['chain']),
            seed=seeds[old_chain['chain_id']], program=str(target), program_sha256=sha(target))
        roles.append(dict(chain=chain, source_chain_id=old_chain['chain_id'], earlier_source_seed=job['source_seed'],
            current_short_sampler_seed=old_chain['seed'], seed_namespace=namespace,
            native_execution_launched=False, scientific_eligibility=False, posterior_qualified=False))
        groups[group].append(chain)
    assert len(models) == len(groups) == 405 and len(roles) == 1620
    assert all(len(v) == 4 and {c['chain'] for c in v} == {1, 2, 3, 4} for v in groups.values())
    assert Counter(r['chain']['prior_label'] for r in roles) == {'broad': 540, 'centered': 540, 'package': 540}
    assert len({r['chain']['effective_input_group'] for r in roles}) == 135
    assert len({alias for r in roles for alias in r['chain']['original_configuration_ids']}) == 324
    assert len({r['chain']['seed'] for r in roles}) == len({r['chain']['chain_id'] for r in roles}) == 1620
    assert {r['chain']['seed'] for r in roles}.isdisjoint(used)
    role_path = root/'future_roles.json'; model_path = root/'model_manifest.json'
    write_json(role_path, roles); write_json(model_path, [models[k] for k in sorted(models)])
    for path in [source, a.validation, a.execution, a.terminal, Path(__file__),
                 Path('scripts/baliphy_joint_node_logger_v5.py'),
                 Path('scripts/baliphy_joint_node_logger_v3.py'),
                 Path('scripts/baliphy_full_node_logger.py'), previous, old_roles, role_path, model_path]:
        pins[str(path)] = sha(path)
    verify({'pins': pins})
    result = dict(status='prepared_full_future_joint_node_logger_v5_cjson_sources_not_launched',
        checked_utc=datetime.now(timezone.utc).isoformat(), source_plan=str(source), source_plan_sha256=sha(source),
        models=405, roles=1620, model_quartets=405, effective_inputs=135, original_configuration_aliases=324,
        model_manifest=str(model_path), future_roles=str(role_path), seed_namespace=namespace,
        forbidden_seed_count=len(used), full_source_grid_reread=True, native_cjson_property_encoding=True,
        production_horizon=None, execution_commands_prepared=False, native_execution_launched=False,
        scientific_eligibility=False, posterior_qualified=False, gpu=False, new_cost_usd=0, pins=pins,
        scope='All405originalmodels with five exactly reversible pure logger edits and1620newseedroles. Native CJSON properties, same-record alignment and full labels retained; model/priors/initializer unchanged. Static full coverage plus prior native software gate does not prove full-grid startup/saved-frame execution or posterior adequacy. Previous failed runs preserved; no old source edits/restarts, GPU or charges.')
    with a.receipt.open('x') as f: f.write(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'pins'}, indent=2))


if __name__ == '__main__': main()
