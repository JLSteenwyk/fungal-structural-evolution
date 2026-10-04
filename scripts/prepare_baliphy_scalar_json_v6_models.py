#!/usr/bin/env python3
"""Prepare the complete future V6 model/seed grid after actual native qualification.

No native commands, horizon or launch are prepared. Existing V5/failed sources
are immutable; production requires V6-aware full-grid audits and sampler gates.
"""
import argparse
from collections import Counter,defaultdict
from datetime import datetime,timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha,write_json
from baliphy_scalar_json_logger_v6c import HELPER,SCHEMA,transform,restore
from baliphy_joint_node_logger_v5 import fresh_seeds
from check_baliphy_scalar_json_logger_v6_v9 import seed_census,SEED_NAMESPACE
from reference_measurement_union_sources import bind,verify


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--validation',type=Path,required=True);p.add_argument('--execution',type=Path,required=True)
    p.add_argument('--transport',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--receipt',type=Path,required=True);a=p.parse_args();assert not a.receipt.exists()
    gate=json.loads(a.validation.read_text());e=json.loads(a.execution.read_text());t=json.loads(a.transport.read_text())
    assert gate['status']=='passed_scalar_json_v6_full_source_and_paired_native_qualification'
    assert (gate['full_programs_checked'],gate['full_roles_checked'],gate['native_runs'],gate['scalar_rows_checked'])==(405,1620,7,63)
    assert len(gate['paired_prior_checks'])==3 and (gate['joint_frames_checked'],gate['serialized_array_checks'])==(9,9)
    assert (gate['native_finite_value_roundtrips'],gate['native_nonfinite_tags_checked'],gate['native_literal_null_tags_checked'])==(17,6,2)
    assert len(gate['malformed_record_cases_rejected_by_both_readers'])==24
    assert gate['native_smallest_subnormal_construction_and_python_readback']==1
    assert gate['exactly_reversible_from_v5'] and gate['context_action_evaluated_once']
    assert e['status']=='exited_zero_with_receipt' and e['exit_code']==0 and not e['timed_out']
    assert e['receipt_sha256']==sha(a.validation)
    assert t['status']=='verified_original_scalar_v6_software_wait_zero' and t['actual_tool_terminal_exit_code']==0
    assert t['invocation_id']==e['invocation_id'] and t['wrapper']==e['wrapper'] and t['entire_terminal_payload_matched']
    assert t['source_hashes'][str(a.validation)]==sha(a.validation)
    verify(gate['source_hashes']);verify(t['source_hashes'])
    pins=dict(t['source_hashes'])
    old_path=Path('metadata/baliphy_joint_node_logger_future_models_20261003_v5.json')
    old=json.loads(old_path.read_text());verify(old['pins']);old_roles=Path(old['future_roles'])
    rows=json.loads(old_roles.read_text());assert len(rows)==1620
    forbidden=seed_census(pins);assert len(forbidden)==gate['forbidden_seed_count']
    assert gate['future_seed_namespace']==SEED_NAMESPACE
    seeds=fresh_seeds([r['chain']['chain_id'] for r in rows],forbidden,SEED_NAMESPACE)
    assert len(seeds)==1620 and set(seeds.values()).isdisjoint(forbidden)
    root=a.output.resolve();root.mkdir(exist_ok=False);(root/'models').mkdir()
    models={};roles=[];groups=defaultdict(list)
    for row in rows:
        original=row['chain'];group=original['effective_input_group']+'-'+original['prior_label']
        source=Path(original['program']);assert sha(source)==original['program_sha256']
        text=source.read_text();changed=transform(text);program=root/'models'/(group+'.hs')
        if group not in models:
            program.write_text(changed);assert restore(program.read_text())==text
            assert 'probeValuesV6' not in changed and 'logEncoderChecks' not in changed
            models[group]=dict(model_input_identity=group,family=original['family'],prior_label=original['prior_label'],
                source_program=str(source),source_program_sha256=sha(source),program=str(program),program_sha256=sha(program),
                reversible_scalar_logger_edit_verified=True,probability_model_changed=False,
                project_scalar_schema=SCHEMA,native_cjson_parameters=True,explicit_nonfinite_context_and_parameter_tags=True)
            bind(pins,source);bind(pins,program)
        else:assert program.read_text()==changed
        chain=dict(original,chain_id=group+'-scalar-cjson-v6-chain'+str(original['chain']),
                   seed=seeds[original['chain_id']],program=str(program),program_sha256=sha(program))
        roles.append(dict(chain=chain,source_v5_chain_id=original['chain_id'],source_v5_seed=original['seed'],
            seed_namespace=SEED_NAMESPACE,native_execution_launched=False,scientific_eligibility=False,posterior_qualified=False))
        groups[group].append(chain)
    assert len(models)==len(groups)==405 and len(roles)==1620
    assert all(len(v)==4 and {c['chain'] for c in v}=={1,2,3,4} for v in groups.values())
    assert Counter(r['chain']['prior_label'] for r in roles)=={'broad':540,'centered':540,'package':540}
    assert len({r['chain']['effective_input_group'] for r in roles})==135
    assert len({alias for r in roles for alias in r['chain']['original_configuration_ids']})==324
    role_path=root/'future_roles.json';model_path=root/'model_manifest.json'
    write_json(role_path,roles);write_json(model_path,[models[k] for k in sorted(models)])
    for q in [a.validation,a.execution,a.transport,old_path,old_roles,HELPER,role_path,model_path,Path(__file__),
              Path('scripts/baliphy_scalar_json_logger_v6c.py'),Path('scripts/read_baliphy_scalar_json_v6b.py'),
              Path('scripts/check_baliphy_scalar_json_logger_v6_v9.py')]:bind(pins,q)
    verify(pins)
    result=dict(status='prepared_full_future_scalar_json_v6_models_not_launched',checked_utc=datetime.now(timezone.utc).isoformat(),
        models=405,roles=1620,model_quartets=405,effective_inputs=135,original_configuration_aliases=324,
        model_manifest=str(model_path),future_roles=str(role_path),seed_namespace=SEED_NAMESPACE,forbidden_seed_count=len(forbidden),
        project_scalar_schema=SCHEMA,exactly_reversible_from_v5=True,explicit_nonfinite_context_and_parameter_tags=True,
        probability_model_changed=False,software_fixture_values_included=False,full_source_grid_reread=True,
        existing_sampler_audit_compatible=False,full_grid_v6_startup_or_sampler_audit_qualified=False,
        production_horizon=None,execution_commands_prepared=False,native_execution_launched=False,
        scientific_eligibility=False,posterior_qualified=False,gpu=False,new_cost_usd=0,pins=pins,
        scope='All405V5 models with exactly reversible scalar rendering changes and1620distinct fresh role seeds. '
        'Native property/sequence/category logs, TSV and model/priors/initializer are preserved. Scalar records '
        'carry a declared V6 schema plus explicit context/parameter nonfinite and literal-null tags. '
        'Existing full-sampler scalar audits require V6-aware integration before any future launch; '
        'full-grid startup/runtime/mixing/likelihood/root/model and posterior adequacy remain unqualified. '
        'No historical artifact repair, old source edits, restart, production horizon, GPU use or paid resources.')
    with a.receipt.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='pins'},indent=2))


if __name__=='__main__':main()
