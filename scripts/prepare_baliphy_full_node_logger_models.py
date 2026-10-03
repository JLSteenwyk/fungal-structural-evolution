#!/usr/bin/env python3
"""Prepare the complete future logger source/seed grid without launching inference."""
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha, write_json
from baliphy_full_node_logger import BEFORE, AFTER, transform, fresh_seeds
from run_baliphy_reference_preflight import verify


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--validation',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--receipt',type=Path,required=True)
    args=parser.parse_args();gate=json.loads(args.validation.read_text())
    assert gate['status']=='passed_full_node_logger_static_and_native_software_qualification'
    assert gate['full_programs_checked']==405 and gate['full_roles_checked']==1620
    assert gate['exact_same_seed_scalar_and_alignment_pairs']==3 and gate['native_fixture_runs']==6
    assert gate['full_node_saved_frames']==9 and gate['full_node_ancestor_frames']==36
    assert gate['original_tip_records_preserved'] and len(gate['malformed_ancestor_records_rejected'])==18
    verify({'pins':gate['source_hashes']})
    source=Path('metadata/baliphy_reference_sampler_qualification_plan_20261003.json')
    native=json.loads(source.read_text());verify(native);jobs=json.loads(Path(native['jobs']).read_text())
    assert len(jobs)==1620
    root=args.output.resolve();root.mkdir(exist_ok=False);(root/'models').mkdir()
    used={j['chain']['seed'] for j in jobs}|{j['source_seed'] for j in jobs}|{20267001,20267002,20267003}
    namespace='fungal-full-node-logger-future-20261003-v1'
    seeds=fresh_seeds([j['chain']['chain_id'] for j in jobs],used,namespace)
    models={};roles=[];groups=defaultdict(list);pins=dict(native['pins'])
    for job in jobs:
        old=job['chain'];group=old['effective_input_group']+'-'+old['prior_label'];path=Path(old['program'])
        assert sha(path)==old['program_sha256'];original=path.read_text();changed=transform(original)
        target=root/'models'/(group+'.hs')
        if group not in models:
            target.write_text(changed)
            assert target.read_text().replace(AFTER,BEFORE)==original
            models[group]=dict(model_input_identity=group,prior_label=old['prior_label'],family=old['family'],
                source_program=str(path),source_program_sha256=sha(path),program=str(target),program_sha256=sha(target),
                reversible_logger_edit_verified=True,probability_model_changed=False)
            pins[str(target)]=sha(target)
        else:assert target.read_text()==changed
        chain=dict(old,chain_id=group+'-full-node-future-chain'+str(old['chain']),seed=seeds[old['chain_id']],
                   program=str(target),program_sha256=sha(target))
        roles.append(dict(chain=chain,source_chain_id=old['chain_id'],earlier_source_seed=job['source_seed'],
            current_short_sampler_seed=old['seed'],seed_namespace=namespace,native_execution_launched=False,
            posterior_qualified=False,scientific_eligibility=False))
        groups[group].append(chain)
    assert len(models)==len(groups)==405 and len(roles)==1620
    assert len({r['chain']['seed'] for r in roles})==len({r['chain']['chain_id'] for r in roles})==1620
    assert all(len(v)==4 and {c['chain'] for c in v}=={1,2,3,4} for v in groups.values())
    assert Counter(r['chain']['prior_label'] for r in roles)=={'broad':540,'centered':540,'package':540}
    assert len({r['chain']['effective_input_group'] for r in roles})==135
    assert len({a for r in roles for a in r['chain']['original_configuration_ids']})==324
    role_path=root/'future_roles.json';model_path=root/'model_manifest.json'
    write_json(role_path,roles);write_json(model_path,[models[k] for k in sorted(models)])
    for p in [source,args.validation,Path(__file__),Path('scripts/baliphy_full_node_logger.py'),role_path,model_path]:pins[str(p)]=sha(p)
    verify({'pins':pins})
    result=dict(status='prepared_full_future_node_logger_sources_and_seed_namespace_not_launched',
        checked_utc=datetime.now(timezone.utc).isoformat(),source_plan=str(source),source_plan_sha256=sha(source),
        models=405,roles=1620,model_quartets=405,effective_inputs=135,original_configuration_aliases=324,
        full_source_grid_reread=True,model_manifest=str(model_path),future_roles=str(role_path),seed_namespace=namespace,
        disjoint_from_prior_real_and_logger_fixture_seeds=True,forbidden_seed_count=len(used),pins=pins,
        generated_source_bytes=sum(Path(r['program']).stat().st_size for r in models.values()),
        planning_state='source_and_seed_preparation_only',production_horizon=None,execution_commands_prepared=False,
        native_execution_launched=False,gpu=False,new_cost_usd=0,scientific_eligibility=False,posterior_qualified=False,
        scope='All405futuremodelprograms x4newdisjointseedroles,135effectiveinputs,324aliases. One exact reversible pure logger label-view edit; full source/manifest/seed inventory checked. Existing sampler/source/attempts untouched. No new production command/horizon, source startup closure, native full-grid logging qualification, memory/mixing/root/model/likelihood/predictor qualification or scientific acceptance. Longer production remains gated on complete original sampler/resource/output closures and a separately justified resource/uncertainty plan.')
    with args.receipt.open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='pins'}),flush=True)


if __name__=='__main__':main()
