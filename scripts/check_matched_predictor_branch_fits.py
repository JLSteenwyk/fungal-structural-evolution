#!/usr/bin/env python3
"""Qualify seven native model roles and complete 931-role workflow accounting."""
import argparse
import copy
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np

from ancestral_chain_attempt import sha
from matched_predictor_branch_fits import ROLES,config,execute,load
from matched_predictor_branch_inputs import verify
from readback_matched_predictor_branch_fits import independent_decode,check_point,summarize


def reject(action):
    try:action()
    except (AssertionError,ValueError,KeyError,FileNotFoundError):return
    raise AssertionError('Altered native export accepted')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True);a=p.parse_args()
    assert not a.output.exists() and not a.receipt.exists();a.output.mkdir(parents=True)
    binary=Path(shutil.which('iqtree3')).resolve();models={k:str(Path('data/structural_models/garg-hochberg-v3')/('Q.3Di.'+k)) for k in ['AF','LLM']}
    model_receipt=Path('metadata/3di_substitution_model_receipt.json');r=json.loads(model_receipt.read_text())
    assert r['status']=='complete_published_3di_model_validation'
    assert all(sha(v)==r['artifacts'][Path(v).name] for v in models.values())
    limits=dict(native_address_space_gib=2,native_cpu_seconds=300,native_wall_seconds=600,native_per_file_limit_mib=8)
    plan=dict(executable=str(binary),models=models,resources=limits)
    key=hashlib.sha256(b'matched-predictor-software-fixture-v1').hexdigest();folder=a.output/'fixture_inputs'/'inputs'/key;folder.mkdir(parents=True)
    rng=np.random.default_rng(202610031);letters=list('ARNDCQEGHILKMFPSTWYV');names=list('abcde')
    base=rng.choice(letters,240);state_base=rng.choice(letters,240)
    content={k:{} for k in ['aa','AlphaFold','ESMFold']}
    for i,t in enumerate(names):
        aa=base.copy();af=state_base.copy();esm=state_base.copy()
        for seq,amount in [(aa,20+i*7),(af,16+i*8),(esm,19+i*9)]:
            positions=rng.choice(len(seq),amount,replace=False);seq[positions]=rng.choice(letters,amount)
        for label,seq in [('aa',aa),('AlphaFold',af),('ESMFold',esm)]:content[label][t]=''.join(seq)
    for label,data in content.items():(folder/(label+'.faa')).write_text(''.join('>'+t+'\n'+data[t]+'\n' for t in names))
    (folder/'topology.nwk').write_text('((a:0.1,b:0.1):0.1,c:0.1,(d:0.1,e:0.1):0.1);\n')
    sc=dict(input_id=key,marker='synthetic_native_fixture',taxa=names,columns=list(range(1,241)),required_observed=50,
        internal_splits=['0x3','0x18'],alignment_sha256={k:sha(folder/(k+'.faa')) for k in content},topology_sha256=sha(folder/'topology.nwk'))
    (folder/'config.json').write_text(json.dumps(sc,indent=2)+'\n')
    out=a.output/'native_fixture';(out/'roles').mkdir(parents=True);positions={t:i for i,t in enumerate(names)};rows=[]
    for role in ROLES:
        row=execute(plan,folder.parent.parent,out,key,sc,role,positions)
        assert row['status']=='native_point_output_integrity_checked_not_model_qualified',row
        assert row['point_estimate'] is not None
        independent=independent_decode(Path(row['native_receipt']).parent,sc,positions);check_point(row['point_estimate'],independent)
        rows.append(row)
    own=['matched_predictor_branch_fits','run_matched_predictor_branch_fits','readback_matched_predictor_branch_fits',
         'check_matched_predictor_branch_fits','ancestral_chain_attempt','run_paired_marker_fits']
    pins={f'scripts/{n}.py':sha(f'scripts/{n}.py') for n in own}
    pins.update({str(p):sha(p) for p in [binary,model_receipt,*map(Path,models.values())]})
    source_plan=dict(input_completion='metadata/matched_predictor_branch_inputs_completed_20261003_v1.json',
        input_plan='metadata/matched_predictor_branch_inputs_plan_20261003_v1.json',pins=pins)
    spp=a.output/'source_plan.json';spp.write_text(json.dumps(source_plan,indent=2)+'\n')
    real,bindings=load(source_plan,spp);templates={row['role']:row for row in rows};mock=[];configs=[]
    for real_key,real_cfg in sorted(real['configs'].items()):
        for role in ROLES:
            cfg=config(plan,real['root'],real_key,real_cfg,role);configs.append(cfg)
            assert '-te' in cfg['command'] and '-keep-ident' in cfg['command'] and cfg['command'][cfg['command'].index('-T')+1]=='1'
            row=copy.deepcopy(templates[role[0]]);row['input_id']=real_key;mock.append(row)
        for model in ['af','af_empirical','llm']:
            pair=[c for c in configs if c['input_id']==real_key and c['role'].endswith('_'+model)]
            assert len(pair)==2 and pair[0]['seed']==pair[1]['seed']
    assert len(mock)==len(configs)==931
    original=summarize(mock,real['configs'],8750);assert original['intact_inputs']==133
    first_keys=sorted(real['configs'])[:2]
    for real_key in first_keys:
        row=next(r for r in mock if r['input_id']==real_key and r['role']=='aa')
        row.update(status='native_unsuccessful_retained',point_estimate=None,exit_code=1,native_status='failed')
    mixed=summarize(mock,real['configs'],8750);assert mixed['unresolved_inputs']==2 and mixed['native_status_counts']['native_unsuccessful_retained']==2
    serialized=a.output/'full_mock_roles.json';serialized.write_text(json.dumps(mock,indent=2)+'\n')
    assert summarize(json.loads(serialized.read_text()),real['configs'],8750)==mixed
    negatives=[]
    for label in ['missing_role','duplicate_role','unknown_role']:
        bad=copy.deepcopy(mock)
        if label=='missing_role':bad.pop()
        elif label=='duplicate_role':bad[-1]=bad[0]
        else:bad[0]['role']='unknown'
        reject(lambda:summarize(bad,real['configs'],8750));negatives.append(label)
    baseline=rows[0]['point_estimate']
    for label in ['negative_branch','changed_split','missing_branch','changed_likelihood','wrong_unit']:
        bad=copy.deepcopy(baseline)
        if label=='negative_branch':bad['branches'][0]['length']=-1
        elif label=='changed_split':bad['branches'][0]['split_mask_hex']='0x999'
        elif label=='missing_branch':bad['branches'].pop()
        elif label=='changed_likelihood':bad['log_likelihood_reported']+=1
        else:bad['branch_unit']='angstrom'
        reject(lambda:check_point(bad,baseline));negatives.append(label)
    bindings.update(pins);verify(bindings)
    result=dict(status='passed_full_matched_predictor_native_fit_software_contracts',checked_utc=datetime.now(timezone.utc).isoformat(),
        full_real_inputs=133,full_real_native_role_configs=931,synthetic_native_roles=7,synthetic_native_point_branch_values=sum(len(r['point_estimate']['branches']) for r in rows),
        independent_dendropy_tree_and_report_readback=True,full_mock_role_serialization_checked=True,
        artificial_failed_roles_retained=2,artificial_unresolved_inputs=2,malformed_cases_rejected=negatives,
        source_hashes=bindings,artifacts={str(p):sha(p) for p in a.output.rglob('*') if p.is_file()},scientific_eligibility=False,
        scope='Seven bounded synthetic native IQTree fits cover AA plus both predictors under three declared structural models; independently decoded complete tree/report branch values. All133real input/seven-role configs and931-row mocked serialization retain two artificial failed roles/unresolved inputs and reject eight altered outputs. Fixed topology, matching predictor seeds, identical source observation masks, one thread and exact input/model/binary pins checked. This is software qualification within full-scope study, not a biological pilot, experimental predictor error or accepted evolutionary effect; full931nativefits and uncertainty/model/framework controls remain pending.')
    with a.receipt.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']},indent=2),flush=True)


if __name__=='__main__':main()
