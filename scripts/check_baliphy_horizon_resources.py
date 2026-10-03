#!/usr/bin/env python3
"""Check full-grid resource accounting, failed traces and false serialized exports."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from baliphy_horizon_resource_inventory import aggregate,collect,trace_summary
from readback_baliphy_horizon_resources import verify_exports
from run_ortholog_pair_guide_comparison import sha


def write(p,value):
    Path(p).parent.mkdir(parents=True,exist_ok=True)
    Path(p).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def fixture(root):
    root=Path(root);bound={}
    def bound_write(p,value):write(p,value);bound[str(p)]=sha(p)
    (root/'input').write_text('synthetic fixture, not a biological alignment\n')
    inp=str(root/'input');chain_list=[];rows=[];details={};groups={}
    header='iter\tRS07:rate\tRS07:meanLength\tASRV.Gamma:alpha\t|A|\t|indels|\n'
    good=header+''.join(f'{i}\t0.02\t2\t1\t100\t2\n' for i in range(1001))
    partial=header+'0\t0.02\t2\t1\t100\t2\n1\t4e17\t55000\tInfinity\t491177\t1351287\n'
    for i in range(1620):
        g=i//12;prior=['broad','centered','package'][(i%12)//4];number=i%4+1
        model=f'group{g:03d}-{prior}';cid=f'{model}-chain{number}'
        chain=dict(chain_id=cid,seed=20000+i,chain=number,effective_input_group=f'group{g:03d}',
            prior_label=prior,original_configuration_ids=[f'family{g}-alias1',f'family{g}-alias2'],
            alignment=inp,tree=inp,program=inp,alignment_sha256=sha(inp),tree_sha256=sha(inp),program_sha256=sha(inp))
        chain_list.append(chain)
        def attempt(kind,success,warnings):
            base=root/kind/cid;folder=base/'attempt-0001';folder.mkdir(parents=True)
            command=['/usr/bin/prlimit','--as='+str((12 if kind=='original' else 48)*2**30),
                '--','/fixture/bali-phy','--seed',str(chain['seed']),'run',inp,'--iterations','1000']
            configuration=dict(command=command,seed=chain['seed'],model_input_identity=model)
            bound_write(base/'configuration.json',configuration)
            encoded=json.dumps(configuration,sort_keys=True,separators=(',',':'),allow_nan=False)
            bound_write(folder/'command.json',command)
            (folder/'stderr.log').write_text('Allocation failed in sample_tri_multi!  Proceeding.\n'*warnings+('' if success else 'std::bad_alloc\n'))
            (folder/'C1.log').write_text(good if success else partial)
            (folder/'C1.P1.fastas').write_text('SIZE ONLY; NEVER DECODED AS A SAMPLE\n')
            for name in ['stderr.log','C1.log','C1.P1.fastas']:bound[str(folder/name)]=sha(folder/name)
            receipt=dict(configuration_sha256=hashlib.sha256(encoded.encode()).hexdigest(),
                exit_code=0 if success else 1,status='exited_zero_pending_scientific_validation' if success else 'failed',
                elapsed_seconds=100+i if success else 10+i,
                artifacts={p.name:sha(p) for p in folder.iterdir()})
            bound_write(folder/'receipt.json',receipt)
            return dict(chain_id=cid,status='all_saved_alignments_and_candidate_nodes_checked' if success else 'failed',
                receipt=str(folder/'receipt.json'),receipt_sha256=sha(folder/'receipt.json'))
        original=attempt('original',i not in [0,4,8],0)
        selected=attempt('recovery',i==0,1) if i in [0,4,8] else original
        rows.append(dict(chain=chain,model_input_identity=model,original_disposition=original,
            selected_disposition=selected,scientific_eligibility=False))
        groups.setdefault(model,dict(chain_ids=[]))['chain_ids'].append(cid)
        details[cid]={}
        if selected['status']=='all_saved_alignments_and_candidate_nodes_checked':
            details[cid]=dict(selected_attempt_receipt_sha256=selected['receipt_sha256'],
                native_files={'C1.P1.fastas':{'bytes_at_observation':Path(selected['receipt']).parent.joinpath('C1.P1.fastas').stat().st_size}})
    overlay=root/'overlay.json';bound_write(overlay,dict(rows=rows,original_samples_concatenated=False))
    rec=root/'recovery.json';bound_write(rec,dict(overlay=str(overlay),overlay_sha256=sha(overlay)))
    prod=root/'producer.json';bound_write(prod,dict(groups=groups))
    archive=root/'archive.json';write(archive,dict(source_hashes=bound,services=[{'scope':'synthetic'},{'scope':'synthetic'}]))
    complete=root/'completion.json';write(complete,dict(status='complete_verified_full_baliphy_recovery_diagnostics',
        full_hash_archive=str(archive),full_hash_archive_sha256=sha(archive),bound_source_hashes=len(bound),
        producer_receipt=str(prod),producer_receipt_sha256=sha(prod)))
    cp=root/'chains.json';write(cp,chain_list)
    det=root/'details.json';write(det,details)
    inv=root/'inventory.json';write(inv,dict(full_chain_details=str(det),full_chain_details_sha256=sha(det)))
    scalar_archive=root/'scalar_archive.json';write(scalar_archive,dict(services=[{'scope':'synthetic'},{'scope':'synthetic'}],source_hashes={}))
    scalar=root/'scalar.json';write(scalar,dict(status='complete_verified_full_independent_baliphy_scalar_comparison',
        full_hash_archive=str(scalar_archive),full_hash_archive_sha256=sha(scalar_archive),bound_source_hashes=0,
        quartets_passing_every_scalar={'250':0,'500':0},quartets_passing_every_length_scalar={'250':0,'500':0}))
    return dict(pins={},diagnostic_completion=str(complete),recovery_completion=str(rec),chain_inputs=str(cp),
        native_inventory=str(inv),scalar_completion=str(scalar),proposed_iterations=10000,seed_namespace='fixture')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    rejected=[]
    with tempfile.TemporaryDirectory() as temporary:
        root=Path(temporary);plan=fixture(root);data,bindings=collect(plan);summary=data['summary']
        assert summary['complete_quartets']==403 and summary['selected_failed_chains']==2
        assert summary['attempt_diagnostics']['bad_alloc_attempts']==5
        assert summary['attempt_diagnostics']['nonfinite_parameter_observations']==5
        assert summary['selected_successful_allocation_warning_chains']==1
        assert summary['sampler_review_quartets']==3
        assert all(r['retained_alignment_samples']==[750,500] and r['expected_saved_alignments']==1001 for r in data['proposed'])
        new={r['fresh_seed'] for r in data['proposed']};old={r['seed'] for r in data['chains']};assert len(new)==1620 and new.isdisjoint(old)
        assert len(bindings)>6000 and all(r['memory_peak_unavailable'] for r in data['attempts'])
        exports=root/'exports';exports.mkdir()
        def export(value):
            for name in ['chains','attempts','quartets','proposed']:
                (exports/(name+'.jsonl')).write_text(''.join(json.dumps(r,allow_nan=False)+'\n' for r in value[name]))
        export(data);verify_exports(exports,data)
        alterations={
            'omit_failed_chain':lambda d:d['chains'].pop(4),
            'omit_recovery_attempt':lambda d:d['attempts'].pop(),
            'erase_allocation_warning':lambda d:next(r for r in d['attempts'] if r['allocation_warning_lines']).update(allocation_warning_lines=0),
            'erase_invalid_parameter':lambda d:d['attempts'][0]['trace'].update(numerical_parameter_review=False),
            'change_runtime':lambda d:d['chains'][0]['selected_attempt'].update(elapsed_worker_seconds=1),
            'change_seed':lambda d:d['proposed'][0].update(fresh_seed=d['chains'][0]['seed']),
            'permit_launch':lambda d:d['proposed'][0].update(production_launch_allowed=True),
            'change_alias':lambda d:d['proposed'][0].update(original_configuration_ids=[]),
            'omit_quartet':lambda d:d['quartets'].pop(),
            'claim_memory_peak':lambda d:d['attempts'][0].update(memory_peak_unavailable=False)}
        for name,change in alterations.items():
            value=deepcopy(data);change(value);export(value)
            try:verify_exports(exports,data)
            except AssertionError:rejected.append(name)
            else:raise AssertionError(name)
        export(data)
        stage=dict(plan,output=str(root/'stage'))
        plan_path=root/'plan.json';write(plan_path,stage)
        for script in ['scripts/prepare_baliphy_horizon_resources.py','scripts/readback_baliphy_horizon_resources.py']:
            subprocess.run([sys.executable,script,'--plan',str(plan_path)],check=True,stdout=subprocess.DEVNULL)
        for script in ['scripts/prepare_baliphy_horizon_resources.py','scripts/readback_baliphy_horizon_resources.py']:
            trial=subprocess.run([sys.executable,script,'--plan',str(plan_path)],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
            assert trial.returncode != 0,'Completed stage restarted: '+script
        for name,text,complete in [('iteration_gap','iter\tRS07:rate\tRS07:meanLength\tASRV.Gamma:alpha\t|A|\t|indels|\n1\t1\t1\t1\t1\t1\n',False),
            ('incomplete_success','iter\tRS07:rate\tRS07:meanLength\tASRV.Gamma:alpha\t|A|\t|indels|\n0\t1\t1\t1\t1\t1\n',True)]:
            q=root/'bad_trace';q.write_text(text)
            try:trace_summary(q,complete,1000)
            except AssertionError:rejected.append(name)
            else:raise AssertionError(name)
    paths=[__file__,'scripts/baliphy_horizon_resource_inventory.py','scripts/prepare_baliphy_horizon_resources.py','scripts/readback_baliphy_horizon_resources.py']
    result=dict(status='passed_full_baliphy_horizon_resource_inventory_software_contracts',
        full_chains=1620,full_quartets=405,effective_input_groups=135,attempts=1623,
        complete_quartets=403,failed_chains=2,fresh_seeds=1620,
        full_producer_and_serialized_reader_exercised=True,completed_restarts_refused=True,
        false_exports_and_malformed_traces_rejected=rejected,source_hashes={p:sha(p) for p in paths},
        scientific_eligibility=False,scope='Full synthetic census/serialization, not native inference, a biological pilot or real completion journals.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result))


if __name__=='__main__':main()
