#!/usr/bin/env python3
"""Reconcile the full original grid and recovered quartet's report inputs.

This is accounting and input-array readback, not an independent reimplementation
of ArviZ diagnostics, native alignment parsing or posterior qualification.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from baliphy_recovery_diagnostic_sources import load,summarize
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def new_group_inputs(info,source,states,bindings):
    members={cid:source['rows'][cid] for cid in info['chain_ids']}
    audits={cid:json.loads(Path(row['selected_disposition']['sample_audit']).read_text()) for cid,row in members.items()}
    for kind in ['scalar','length']:
        receipt=json.loads(Path(info[kind]['receipt']).read_text())
        assert receipt['status']=='complete_recovery_'+kind+'_screens_not_posterior_qualification'
        assert receipt['source_attempts']==[members[cid]['selected_disposition'] for cid in info['chain_ids']]
        for cutoff,value in receipt['outputs'].items():
            dp=Path(value['path']);report=json.loads(dp.read_text())
            manifests=[p for p in report['pins'] if Path(p).name=='manifest-'+cutoff+'.json']
            assert len(manifests)==1;mp=Path(manifests[0]);bind(bindings,mp,report['pins'][str(mp)])
            manifest=json.loads(mp.read_text());chains=manifest['chains']
            assert [r['chain_id'] for r in chains]==info['chain_ids'] and len(chains)==4
            assert manifest['discard_through_iteration']==int(cutoff)
            assert set(manifest['variables'])==set(report['variables'])
            expected=list(range(1001)) if kind=='scalar' else list(range(0,1001,10))
            assert manifest['expected_iterations']==expected
            for chain in chains:
                cid=chain['chain_id'];row=members[cid]
                assert chain['seed']==row['chain']['seed'] and chain['model_input_identity']==row['model_input_identity']
                bind(bindings,chain['log'],chain['log_sha256'])
                if kind=='scalar':
                    assert chain['log']==audits[cid]['scalar_log'] and chain['log_sha256']==audits[cid]['scalar_log_sha256']
                    with Path(chain['log']).open() as f:header=next(csv.reader(f,delimiter='\t'))
                    assert manifest['variables']==[v for v in header if v not in {'iter','scale','scale1','scale*|T|','scale1*|T|','|T|'}]
                else:
                    assert manifest['variables']==['length_level%d'%i for i in range(4)]
                    indexed={(int(s['iteration']),int(s['level'])):int(s['ungapped_length']) for s in audits[cid]['candidate_samples']}
                    assert len(indexed)==404
                    with Path(chain['log']).open() as f:traces=list(csv.DictReader(f,delimiter='\t'))
                    assert [int(r['iter']) for r in traces]==expected
                    for trace in traces:
                        assert set(trace)=={'iter',*manifest['variables']}
                        for level in range(4):assert int(trace['length_level%d'%level])==indexed[int(trace['iter']),level]
    category=json.loads(Path(info['categorical']['receipt']).read_text())
    manifest=json.loads(Path(category['manifest']).read_text())
    with np.load(manifest['arrays'],allow_pickle=False) as quartet:
        assert np.array_equal(quartet['iterations'],np.arange(0,1001,10))
        assert quartet['values'].shape[0]==quartet['unanchored_residue_counts'].shape[0]==4
        assert [r['chain_id'] for r in manifest['chains']]==info['chain_ids']
        for i,cid in enumerate(info['chain_ids']):
            rp=Path(states[cid]['receipt']);report=json.loads(rp.read_text());summary=report['summaries'][0]
            fp=rp.parent/cid/'states.npz';bind(bindings,fp,summary['artifacts'][str(fp)])
            with np.load(fp,allow_pickle=False) as saved:
                assert np.array_equal(quartet['values'][i],saved['states'])
                assert np.array_equal(quartet['unanchored_residue_counts'][i],saved['unanchored_residue_counts'])
            chain=manifest['chains'][i]
            assert chain['seed']==members[cid]['chain']['seed']
            assert chain['log']==audits[cid]['scalar_log'] and chain['log_sha256']==audits[cid]['scalar_log_sha256']


def run(path,output):
    plan=json.loads(path.read_text());source,bindings=load(plan,path);original=dict(bindings)
    root=Path(plan['output']);rp=root/'receipt.json';receipt=json.loads(rp.read_text());bind(bindings,rp)
    assert receipt['status']=='complete_full_baliphy_recovery_diagnostics_pending_accounting_readback'
    assert receipt['plan_sha256']==sha(path) and receipt['scientific_eligibility'] is False
    for p,d in original.items():assert receipt['source_hashes'][p]==d
    for p,d in receipt['source_hashes'].items():bind(bindings,p,d)
    for name,d in receipt['artifacts'].items():bind(bindings,root/name,d)
    for group,info in receipt['groups'].items():
        if info['origin']=='compute_complete_recovery_quartet':new_group_inputs(info,source,receipt['states'],bindings)
    summary=summarize(receipt['groups'],receipt['states'],source,bindings)
    assert all(receipt[k]==v for k,v in summary.items())
    verify(bindings)
    result=dict(status='passed_full_baliphy_recovery_diagnostic_accounting_and_new_quartet_inputs',
        plan_sha256=sha(path),producer_receipt_sha256=sha(rp),**summary,
        source_hashes=bindings,scientific_eligibility=False,scope=plan['scope'])
    with Path(output).open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(summary),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.plan,a.output)
