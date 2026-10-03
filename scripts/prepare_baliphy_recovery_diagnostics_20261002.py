#!/usr/bin/env python3
"""Finish full scalar/length/category diagnostics with whole-attempt recovery.

Closed unchanged quartets are reused; every source setting and unresolved chain
is retained. New state projections and reports use the pinned original tools.
"""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys
import time
import numpy as np
from Bio import SeqIO
from baliphy_recovery_diagnostic_sources import load,document,summarize
from audit_baliphy_memory_recovery_20261002_v2 import CHECKED
from prepare_ancestral_length_diagnostics import length_rows
from read_ancestral_state_quartet import validate_arrays
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def write(path,value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def command(cmd,log):
    # These are bounded postprocessing children, never native sampling jobs.
    with log.open('x') as out:
        subprocess.run(cmd,check=True,stdout=out,stderr=subprocess.STDOUT)
    write(log.with_suffix('.command.json'),dict(command=cmd,exit_code=0,log_sha256=sha(log)))


def new_states(plan,source,root):
    states={}
    for cid,row in source['rows'].items():
        selected=row['selected_disposition']
        if selected['status']!=CHECKED:
            states[cid]=dict(status='unresolved_failed_native_chain_retained',native_disposition=selected)
        elif row['selected_disposition']==row['original_disposition']:
            old=source['reports']['states']['chains'][cid]
            assert old['status']=='state_trace_complete_not_posterior_qualification'
            states[cid]=dict(status='reused_closed_original_state_trace',receipt=old['receipt'],receipt_sha256=old['receipt_sha256'])
        else:
            folder=root/'state-extraction'/cid;folder.mkdir(parents=True,exist_ok=False)
            snapshot=folder/'snapshot.json'
            write(snapshot,dict(status='verified_frozen_terminal_chain_subset',rows=[dict(chain=cid,
                seed=row['chain']['seed'],audit=selected['sample_audit'],audit_sha256=selected['sample_audit_sha256'])]))
            output=folder/'traces'
            command([sys.executable,'scripts/prepare_ancestral_state_traces.py','--snapshot',str(snapshot),
                '--inputs',plan['chain_inputs'],'--output',str(output)],folder/'extract.log')
            rp=output/'receipt.json';receipt=json.loads(rp.read_text())
            assert receipt['status']=='complete_independently_checked_anchored_state_traces' and receipt['chains']==1
            states[cid]=dict(status='new_recovery_whole_chain_state_trace',receipt=str(rp),receipt_sha256=sha(rp))
    return states


def scalar_and_length(plan,members,folder):
    scalar_chains=[];length_chains=[];audits=[];nodes=None;headers=[]
    expected=list(range(0,1001,10));length_folder=folder/'length';length_folder.mkdir()
    for row in members:
        chain=row['chain'];selected=row['selected_disposition'];audit=json.loads(Path(selected['sample_audit']).read_text())
        base=dict(chain_id=chain['chain_id'],seed=chain['seed'],model_input_identity=row['model_input_identity'])
        scalar_chains.append(dict(**base,log=audit['scalar_log'],log_sha256=audit['scalar_log_sha256']))
        with Path(audit['scalar_log']).open() as f:headers.append(next(csv.reader(f,delimiter='\t')))
        traces,current=length_rows(audit['candidate_samples'],expected)
        if nodes is not None:assert current==nodes
        nodes=current
        fp=length_folder/(chain['chain_id']+'.tsv')
        with fp.open('x') as f:
            writer=csv.DictWriter(f,list(traces[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(traces)
        length_chains.append(dict(**base,log=str(fp),log_sha256=sha(fp)))
        audits.append(selected)
    assert all(header==headers[0] for header in headers)
    fixed={'iter','scale','scale1','scale*|T|','scale1*|T|','|T|'};assert fixed<=set(headers[0])
    variables=[name for name in headers[0] if name not in fixed]
    assert {'prior','likelihood','posterior','ASRV.Gamma:alpha','RS07:rate','RS07:meanLength'}<=set(variables)
    assert len([v for v in variables if v.startswith('F:pi[')])==20
    result={}
    for kind,chains,grid,names in [('scalar',scalar_chains,list(range(1001)),variables),
                                   ('length',length_chains,expected,['length_level%d'%i for i in range(4)])]:
        target=folder/kind;target.mkdir(exist_ok=True);outputs={}
        for cutoff in [250,500]:
            mp=target/('manifest-%d.json'%cutoff);write(mp,dict(chains=chains,variables=names,expected_iterations=grid,discard_through_iteration=cutoff))
            destination=target/('discard-%d'%cutoff)
            command([plan['diagnostic_python'],'scripts/ancestral_chain_diagnostics.py','--manifest',str(mp),
                     '--output',str(destination)],target/('diagnose-%d.log'%cutoff))
            dp=destination/'diagnostics.json';outputs[str(cutoff)]=dict(path=str(dp),sha256=sha(dp))
        rp=target/'receipt.json';write(rp,dict(status='complete_recovery_'+kind+'_screens_not_posterior_qualification',
            group=members[0]['model_input_identity'],source_attempts=audits,outputs=outputs))
        result[kind]=dict(receipt=str(rp),receipt_sha256=sha(rp))
    return result,scalar_chains


def categorical(plan,members,states,chains,folder):
    values=[];free=[];common=None;evidence={};expected=list(range(0,1001,10))
    for row in members:
        chain=row['chain'];info=states[chain['chain_id']];rp=Path(info['receipt'])
        assert sha(rp)==info['receipt_sha256'];report=json.loads(rp.read_text());assert len(report['summaries'])==1
        summary=report['summaries'][0];assert summary['source_audit_sha256']==row['selected_disposition']['sample_audit_sha256']
        ap=rp.parent/chain['chain_id']/'states.npz';cp=ap.parent/'coordinates.json'
        assert sha(ap)==summary['artifacts'][str(ap)] and sha(cp)==summary['artifacts'][str(cp)]
        coordinates=json.loads(cp.read_text());assert coordinates['input_alignment_sha256']==chain['alignment_sha256']==sha(chain['alignment'])
        records=list(SeqIO.parse(chain['alignment'],'fasta'));assert len(records)==len({r.id for r in records})
        assert coordinates['tips']==[dict(tip=r.id,length=len(str(r.seq).replace('-',''))) for r in sorted(records,key=lambda r:r.id)]
        audit=json.loads(Path(row['selected_disposition']['sample_audit']).read_text())
        assert coordinates['nodes']==sorted({r['source_node'] for r in audit['candidate_samples']})
        if common is None:common=coordinates
        assert common==coordinates
        with np.load(ap,allow_pickle=False) as saved:
            values.append(validate_arrays(saved,coordinates,expected));free.append(saved['unanchored_residue_counts'].copy())
        for p in [rp,ap,cp,Path(row['selected_disposition']['sample_audit']),Path(row['selected_disposition']['receipt'])]:evidence[str(p)]=sha(p)
    target=folder/'categorical';target.mkdir();arrays=target/'quartet.npz'
    np.savez_compressed(arrays,values=np.stack(values),iterations=np.asarray(expected),unanchored_residue_counts=np.stack(free))
    manifest=target/'manifest.json';write(manifest,dict(status='provenance_checked_state_quartet',chains=chains,
        coordinates=common,expected_iterations=expected,arrays=str(arrays),arrays_sha256=sha(arrays),evidence=evidence))
    output=target/'reports'
    command([plan['diagnostic_python'],'scripts/report_ancestral_state_quartet.py','--manifest',str(manifest),
             '--output',str(output)],target/'report.log')
    rp=output/'receipt.json';report=json.loads(rp.read_text())
    assert report['status']=='both_cutoff_categorical_reports_complete_not_posterior_qualification'
    receipt=target/'receipt.json';write(receipt,dict(status='verified_quartet_categorical_reports_complete_not_posterior_qualification',
        group=members[0]['model_input_identity'],report=str(rp),report_sha256=sha(rp),manifest=str(manifest),manifest_sha256=sha(manifest)))
    return dict(receipt=str(receipt),receipt_sha256=sha(receipt))


def run(path):
    started=time.monotonic();plan=json.loads(path.read_text());source,bindings=load(plan,path)
    root=Path(plan['output']);root.mkdir(exist_ok=False);states=new_states(plan,source,root);groups={}
    for group,route in source['routes'].items():
        info=dict(origin=route['status'],chain_ids=route['chain_ids'],scientific_eligibility=False)
        if route['status']=='unresolved_failed_native_chain_retained':info['status']=route['status']
        else:
            info['status']='complete_scalar_length_category_screens_not_posterior_qualification'
            if route['status']=='reuse_closed_original_unchanged_quartet':
                for kind in ['scalar','length','categorical']:
                    old=source['reports'][kind]['groups'][group]
                    info[kind]=dict(receipt=old['receipt'],receipt_sha256=old['receipt_sha256'])
            else:
                folder=root/'groups'/group;folder.mkdir(parents=True,exist_ok=False)
                members=[source['rows'][cid] for cid in route['chain_ids']]
                reports,chains=scalar_and_length(plan,members,folder);info.update(reports)
                info['categorical']=categorical(plan,members,states,chains,folder)
            print('recovery_diagnostic_quartet',group,route['status'],flush=True)
        groups[group]=info
    summary=summarize(groups,states,source,bindings);verify(bindings)
    # Bind all freshly written manifests, logs, arrays and reports, including
    # completed postprocessing commands. No large data enter Git history.
    artifacts={str(p.relative_to(root)):sha(p) for p in root.rglob('*') if p.is_file()}
    receipt=dict(status='complete_full_baliphy_recovery_diagnostics_pending_accounting_readback',
        plan_sha256=sha(path),groups=groups,states=states,**summary,artifacts=artifacts,
        source_hashes=bindings,elapsed_seconds=time.monotonic()-started,scientific_eligibility=False,scope=plan['scope'])
    write(root/'receipt.json',receipt);print(json.dumps(summary),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True)
    run(p.parse_args().plan)
