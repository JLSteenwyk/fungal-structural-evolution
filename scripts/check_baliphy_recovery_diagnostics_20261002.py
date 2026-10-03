#!/usr/bin/env python3
"""Check full original diagnostic accounting and synthetic recovered inputs."""
import argparse
import copy
import csv
from datetime import datetime,timezone
import json
from pathlib import Path
import tempfile
import numpy as np
from audit_baliphy_memory_recovery_20261002_v2 import build_overlay,CHECKED
from baliphy_recovery_diagnostic_sources import dispatch,summarize
from prepare_baliphy_recovery_diagnostics_20261002 import scalar_and_length,categorical,write
from readback_baliphy_recovery_diagnostics_20261002 import new_group_inputs
from run_ortholog_pair_guide_comparison import sha


PLANS=dict(scalar='metadata/baliphy_independent_chain_diagnostic_plan_20260927.json',
    states='metadata/baliphy_full_anchored_states_plan_20260928_v2.json',
    length='metadata/baliphy_candidate_length_diagnostic_plan_20260927.json',
    categorical='metadata/ancestral_categorical_queue_plan_20260928.json')


def baseline(jobs,original,reports):
    indexed={r['chain_id']:r for r in original}
    rows={j['chain']['chain_id']:dict(chain=j['chain'],selected_disposition=indexed[j['chain']['chain_id']],
        original_disposition=indexed[j['chain']['chain_id']]) for j in jobs}
    states={};groups={};routes={}
    for cid,row in rows.items():
        old=reports['states']['chains'][cid]
        states[cid]=dict(status='reused_closed_original_state_trace',receipt=old['receipt'],receipt_sha256=old['receipt_sha256']) if 'receipt' in old else dict(
            status='unresolved_failed_native_chain_retained',native_disposition=row['selected_disposition'])
    for group in reports['scalar']['groups']:
        ids=sorted(j['chain']['chain_id'] for j in jobs if j['config']['model_input_identity']==group)
        ready=all(rows[c]['selected_disposition']['status']==CHECKED for c in ids)
        origin='reuse_closed_original_unchanged_quartet' if ready else 'unresolved_failed_native_chain_retained'
        routes[group]=dict(status=origin,chain_ids=ids)
        info=dict(status='complete_scalar_length_category_screens_not_posterior_qualification' if ready else origin,
                  origin=origin,chain_ids=ids,scientific_eligibility=False)
        if ready:
            for kind in ['scalar','length','categorical']:
                old=reports[kind]['groups'][group];info[kind]=dict(receipt=old['receipt'],receipt_sha256=old['receipt_sha256'])
        groups[group]=info
    source=dict(rows=rows,routes=routes,reports=reports,summary=dict(selected_checked_chains=1617,complete_quartets=402))
    summary=summarize(groups,states,source,{})
    post=json.loads(Path('metadata/baliphy_full_postprocessing_completed_20261002_v2.json').read_text())
    horizon=json.loads(Path('metadata/baliphy_initial_horizon_completed_20261002.json').read_text())
    for k in ['state_observations','per_chain_anchor_occurrences','unanchored_residue_observations','length_status_counts',
        'quartets_passing_every_length_scalar','categorical_coordinate_counts','categorical_pattern_counts','coordinate_status_counts','pattern_status_counts']:
        assert summary[k]==post[k],k
    assert summary['scalar_status_counts']==horizon['scalar_variable_counts'] and summary['quartets_passing_every_scalar']==horizon['quartets_passing_every_scalar']
    # A source-hash-consistent foreign report cannot replace an unchanged group.
    names=[g for g,r in groups.items() if 'scalar' in r];bad=copy.deepcopy(groups)
    bad[names[0]]['scalar']=copy.deepcopy(groups[names[1]]['scalar'])
    try:summarize(bad,states,source,{})
    except AssertionError:pass
    else:raise AssertionError('Foreign original scalar report reused')
    return summary


def fixture(root,python):
    rng=np.random.default_rng(793);alignment=root/'alignment.faa';alignment.write_text('>t0\nAA\n>t1\nAAA\n>t2\nAAA\n')
    nodes=['n0','n1','n2','n3'];tips=[dict(tip='t0',length=2),dict(tip='t1',length=3),dict(tip='t2',length=3)]
    coordinates=dict(alphabet='ACDEFGHIKLMNPQRSTVWYX-',nodes=nodes,tips=tips,input_alignment=str(alignment),input_alignment_sha256=sha(alignment))
    fixed=['iter','scale','scale1','scale*|T|','scale1*|T|','|T|'];variables=['prior','likelihood','posterior','ASRV.Gamma:alpha','RS07:rate','RS07:meanLength']+['F:pi[%d]'%i for i in range(20)]
    members=[];states={}
    for ordinal in range(4):
        cid='synthetic-chain%d'%ordinal;folder=root/cid;folder.mkdir();log=folder/'C1.log'
        with log.open('x') as f:
            writer=csv.DictWriter(f,fixed+variables,delimiter='\t',lineterminator='\n');writer.writeheader()
            for iteration in range(1001):writer.writerow(dict(zip(fixed+variables,[iteration]+[1]*5+rng.normal(size=len(variables)).tolist())))
        samples=[dict(iteration=i,level=level,source_node=nodes[level],runtime_node=nodes[level],ungapped_length=level+20+(i//10)%3) for i in range(0,1001,10) for level in range(4)]
        audit=folder/'audit.json';write(audit,dict(scalar_log=str(log),scalar_log_sha256=sha(log),candidate_samples=samples))
        attempt=folder/'native_receipt.json';write(attempt,dict(scope='Synthetic input only'))
        selected=dict(receipt=str(attempt),receipt_sha256=sha(attempt),sample_audit=str(audit),sample_audit_sha256=sha(audit))
        chain=dict(chain_id=cid,seed=2000+ordinal,alignment=str(alignment),alignment_sha256=sha(alignment))
        members.append(dict(chain=chain,selected_disposition=selected,model_input_identity='synthetic-group'))
        reportdir=folder/'state-report';target=reportdir/cid;target.mkdir(parents=True)
        values=rng.integers(0,22,size=(101,4,8),dtype=np.uint8);iterations=np.arange(0,1001,10);free=np.zeros((101,4),dtype=np.int32)
        counts={f'counts_after_{cutoff}':np.stack([np.count_nonzero(values[iterations>cutoff]==state,axis=0) for state in range(22)],axis=-1).astype(np.uint16) for cutoff in [250,500]}
        arrays=target/'states.npz';np.savez_compressed(arrays,states=values,iterations=iterations,unanchored_residue_counts=free,**counts)
        cp=target/'coordinates.json';write(cp,coordinates);rp=reportdir/'receipt.json'
        write(rp,dict(summaries=[dict(source_audit_sha256=sha(audit),artifacts={str(arrays):sha(arrays),str(cp):sha(cp)})]))
        states[cid]=dict(receipt=str(rp),receipt_sha256=sha(rp))
    output=root/'new-quartet';output.mkdir();reports,chains=scalar_and_length(dict(diagnostic_python=python),members,output)
    reports['categorical']=categorical(dict(diagnostic_python=python),members,states,chains,output)
    info=dict(chain_ids=[r['chain']['chain_id'] for r in members],**reports)
    source=dict(rows={r['chain']['chain_id']:r for r in members});new_group_inputs(info,source,states,{})
    rejected=[]
    # These mutations target manifest identity, sampling horizon, selected
    # native length traces, and complete assembled state/free arrays.
    for kind in ['scalar','length']:
        receipt=json.loads(Path(info[kind]['receipt']).read_text());dp=Path(receipt['outputs']['250']['path'])
        report=json.loads(dp.read_text());mp=Path(next(p for p in report['pins'] if Path(p).name=='manifest-250.json'))
        original=mp.read_bytes()
        for change in ['seed','grid']:
            altered=json.loads(original)
            if change=='seed':altered['chains'][0]['seed']+=1
            else:altered['expected_iterations']=altered['expected_iterations'][:-1]
            mp.write_text(json.dumps(altered))
            try:new_group_inputs(info,source,states,{})
            except AssertionError:rejected.append(kind+'_'+change)
            else:raise AssertionError('Changed manifest accepted')
            mp.write_bytes(original)
    category=json.loads(Path(info['categorical']['receipt']).read_text());manifest=json.loads(Path(category['manifest']).read_text());fp=Path(manifest['arrays'])
    original=fp.read_bytes()
    for key in ['values','unanchored_residue_counts']:
        with np.load(fp,allow_pickle=False) as saved:arrays={k:saved[k].copy() for k in saved.files}
        arrays[key].flat[0]+=1;np.savez_compressed(fp,**arrays)
        try:new_group_inputs(info,source,states,{})
        except AssertionError:rejected.append('assembled_'+key)
        else:raise AssertionError('Changed state-array assembly accepted')
        fp.write_bytes(original)
    return rejected


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    source=json.loads(Path('metadata/baliphy_independent_chain_plan_20260927.json').read_text());jobs=json.loads(Path(source['jobs']).read_text())
    original=json.loads((Path(source['output'])/'receipt.json').read_text())
    retry_plan=json.loads(Path('metadata/baliphy_memory_recovery_plan_20261002.json').read_text());retry=json.loads((Path(retry_plan['output'])/'receipt.json').read_text())
    rows,quartets,summary=build_overlay(jobs,original['chains'],retry['recovery_chains'],retry_plan['failed_chain_ids'])
    routes=dispatch(dict(rows=rows,quartets=quartets));counts={status:sum(r['status']==status for r in routes.values()) for status in {r['status'] for r in routes.values()}}
    assert counts==dict(reuse_closed_original_unchanged_quartet=402,compute_complete_recovery_quartet=1,unresolved_failed_native_chain_retained=2)
    reports={k:json.loads((Path(json.loads(Path(v).read_text())['output'])/'receipt.json').read_text()) for k,v in PLANS.items()}
    old_summary=baseline(jobs,original['chains'],reports)
    python=json.loads(Path(PLANS['scalar']).read_text())['python']
    with tempfile.TemporaryDirectory(prefix='fungal-recovery-diagnostic-contract-') as d:rejected=fixture(Path(d),python)
    paths=[Path(__file__),Path('scripts/baliphy_recovery_diagnostic_sources.py'),Path('scripts/prepare_baliphy_recovery_diagnostics_20261002.py'),Path('scripts/readback_baliphy_recovery_diagnostics_20261002.py'),*map(Path,PLANS.values())]
    result=dict(status='passed_full_baliphy_recovery_diagnostic_accounting_and_input_contracts',checked_utc=datetime.now(timezone.utc).isoformat(),
        full_original_grid=summary,quartet_routes=counts,complete_original_diagnostic_summary_reproduced=old_summary,
        foreign_unchanged_quartet_report_rejected=True,synthetic_whole_quartet_input_readback_passed=True,
        malformed_new_inputs_rejected=rejected,source_hashes={str(p):sha(p) for p in paths},
        scope='Full original402-quartet accounting reproduced, full405-quartet recovery routes retained, and synthetic four-chain manifests/length/state-array assemblies checked. No native sampling, biological pilot, new production diagnostic metric claim or posterior qualification.')
    write(a.output,result);print(json.dumps(result))


if __name__=='__main__':main()
