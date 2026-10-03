#!/usr/bin/env python3
"""Describe the complete closed initial-horizon resource survey; no inference."""
import argparse
from collections import Counter,defaultdict
import csv
from datetime import datetime,timezone
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--completion',type=Path,default=Path('metadata/baliphy_horizon_resources_completed_20261003.json'))
    p.add_argument('--prefix',type=Path,default=Path('docs/figures/baliphy_horizon_resources_20261003'))
    p.add_argument('--table',type=Path,default=Path('docs/tables/baliphy_horizon_resources_20261003.tsv'))
    p.add_argument('--receipt',type=Path,default=Path('metadata/baliphy_horizon_resource_figure_20261003.json'))
    a=p.parse_args();complete=json.loads(a.completion.read_text())
    assert complete['status']=='complete_verified_full_baliphy_horizon_resource_inventory'
    bindings={str(a.completion):sha(a.completion)}
    for key in ['full_hash_archive','producer_receipt','independent_readback']:
        assert sha(complete[key])==complete[key+'_sha256'];bindings[complete[key]]=complete[key+'_sha256']
    producer=json.loads(Path(complete['producer_receipt']).read_text());root=Path(complete['producer_receipt']).parent
    for name,d in producer['artifacts'].items():assert sha(root/name)==d;bindings[str(root/name)]=d
    chains=[json.loads(line) for line in (root/'chains.jsonl').open()]
    quartets=[json.loads(line) for line in (root/'quartets.jsonl').open()]
    assert len(chains)==1620 and len(quartets)==405
    grouped=defaultdict(list)
    for chain in chains:grouped[chain['model_input_identity']].append(chain)
    events=[];table=[];hours=defaultdict(float);counts=Counter();failed=[]
    for group in sorted(grouped):
        members=grouped[group];first=members[0];success=[r for r in members if r['selected_status']=='all_saved_alignments_and_candidate_nodes_checked']
        family=first['family'];prior=first['prior_label'];seconds=sum(r['selected_attempt']['elapsed_worker_seconds'] for r in success)
        hours[family,prior]+=seconds/3600
        for r in members:
            attempt=r['selected_attempt']
            for parameter,value in attempt['trace']['parameters'].items():
                events.extend(dict(chain_id=r['chain_id'],parameter=parameter,**v) for v in value['nonfinite'])
            if r in success:
                counts[prior,'nonfinite']+=any(v['nonfinite'] for v in attempt['trace']['parameters'].values())
                counts[prior,'warning']+=bool(attempt['allocation_warning_lines'])
            else:failed.append(r)
        table.append(dict(group=group,family=family,prior=prior,effective_input_group=first['effective_input_group'],
            chain_ids=json.dumps([r['chain_id'] for r in members]),checked_chains=len(success),failed_chains=4-len(success),
            successful_elapsed_worker_seconds=seconds,successful_elapsed_worker_hours=seconds/3600,
            successful_nonfinite_parameter_chains=sum(any(v['nonfinite'] for v in r['selected_attempt']['trace']['parameters'].values()) for r in success),
            successful_allocation_warning_chains=sum(bool(r['selected_attempt']['allocation_warning_lines']) for r in success)))
    assert len(table)==405 and sum(r['checked_chains'] for r in table)==1618 and len(failed)==2
    assert len(events)==complete['attempt_diagnostics']['nonfinite_parameter_observations']==336
    assert sum(counts[p,'nonfinite'] for p in ['broad','centered','package'])==134
    assert sum(counts[p,'warning'] for p in ['broad','centered','package'])==4
    assert all(e['parameter']=='ASRV.Gamma:alpha' and 1<=e['iteration']<=14 for e in events)
    assert abs(sum(hours.values())-complete['observed_successful_elapsed_worker_seconds']/3600)<1e-9
    artifacts=[a.table,a.receipt,*[Path(str(a.prefix)+suffix) for suffix in ['.png','.svg','.pdf']]]
    assert not any(q.exists() for q in artifacts),'Figure exports are immutable; use a new prefix/table/receipt'
    a.table.parent.mkdir(parents=True,exist_ok=True);a.prefix.parent.mkdir(parents=True,exist_ok=True)
    with a.table.open('x') as f:
        writer=csv.DictWriter(f,fieldnames=list(table[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(table)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    colors={'broad':'#386cb0','centered':'#f1a340','package':'#1b9e77'};priors=list(colors);families=sorted({r['family'] for r in chains});y=np.arange(len(families));left=np.zeros(len(y))
    for prior in priors:
        values=np.asarray([hours[f,prior] for f in families]);axes[0,0].barh(y,values,left=left,color=colors[prior],label=prior);left+=values
    axes[0,0].set_yticks(y,families);axes[0,0].invert_yaxis();axes[0,0].set_xlabel('Successful worker wall time (hours)')
    axes[0,0].set_title('A. All 13 families; 1,618 completed chains');axes[0,0].legend(frameon=False,fontsize=8)
    x=np.arange(3)
    for offset,kind,label,color in [(-.18,'nonfinite','Nonfinite alpha, iterations 1–14','#b2182b'),(.18,'warning','Allocation warning','#8073ac')]:
        values=[counts[p,kind] for p in priors];bars=axes[0,1].bar(x+offset,values,width=.36,label=label,color=color);axes[0,1].bar_label(bars,padding=3)
    axes[0,1].set_xticks(x,priors);axes[0,1].set_ylabel('Distinct completed chains');axes[0,1].set_title('B. Early numerical and allocation reviews')
    axes[0,1].legend(frameon=False,fontsize=8);axes[0,1].set_ylim(0,max(counts.values())*1.25)
    count=Counter(e['iteration'] for e in events);iterations=list(range(1,15))
    axes[1,0].bar(iterations,[count[i] for i in iterations],color='#b2182b');axes[1,0].set_xticks(iterations)
    axes[1,0].set_xlabel('Logged iteration');axes[1,0].set_ylabel('Nonfinite alpha observations');axes[1,0].set_title('C. All 336 selected-attempt events')
    for row in failed:
        trace=row['selected_attempt']['trace'];log=Path(row['selected_attempt']['scalar_log'])
        assert sha(log)==row['selected_attempt']['scalar_log_sha256'];bindings[str(log)]=sha(log)
        with log.open() as f:records=list(csv.DictReader(f,delimiter='\t'))
        dataset='whole' if '-whole-' in row['original_configuration_ids'][0] else 'domain'
        axes[1,1].plot([int(r['iter']) for r in records],[float(r['|A|']) for r in records],marker='o',label=f'{dataset}, {row["prior_label"]}')
        assert int(records[-1]['iter'])==trace['last']['iteration']
    axes[1,1].set_yscale('log');axes[1,1].set_xticks([0,1,2]);axes[1,1].set_xlabel('Logged iteration');axes[1,1].set_ylabel('Logged alignment columns')
    axes[1,1].set_title('D. Two failed OG0000972 chains');axes[1,1].legend(frameon=False,fontsize=8)
    fig.suptitle('Initial ancestral sampling: measured resources and unresolved numerical behavior',fontsize=13)
    for suffix in ['.png','.svg','.pdf']:fig.savefig(str(a.prefix)+suffix,dpi=180)
    plt.close(fig)
    for path,d in bindings.items():assert sha(path)==d,path
    result=dict(status='complete_descriptive_full_horizon_resource_figure',checked_utc=datetime.now(timezone.utc).isoformat(),
        source_hashes={**bindings,__file__:sha(__file__)},full_chains=1620,full_quartets=405,families=len(families),
        successful_nonfinite_parameter_chains=134,successful_allocation_warning_chains=4,nonfinite_alpha_observations=336,
        nonfinite_iteration_range=[1,14],nonfinite_observations_after_first_burn_in_cutoff=0,
        artifacts={str(q):sha(q) for q in artifacts if q!=a.receipt},scientific_eligibility=False,
        scope='All405dependent model quartets retained; time is worker wall duration, not CPU hours or a finish forecast. '
        'Warnings and nonfinite initial tokens are retained, not proof of a particular cause or stationary-posterior corruption. '
        'Both selected failures plotted from complete partial raw logs; no extrapolated memory bound or accepted posterior/root/model.')
    with a.receipt.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']}))


if __name__=='__main__':main()
