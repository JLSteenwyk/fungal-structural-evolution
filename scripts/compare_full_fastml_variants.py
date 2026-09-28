#!/usr/bin/env python3
"""Compare every paired FastML fit and native marginal without clipping."""
import argparse
from collections import defaultdict,Counter
import json
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
from Bio import Phylo
from ancestral_chain_attempt import sha,write_json


def node_splits(path):
    tree=Phylo.read(path,'newick')
    if tree.root.name is None:tree.root.name=tree.root.comment
    mapping={}
    for node in tree.find_clades():
        assert node.name is not None and node.name not in mapping
        mapping[node.name]=tuple(sorted(x.name for x in node.get_terminals()))
    return mapping


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    plan=json.loads(args.plan.read_text());producer=json.loads(Path(plan['producer_plan']).read_text())
    for p in [plan,producer]:
        for path,h in p['pins'].items():assert sha(path)==h,path
    for unit in ['fungal-full-fastml-variants-20260928','fungal-full-fastml-roundoff-readback-20260928']:
        state=dict(l.split('=',1) for l in subprocess.check_output(['systemctl','--user','show',unit,'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
        assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0'),(unit,state)
    root=Path(plan['output']);receipt=json.loads((root/'receipt.json').read_text())
    assert receipt['plan_sha256']==sha(args.plan) and not receipt['unresolved']
    assert len(receipt['readbacks'])==len(producer['jobs'])==312
    groups=defaultdict(dict);artifacts=0
    for item in producer['jobs']:
        entry=receipt['readbacks'][item['id']];assert sha(entry['path'])==entry['sha256']
        r=json.loads(Path(entry['path']).read_text());assert r['variant']==item['variant'] and r['job_id']==item['job']['job_id']
        for path,h in r['evidence'].items():assert sha(path)==h;artifacts+=1
        native=json.loads((Path(producer['output'])/item['id']/'readback.json').read_text())
        assert native['job']==item['job'] and native['variant']==item['variant']
        folder=None
        if item['job']['character_count']:
            rp=Path(native['attempt_receipt']);assert sha(rp)==native['attempt_receipt_sha256']
            ar=json.loads(rp.read_text());assert ar['exit_code']==0
            assert json.loads((rp.parent.parent/'configuration.json').read_text())==item['config']
            for name,h in ar['artifacts'].items():assert sha(rp.parent/name)==h;artifacts+=1
            assert r['status']=='integrity_and_independent_replay_complete_not_fit_qualification'
            folder=rp.parent/'RESULTS'
        else:assert r['status']=='no_coded_characters'
        assert item['variant'] not in groups[r['job_id']]
        groups[r['job_id']][item['variant']]=(item['job'],r,folder)
    args.output.mkdir(parents=True,exist_ok=False)
    rows=[];total_rows=0;crossings=0;changed=0;maximum=0.;empty=[]
    for job_id,pair in sorted(groups.items()):
        assert set(pair)=={'precision-only','precision-cache-refresh'}
        j,a,af=pair['precision-only'];k,b,bf=pair['precision-cache-refresh'];assert j==k
        if not j['character_count']:
            empty.append(job_id);continue
        assert node_splits(af/'TheTree.INodes.ph')==node_splits(bf/'TheTree.INodes.ph')
        tables=[]
        for folder,r in [(af,a),(bf,b)]:
            table=pd.read_csv(folder/'AncestralReconstructPosterior.txt',sep='\t')
            assert len(table)==r['probability_rows'] and table.State.eq(1).all()
            assert not table.duplicated(['Node','POS','State']).any()
            table=table.set_index(['Node','POS','State']).sort_index()
            assert np.isfinite(table.Prob).all()
            tables.append(table)
        assert tables[0].index.equals(tables[1].index)
        x=tables[0].Prob.to_numpy();y=tables[1].Prob.to_numpy();delta=np.abs(y-x)
        flip=(x>=.5)!=(y>=.5)
        row=dict(job_id=job_id,character_count=j['character_count'],probability_rows=len(x),
            maximum_native_probability_difference=float(delta.max()),native_rows_changed=int(np.count_nonzero(delta)),
            native_rows_difference_over_1e_minus6=int(np.count_nonzero(delta>1e-6)),
            native_rows_difference_over_1e_minus3=int(np.count_nonzero(delta>1e-3)),presence_threshold_crossings=int(flip.sum()),
            maximum_fitted_parameter_absolute_difference=max(abs(a['fitted_parameters'][p]-b['fitted_parameters'][p]) for p in ['alpha','gain','loss']),
            precision_only_replay_minus_reported=a['likelihood_difference'],cache_refresh_replay_minus_reported=b['likelihood_difference'],
            reported_likelihood_change=b['reported_log_likelihood']-a['reported_log_likelihood'],
            replay_likelihood_change=b['replay_log_likelihood']-a['replay_log_likelihood'],
            precision_only_maximum_replay_probability_error=a['maximum_probability_difference'],
            cache_refresh_maximum_replay_probability_error=b['maximum_probability_difference'])
        rows.append(row);total_rows+=len(x);crossings+=int(flip.sum());changed+=int(np.count_nonzero(delta));maximum=max(maximum,float(delta.max()))
    assert len(groups)==156 and len(rows)==153 and len(empty)==3 and total_rows==8058340
    pd.DataFrame(rows).to_csv(args.output/'paired_fit_comparison.tsv',sep='\t',index=False)
    result=dict(status='complete_full_paired_native_probability_comparison',paired_inputs=len(groups),nonempty_pairs=len(rows),empty_inputs=empty,
        paired_probability_rows=total_rows,native_rows_changed=changed,maximum_native_probability_difference=maximum,presence_threshold_crossings=crossings,
        fits_with_native_probability_change=sum(r['native_rows_changed']>0 for r in rows),
        fits_with_parameter_change=sum(r['maximum_fitted_parameter_absolute_difference']>0 for r in rows),
        precision_only_maximum_absolute_replay_likelihood_difference=max(abs(r['precision_only_replay_minus_reported']) for r in rows),
        cache_refresh_maximum_absolute_replay_likelihood_difference=max(abs(r['cache_refresh_replay_minus_reported']) for r in rows),
        fits_cache_refresh_absolute_likelihood_difference_over_1e_minus5=sum(abs(r['cache_refresh_replay_minus_reported'])>1e-5 for r in rows),
        artifact_hashes_checked=artifacts,source_readback_receipt_sha256=sha(root/'receipt.json'),
        pins={str(p):sha(p) for p in [args.plan,Path(plan['producer_plan']),Path(__file__)]},
        artifacts={'paired_fit_comparison.tsv':sha(args.output/'paired_fit_comparison.tsv')},
        scope='All paired native marginals compared at matching node-descendant identities and character positions. Thresholds are descriptive, not fit-acceptance gates. No raw clipping. Cache refresh sensitivity does not establish optimizer convergence, correct ascertainment or joint indel/sequence posterior adequacy.')
    write_json(args.output/'receipt.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':main()
