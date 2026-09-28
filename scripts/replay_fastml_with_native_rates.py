#!/usr/bin/env python3
"""Separate fixed-parameter replay using rates from the frozen native library."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from Bio import Phylo,SeqIO
from ancestral_chain_attempt import sha,write_json
from replay_fastml_indel_probabilities import infer


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);args=ap.parse_args()
    plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        assert sha(args.plan)==ph
        for path,h in plan['pins'].items():assert sha(path)==h,path
    verify()
    producer=json.loads(Path(plan['producer_plan']).read_text());readroot=Path(plan['readback'])
    rr=json.loads((readroot/'receipt.json').read_text());assert not rr['unresolved']
    rate_root=Path(plan['rates']);rate_receipt=json.loads((rate_root/'receipt.json').read_text())
    for name,h in rate_receipt['artifacts'].items():assert sha(rate_root/name)==h
    rates={r['id']:r for r in map(json.loads,(rate_root/'rate_comparison.jsonl').read_text().splitlines())}
    assert len(rates)==306
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);rows=[];empty=[]
    for item in producer['jobs']:
        key=item['id'];entry=rr['readbacks'][key];assert sha(entry['path'])==entry['sha256']
        original=json.loads(Path(entry['path']).read_text())
        for path,h in original['evidence'].items():assert sha(path)==h
        if not item['job']['character_count']:
            assert original['status']=='no_coded_characters';empty.append(key);continue
        assert original['status']=='integrity_and_independent_replay_complete_not_fit_qualification'
        source=Path(producer['output'])/key/'readback.json';native=json.loads(source.read_text())
        rp=Path(native['attempt_receipt']);assert sha(rp)==native['attempt_receipt_sha256'];attempt=json.loads(rp.read_text());assert attempt['exit_code']==0
        for name,h in attempt['artifacts'].items():assert sha(rp.parent/name)==h
        folder=rp.parent/'RESULTS';r=rates[key];assert r['alpha']==original['fitted_parameters']['alpha']
        tree=Phylo.read(folder/'TheTree.INodes.ph','newick');assert tree.root.name is None and tree.root.comment=='N1';tree.root.name='N1'
        job=item['job'];observed={s.id:str(s.seq) for s in SeqIO.parse(job['characters'],'fasta')}
        assert len(observed)==job['proteins'] and set(observed)=={n.name for n in tree.get_terminals()}
        assert {len(v) for v in observed.values()}=={job['character_count']}
        pars=original['fitted_parameters'];rate=np.array(r['native_rates']);assert rate.shape==(4,) and np.isfinite(rate).all() and (rate>0).all()
        post,ll=infer(tree,{k:v+'0' for k,v in observed.items()},pars['gain'],pars['loss'],rate)
        assert np.isfinite(ll).all() and ll[-1]<0
        corrected=float((ll[:-1]-np.log(-np.expm1(ll[-1]))).sum());maximum=0.;seen=set()
        with (folder/'AncestralReconstructPosterior.txt').open() as f:
            for row in csv.DictReader(f,delimiter='\t'):
                node,pos,value=row['Node'],int(row['POS']),float(row['Prob'])
                assert row['State']=='1' and node in post and 1<=pos<=job['character_count'] and np.isfinite(value)
                assert (node,pos) not in seen;seen.add((node,pos))
                difference=abs(value-float(post[node][pos-1]));assert np.isfinite(difference);maximum=max(maximum,difference)
        assert len(seen)==original['probability_rows']
        result=dict(id=key,job_id=job['job_id'],variant=item['variant'],probability_rows=len(seen),
            native_rate_replay_log_likelihood=corrected,reported_log_likelihood=original['reported_log_likelihood'],
            native_rate_likelihood_difference=corrected-original['reported_log_likelihood'],
            scipy_rate_likelihood_difference=original['likelihood_difference'],
            native_rate_maximum_probability_difference=maximum,scipy_rate_maximum_probability_difference=original['maximum_probability_difference'],
            source_readback_sha256=entry['sha256'],source_attempt_sha256=sha(rp),native_rates=rate.tolist(),raw_values_clipped=False)
        write_json(out/(key+'.json'),result);rows.append(result)
        print('Native-rate replays',len(rows),'/306',flush=True)
    verify();assert len(rows)==306 and len(empty)==6
    write_json(out/'receipt.json',dict(status='complete_native_rate_fixed_parameter_replay_not_fit_qualification',plan_sha256=ph,fits=len(rows),empty_inputs=empty,
        maximum_absolute_likelihood_difference=max(abs(r['native_rate_likelihood_difference']) for r in rows),
        maximum_probability_difference=max(r['native_rate_maximum_probability_difference'] for r in rows),
        artifacts={p.name:sha(p) for p in out.iterdir()},
        scope='Uses native gamma discretization with independent pruning implementation. Preserves SciPy-rate replay and raw probabilities. Does not qualify optimizer stationarity or mask-conditional ascertainment, joint indel uncertainty, or ancestral ensembles.'))


if __name__=='__main__':main()
