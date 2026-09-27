#!/usr/bin/env python3
"""Retain both endpoint mappings and prespecified metadata-caliper sensitivities for all edges."""
import argparse,csv,json,time
from collections import Counter
from pathlib import Path
import psutil
from background_match_covariates import differences,BANDS
from run_ortholog_pair_guide_comparison import sha

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for path,h in plan['pins'].items():assert sha(path)==h,path
    verify()
    if 'producer' in plan:
        dep=plan['producer']
        while psutil.pid_exists(dep['pid']):
            try:
                proc=psutil.Process(dep['pid'])
                if abs(proc.create_time()-dep['created'])>.01 or proc.status()==psutil.STATUS_ZOMBIE:break
                assert proc.cmdline()==dep['cmdline']
            except psutil.NoSuchProcess:break
            time.sleep(30)
    verify();root=Path(plan['source']);rp=root/'receipt.json';r=json.loads(rp.read_text());proof=Path(plan['readback']);audit=json.loads(proof.read_text())
    assert audit['status']=='passed_full_background_match_graph_readback' and audit['producer_receipt_sha256']==sha(rp)
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    def nodes(name):return {r['node_id']:r for line in open(root/name) for r in [json.loads(line)]}
    targets=nodes('target_nodes.jsonl');backgrounds=nodes('background_nodes.jsonl');counts=Counter();n=0
    out=Path(plan['output']);out.mkdir(exist_ok=False)
    with open(root/'edges.tsv') as src,(out/'edge_covariates.tsv').open('w') as dest:
        writer=None
        for edge in csv.DictReader(src,delimiter='\t'):
            target=targets[edge['target_id']];background=backgrounds[edge['background_id']];values=differences(target,background);row=dict(edge,**values)
            if writer is None:writer=csv.DictWriter(dest,list(row),delimiter='\t',lineterminator='\n');writer.writeheader()
            writer.writerow(row);n+=1;prefix=target['guide']+'|'+edge['policy'];counts[prefix+'|edges']+=1
            for band in BANDS:counts[prefix+'|within_'+band]+=values['within_'+band]
            for kind in ['genes','models','sequences']:counts[prefix+'|shared_'+kind]+=values['shared_'+kind]>0
            if n%200000==0:print('Covariate edges',n,flush=True)
    assert n==r['edges'];verify()
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    result=dict(status='complete_background_match_covariates_pending_readback',plan_sha256=ph,graph_receipt_sha256=sha(rp),graph_readback_sha256=sha(proof),edges=n,bands=BANDS,counts=dict(counts),artifacts={'edge_covariates.tsv':sha(out/'edge_covariates.tsv')},scope='All candidate edges, both endpoint mappings. Each caliper band requires one mapping to satisfy all endpoint-wise maximum length-ratio, mean-pLDDT and low-confidence-fraction differences jointly; component minima are never combined across mappings. Bands are descriptive sensitivity settings, not validated reliability cutoffs or selected matches. Shared genes/models/sequences retained as dependency flags. No structural outcome used, confidence balancing may change the estimand and does not remove prediction circularity.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
