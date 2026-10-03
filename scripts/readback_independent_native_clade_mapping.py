#!/usr/bin/env python3
"""Rebuild all source/runtime clades and compare every serialized node pair."""
import argparse
import fcntl
import json
from pathlib import Path

from independent_native_clade_mapping import load,replay,summarize
from prepare_independent_native_clade_mapping import marker,artifacts
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def run(path,output):
    plan=json.loads(path.read_text());root=Path(plan['output']);lock=(root/'readback.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);assert not (root/'readback.json').exists(),'Completed original reader cannot restart'
    source,bindings=load(plan,path);initial=dict(bindings);rp=root/'receipt.json';receipt=json.loads(rp.read_text());bind(bindings,rp)
    assert receipt['status']=='complete_full_independent_native_clade_mapping_pending_readback'
    assert receipt['plan_sha256']==sha(path) and receipt['source_hashes']==initial and receipt['scientific_eligibility'] is False
    for name,d in receipt['artifacts'].items():bind(bindings,root/name,d)
    verify(bindings);stage=marker(path,source);assert json.loads((root/'stage_plan.json').read_text())==stage
    entries=json.loads((root/'chain_manifest.json').read_text());assert [r['chain_id'] for r in entries]==sorted(source['details'])
    results={}
    for number,(cid,entry) in enumerate(zip(sorted(source['details']),entries),1):
        assert entry['path']=='chains/'+cid+'.json';cp=root/entry['path'];assert sha(cp)==entry['sha256']
        expected=replay(source,cid);assert json.loads(cp.read_text())==dict(stage=stage,chain_id=cid,result=expected)
        results[cid]=expected;print('readback_independent_native_clade_chains',number,'/',len(entries),flush=True)
    summary=summarize(results,source,plan);assert all(receipt[k]==v for k,v in summary.items())
    assert artifacts(root,entries)==receipt['artifacts'];verify(bindings)
    result=dict(status='passed_full_independent_native_clade_mapping_serialized_readback',plan_sha256=sha(path),
        producer_receipt_sha256=sha(rp),**summary,source_hashes=bindings,scientific_eligibility=False,
        third_independent_tree_decoder=False,scope=plan['scope'])
    with Path(output).open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(summary),flush=True);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.plan,a.output)
