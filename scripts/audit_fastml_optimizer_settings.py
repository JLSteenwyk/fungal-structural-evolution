#!/usr/bin/env python3
"""Inventory effective optimization limits and native termination diagnostics."""
import argparse
from collections import Counter
import json
from pathlib import Path
import re
import pandas as pd
from ancestral_chain_attempt import sha,write_json


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    plan=json.loads(args.plan.read_text());root=Path(plan['output']);rr=json.loads((root/'receipt.json').read_text())
    assert rr['plan_sha256']==sha(args.plan)
    for path,h in plan['pins'].items():assert sha(path)==h
    rows=[];empty=[];source_bindings={}
    for item in plan['jobs']:
        path=root/item['id']/'readback.json';assert sha(path)==rr['readbacks'][str(path)]
        r=json.loads(path.read_text());assert r['job']==item['job'] and r['variant']==item['variant']
        if not item['job']['character_count']:
            assert r['status']=='no_coded_characters';empty.append(item['id']);continue
        rp=Path(r['attempt_receipt']);assert sha(rp)==r['attempt_receipt_sha256'];attempt=json.loads(rp.read_text());assert attempt['exit_code']==0
        sources=[rp.parent/'stdout.log',rp.parent/'RESULTS/log.txt']
        for p in sources:assert sha(p)==attempt['artifacts'][str(p.relative_to(rp.parent))];source_bindings[str(p)]=sha(p)
        stdout=sources[0].read_text();log=sources[1].read_text()
        values={name:value for name,value in re.findall(r'^(_\w+)\s+\((?:Float|Int|Str)\)\s+(\S+)',stdout,re.M)}
        keys=['_optimizationLevel','_maxNumOfIterations','_maxNumOfIterationsModel','_epsilonOptimizationModel','_epsilonOptimizationIterationCycle','_performOptimizationsBBL','_isInitGainLossByEmpiricalFreq']
        assert all(k in values for k in keys),(item['id'],values.keys())
        starts=re.findall(r'optimization starting- epsilonOptParam=([^ ]+) epsilonOptIter=\s*([^,]+), MaxNumIterations=(\d+)',log)
        assert starts,item['id']
        rows.append(dict(id=item['id'],variant=item['variant'],**{k.lstrip('_'):values[k] for k in keys},
                         model_optimizer_calls=len(starts),maximum_logged_model_iterations=max(int(x[2]) for x in starts),
                         model_iteration_limit_messages=log.count('Too many iterations in optimizeGainLossModel'),
                         outer_iteration_limit_messages=log.count('Too many iterations in gainLossOptimizer'),
                         source_attempt_sha256=sha(rp)))
    assert len(rows)==306 and len(empty)==6
    args.output.mkdir(parents=True,exist_ok=False);pd.DataFrame(rows).to_csv(args.output/'optimizer_settings.tsv',sep='\t',index=False)
    byvariant={}
    for variant in ['precision-only','precision-cache-refresh']:
        subset=[r for r in rows if r['variant']==variant]
        byvariant[variant]=dict(fits=len(subset),optimization_levels=dict(Counter(r['optimizationLevel'] for r in subset)),
            outer_iteration_limits=dict(Counter(r['maxNumOfIterations'] for r in subset)),model_iteration_limits=dict(Counter(r['maxNumOfIterationsModel'] for r in subset)),
            fits_with_model_iteration_limit_message=sum(r['model_iteration_limit_messages']>0 for r in subset),
            fits_with_outer_iteration_limit_message=sum(r['outer_iteration_limit_messages']>0 for r in subset),
            logged_model_iteration_limits=dict(Counter(str(r['maximum_logged_model_iterations']) for r in subset)))
    source=Path('data/software_audits/fastml-3.11/source/FastML.v3.11/programs/gainLoss/gainLossOptions.cpp')
    text=source.read_text();start=text.index('void gainLossOptions::updateOptimizationLevel(');end=text.index('\n}',start)+2
    (args.output/'optimization_level_source_excerpt.txt').write_text(text[start:end]+'\n')
    result=dict(status='complete_full_native_optimizer_setting_audit',variants=byvariant,empty_inputs=empty,
        pins={str(p):sha(p) for p in [args.plan,source,Path(__file__),root/'receipt.json']},
        source_logs=source_bindings,artifacts={p.name:sha(p) for p in args.output.iterdir()},
        scope='Native effective settings and actual optimizer log diagnostics for all306 fits. A low setting forces one outer/model iteration. Limit messages establish incomplete iteration budgets, not a quantified optimum gap. All outputs remain unqualified; a future refinement must explicitly retain a non-overriding optimization level and verify effective settings.')
    write_json(args.output/'receipt.json',result);print(json.dumps(byvariant,indent=2))


if __name__=='__main__':main()
