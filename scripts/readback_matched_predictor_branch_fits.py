#!/usr/bin/env python3
"""Independently decode every native point-fit tree/report and retain failures."""
import argparse
from collections import Counter
import json
import math
from pathlib import Path
import re

import dendropy

from matched_predictor_branch_fits import load,config,ROLES,SUMMARY_FIELDS
from matched_predictor_branch_inputs import verify
from ancestral_chain_attempt import sha
from run_matched_predictor_branch_fits import STATUS as PRODUCER_STATUS

STATUS='passed_full_matched_predictor_native_point_fit_independent_readback'


def graph(tree,taxa,positions):
    assert sorted(n.taxon.label for n in tree.leaf_node_iter())==sorted(taxa)
    full=sum(2**positions[t] for t in taxa);descendants={};lengths={}
    for node in tree.postorder_node_iter():
        descendants[node]=2**positions[node.taxon.label] if node.is_leaf() else sum(descendants[c] for c in node.child_node_iter())
    for node in tree.preorder_node_iter():
        if node is tree.seed_node:continue
        a=descendants[node];b=full^a
        key=min([a,b],key=lambda x:(x.bit_count(),x))
        assert node.edge.length is not None and math.isfinite(node.edge.length) and node.edge.length>=0
        lengths[key]=lengths.get(key,0.)+node.edge.length
    assert len(lengths)==2*len(taxa)-3
    return lengths,full


def independent_decode(folder,source_config,positions):
    tree=dendropy.Tree.get(path=folder/'fit.treefile',schema='newick',rooting='force-unrooted',preserve_underscores=True)
    lengths,full=graph(tree,source_config['taxa'],positions)
    internal={hex(k) for k in lengths if min(k.bit_count(),(full^k).bit_count())>=2}
    assert internal==set(source_config['internal_splits'])
    report=(folder/'fit.iqtree').read_text();log=(folder/'fit.log').read_text()
    line=next(line for line in report.splitlines() if line.startswith('Input data:'))
    numbers=re.findall(r'\d+',line);assert tuple(map(int,numbers[:2]))==(len(source_config['taxa']),len(source_config['columns']))
    value=float(next(line.split(':',1)[1].split()[0] for line in report.splitlines() if line.startswith('Log-likelihood of the tree:')))
    assert math.isfinite(value)
    serialized=[line for line in report.splitlines() if line.startswith('(') and line.strip().endswith(';')];assert len(serialized)==1
    second=dendropy.Tree.get(data=serialized[0],schema='newick',rooting='force-unrooted',preserve_underscores=True)
    reported,_=graph(second,source_config['taxa'],positions)
    assert set(lengths)==set(reported)
    assert all(math.isclose(lengths[k],reported[k],rel_tol=1e-8,abs_tol=1e-8) for k in lengths)
    notes=[line.strip() for line in (report+'\n'+log).splitlines() if 'WARNING' in line.upper() or 'NOTE:' in line.upper()]
    return dict(log_likelihood_reported=value,branches=[dict(split_mask_hex=hex(k),length=lengths[k],
        internal=min(k.bit_count(),(full^k).bit_count())>=2) for k in sorted(lengths)],diagnostic_messages=notes,
        branch_unit='expected_state_substitutions_per_site',zero_branch_values=sum(x==0 for x in lengths.values()),scientific_eligibility=False)


def check_point(actual,expected):
    assert set(actual)==set(expected)
    for k in actual:
        if k=='branches':
            assert len(actual[k])==len(expected[k])
            for a,b in zip(actual[k],expected[k]):
                assert set(a)==set(b) and a['split_mask_hex']==b['split_mask_hex'] and a['internal']==b['internal']
                assert math.isclose(a['length'],b['length'],rel_tol=1e-10,abs_tol=1e-12)
        else:assert actual[k]==expected[k]


def summarize(rows,input_keys,input_cases):
    by_input={k:[] for k in input_keys}
    assert len(rows)==len(input_keys)*7
    assert {(r['input_id'],r['role']) for r in rows}=={(k,r[0]) for k in input_keys for r in ROLES}
    for row in rows:by_input[row['input_id']].append(row)
    intact=sum(all(r['point_estimate'] is not None for r in values) for values in by_input.values())
    return dict(input_comparison_cases=input_cases,unique_ready_inputs=len(input_keys),native_roles=len(rows),
        native_status_counts=dict(Counter(r['status'] for r in rows)),intact_inputs=intact,unresolved_inputs=len(input_keys)-intact,
        point_branch_values=sum(len(r['point_estimate']['branches']) for r in rows if r['point_estimate'] is not None),
        native_seconds_sum=sum(r['elapsed_seconds'] for r in sorted(rows,key=lambda r:(r['input_id'],r['role']))))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    plan=json.loads(a.plan.read_text());source,bindings=load(plan,a.plan);root=Path(plan['output'])
    rp=root/'receipt.json';producer=json.loads(rp.read_text())
    assert producer['status']==PRODUCER_STATUS and producer['plan_sha256']==sha(a.plan)
    verify(producer['source_hashes']);assert all(sha(root/p)==v for p,v in producer['artifacts'].items())
    rows=[];review=[]
    for key,sc in sorted(source['configs'].items()):
        for role in ROLES:
            row=json.loads((root/'roles'/(key+'-'+role[0]+'.json')).read_text())
            assert (row['input_id'],row['role'])==(key,role[0])
            folder=root/'native'/key/role[0];expected=config(plan,source['root'],key,sc,role)
            assert json.loads((folder/'configuration.json').read_text())==expected
            assert [p.name for p in folder.glob('attempt-*') if p.is_dir()]==['attempt-0001']
            attempt=folder/'attempt-0001';native=json.loads((attempt/'receipt.json').read_text())
            assert row['native_receipt']==str((attempt/'receipt.json').resolve()) and row['native_receipt_sha256']==sha(attempt/'receipt.json')
            assert row['exit_code']==native['exit_code'] and row['native_status']==native['status']
            assert row['elapsed_seconds']==native['elapsed_seconds'] and row['seed']==expected['seed']
            assert json.loads((attempt/'command.json').read_text())==[x.replace('{attempt}',str(attempt.resolve())) for x in expected['command']]
            assert all(sha(attempt/p)==v for p,v in native['artifacts'].items())
            if native['exit_code']!=0:
                assert row['status']=='native_unsuccessful_retained' and row['point_estimate'] is None
            elif row['point_estimate'] is not None:
                check_point(row['point_estimate'],independent_decode(attempt,sc,source['axes']['positions']))
                assert row['status']=='native_point_output_integrity_checked_not_model_qualified'
            else:
                assert row['status']=='native_output_requires_review';review.append(dict(input_id=key,role=role[0]))
            assert row['scientific_eligibility'] is False;rows.append(row)
        print('independent_matched_predictor_input_fits_verified',key,'roles',len(rows),flush=True)
    summary=summarize(rows,source['configs'],source['completion']['comparison_cases'])
    assert all(summary[k]==producer[k] for k in SUMMARY_FIELDS)
    assert len(list((root/'roles').glob('*.json')))==len(rows)
    verify(bindings);verify(producer['source_hashes']);assert all(sha(root/p)==v for p,v in producer['artifacts'].items())
    bindings[str(rp)]=sha(rp)
    result=dict(status=STATUS,plan_sha256=sha(a.plan),producer_receipt_sha256=sha(rp),**summary,
        native_output_reviews_retained=review,source_hashes=bindings,scientific_eligibility=False,scope=plan['scope'])
    with (root/'readback.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
