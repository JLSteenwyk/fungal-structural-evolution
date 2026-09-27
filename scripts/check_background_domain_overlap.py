#!/usr/bin/env python3
"""Check exact reusable interval identities; overlap does not certify native results."""
import argparse,csv,json
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--background',type=Path,required=True);p.add_argument('--existing',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
def load(root):
    r=json.loads((root/'receipt.json').read_text());bindings={str(root/'receipt.json'):sha(root/'receipt.json')}
    for n in ['intervals.jsonl','domain_pairs.tsv']:
        assert sha(root/n)==r['artifacts'][n];bindings[str(root/n)]=sha(root/n)
    intervals={}
    for line in (root/'intervals.jsonl').read_text().splitlines():
        row=json.loads(line);assert row['interval_id'] not in intervals;intervals[row['interval_id']]=row
    pairs={}
    for row in csv.DictReader((root/'domain_pairs.tsv').open(),delimiter='\t'):
        assert row['domain_pair_key'] not in pairs;pairs[row['domain_pair_key']]=(row['interval_a'],row['interval_b'])
    assert len(intervals)==r['unique_intervals'] and len(pairs)==r['unique_interval_pairs']
    return intervals,pairs,bindings
bi,bp,bb=load(a.background);ei,ep,eb=load(a.existing);shared_intervals=sorted(bi.keys()&ei.keys());shared_pairs=sorted(bp.keys()&ep.keys())
assert all(bi[k]==ei[k] for k in shared_intervals) and all(bp[k]==ep[k] for k in shared_pairs)
assert all(x in shared_intervals for k in shared_pairs for x in bp[k])
for path,h in (bb|eb).items():assert sha(path)==h
result=dict(status='passed_exact_background_domain_inventory_overlap',background_intervals=len(bi),background_pairs=len(bp),shared_intervals=len(shared_intervals),shared_pairs=len(shared_pairs),new_pairs=len(bp)-len(shared_pairs),shared_interval_ids=shared_intervals,shared_pair_keys=shared_pairs,source_bindings=bb|eb,scope='Exact interval descriptors and endpoint identities agree for every overlap. This does not certify existing coordinate serialization or native results; any computation reuse still requires successful independent numeric checks. Full background inventory remains intact.')
a.output.write_text(json.dumps(result,indent=2)+'\n');print({k:v for k,v in result.items() if k not in ['shared_interval_ids','shared_pair_keys','source_bindings','scope']})
