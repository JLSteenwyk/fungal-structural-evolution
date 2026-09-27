#!/usr/bin/env python3
"""Snapshot all available independently verified ancestor exclusion constraints."""
import csv,gzip,json
from collections import defaultdict
from pathlib import Path
import numpy as np
from prepare_case_ancestral_neighborhoods import sha


def main():
    root=Path('results/ancestral/indel-compatibility-constraints-20260927-v1');rp=root/'receipt.json';receipt=json.loads(rp.read_text());cp=root/'mutually_exclusive_gap_sets.jsonl.gz';assert sha(cp)==receipt['artifacts'][cp.name]
    sets=defaultdict(list)
    with gzip.open(cp,'rt') as h:
        for line in h:
            row=json.loads(line);sets[row['encoding']].append(row)
    mp=Path('results/ancestral/indel-candidate-node-mapping-20260927-v1/candidate_nodes.tsv')
    mappings=defaultdict(list)
    with mp.open() as h:
        for row in csv.DictReader(h,delimiter='\t'):
            if row['guide']=='profile' and row['status']=='matched_unrooted_vertex':mappings[row['job_id']].append(row)
    rows=[];pins={str(rp):sha(rp),str(cp):sha(cp),str(mp):sha(mp)}
    for ap in sorted(Path('results/ancestral/stable-indel-posterior-audit-20260927-v1').glob('*.json')):
        if ap.name=='receipt.json':continue
        audit=json.loads(ap.read_text());jid=audit['job_id'];pins[str(ap)]=sha(ap)
        if audit['status']=='empty_input_verified':continue
        assert audit['status']=='all_node_gap_probabilities_independently_verified'
        folder=Path('results/ancestral/stable-indel-ancestors-20260927-v1')/jid;pr=folder/'receipt.json';assert sha(pr)==audit['source_receipt_sha256'];r=json.loads(pr.read_text());assert sha(folder/'posterior.npz')==r['artifacts']['posterior.npz'];pins[str(pr)]=sha(pr)
        encoding=jid.removesuffix('-all_taxa').removesuffix('-observed_mask')
        with np.load(folder/'posterior.npz',allow_pickle=False) as data:
            names=list(data['node_names'])
            for target in mappings[jid]:
                p=data['probability_gap'][names.index(target['posterior_node'])];constraints=sets[encoding]
                sums=[float(p[np.array(s['characters'])-1].sum()) for s in constraints]
                largest=int(np.argmax(sums)) if sums else None
                rows.append(dict(job_id=jid,level=int(target['level']),constraint_sets=len(sums),violating_sets=sum(x>1+1e-9 for x in sums),maximum_probability_sum=max(sums,default=0),witness_characters=constraints[largest]['characters'] if largest is not None else []))
    result=dict(status='available_verified_posterior_compatibility_snapshot',candidate_node_models=len(rows),candidate_nodes_with_violations=sum(r['violating_sets']>0 for r in rows),maximum_probability_sum=max((r['maximum_probability_sum'] for r in rows),default=0),rows=rows,pins=pins,script_sha256=sha(__file__),scope='Only available independently verified production outputs. Sums>1+1e-9 violate a necessary exclusion constraint for any joint distribution over exact gap runs. Redundant sets and model replicates are not independent biological events. No probabilities projected or sequences assembled; complete snapshot/final assembly remains pending.')
    Path('metadata/indel_posterior_compatibility_snapshot_20260927.json').write_text(json.dumps(result,indent=2)+'\n');print({k:v for k,v in result.items() if k not in ['rows','pins']})

if __name__=='__main__':main()
