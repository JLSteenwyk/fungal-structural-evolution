#!/usr/bin/env python3
"""Independently read every emitted matched-site record against saved arrays."""
import csv,gzip,hashlib,json,math
from collections import defaultdict
from pathlib import Path
import numpy as np


def main():
    root=Path('results/ancestral/refined-whole-domain-probability-comparison-20260927-v1')
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    receipt=json.loads((root/'receipt.json').read_text())
    for p,h in receipt['artifacts'].items():assert sha(root/p)==h
    folders={'whole':Path('results/ancestral/refined-whole-ancestors-20260927-v2'),'domain':Path('results/ancestral/refined-domain-ancestors-20260927-v1')}
    arrays={};aa='ARNDCQEGHILKMFPSTWYV'
    def vector(label,job,level,col):
        key=label,job
        if key not in arrays:
            with np.load(folders[label]/(job+'.npz')) as z:
                assert ''.join(z['amino_acids'].tolist())==aa and z['levels'].tolist()==[0,1,2]
                arrays[key]=z['posterior'].copy()
        return arrays[key][level,col-1].tolist()
    valid={}
    corr=Path('results/ancestral/whole-domain-alignment-correspondence-20260927-v1/column_correspondence.tsv.gz')
    with gzip.open(corr,'rt') as h:
        for r in csv.DictReader(h,delimiter='\t'):
            if r['status']=='matched_projected_column':valid[tuple(r[k] for k in ['family','boundary','whole_method','domain_method','whole_column','domain_column'])]=r['coordinate_signature_sha256']
    totals=defaultdict(lambda:[0,0,0,0.]);seen=set();counts=[0,0,0];maximum=0.
    with gzip.open(root/'matched_site_comparisons.tsv.gz','rt') as h:
        for r in csv.DictReader(h,delimiter='\t'):
            key=tuple(r[k] for k in ['family','boundary','whole_method','domain_method','whole_column','domain_column'])
            assert valid[key]==r['coordinate_signature_sha256']
            level=int(r['level']);pair=(r['whole_job'],r['domain_job'],level)
            identity=pair+(r['whole_column'],r['domain_column']);assert identity not in seen;seen.add(identity)
            x=vector('whole',r['whole_job'],level,int(r['whole_column']));y=vector('domain',r['domain_job'],level,int(r['domain_column']))
            ix=max(range(20),key=lambda i:x[i]);iy=max(range(20),key=lambda i:y[i]);tv=math.fsum(abs(a-b) for a,b in zip(x,y))/2
            assert r['whole_map']==aa[ix] and r['domain_map']==aa[iy]
            assert abs(float(r['whole_map_probability'])-x[ix])<1e-12 and abs(float(r['domain_map_probability'])-y[iy])<1e-12
            assert abs(float(r['total_variation'])-tv)<1e-12
            change=ix!=iy;high=change and x[ix]>=.9 and y[iy]>=.9
            assert (r['opposing_at_least_090']=='True')==high
            counts[0]+=1;counts[1]+=change;counts[2]+=high;maximum=max(maximum,tv)
            t=totals[pair];t[0]+=1;t[1]+=change;t[2]+=high;t[3]+=tv
    with (root/'node_summary.tsv').open() as h:
        summaries=list(csv.DictReader(h,delimiter='\t'))
    assert len(summaries)==1872
    for r in summaries:
        t=totals[r['whole_job'],r['domain_job'],int(r['level'])]
        assert t[:3]==[int(r[k]) for k in ['matched_sites','map_disagreements','opposing_at_least_090']]
        assert math.isclose(t[3],float(r['total_variation_sum']),abs_tol=1e-9)
    assert counts==[receipt[k] for k in ['matched_node_sites','map_disagreements','opposing_at_least_090']]
    assert abs(maximum-receipt['maximum_total_variation'])<1e-12
    output=dict(status='passed_all_saved_site_probability_readbacks',site_records=counts[0],map_disagreements=counts[1],opposing_at_least_090=counts[2],node_summaries=len(summaries),source_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(__file__),scope='Full serialized-site and node-summary readback against saved probability arrays and exact coordinate correspondence; this does not independently refit models or establish biological validity.')
    Path('metadata/refined_whole_domain_probability_readback_completed_20260927.json').write_text(json.dumps(output,indent=2)+'\n');print(json.dumps(output))


if __name__=='__main__':main()
