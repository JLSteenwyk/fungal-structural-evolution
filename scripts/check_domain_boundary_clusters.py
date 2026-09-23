#!/usr/bin/env python3
"""Exercise boundary cluster categories and malformed paired identities."""
import json
from pathlib import Path
from compare_domain_boundary_clusters import compare
from catalog_whole_proteome_structures import sha


def main():
    pairs=[('a','a','identical_interval'),('a','b','distinct_intervals_same_cluster'),('a','c','distinct_intervals_different_clusters'),('x','a','alignment_unclustered'),('a','x','envelope_unclustered'),('x','y','both_unclustered')]
    links=[dict(model_key='m',hit_id=str(i),boundary=boundary,interval_id=interval) for i,(a,e,_) in enumerate(pairs) for boundary,interval in [('alignment',a),('envelope',e)]]
    rows=compare(links,{'a':'a','b':'a','c':'c'})
    assert [r['disposition'] for r in rows]==[p[2] for p in pairs]
    rejected=0
    for invalid in [links+links[:1],links[1:],[dict(links[0],boundary='unknown')]+links[1:]]:
        try:compare(invalid,{'a':'a','b':'a','c':'c'})
        except ValueError:rejected+=1
    assert rejected==3
    result={'status':'passed_all_six_boundary_dispositions_and_identity_rejections','categories':6,'rejected_cases':rejected,'script_sha256':sha(__file__),'controller_sha256':sha('scripts/compare_domain_boundary_clusters.py'),'scope':'Synthetic pair classification and identity fixtures. Production clustering and full-data output remain unvalidated.'}
    Path('metadata/domain_boundary_cluster_fixture_checks.json').write_text(json.dumps(result,indent=2)+'\n');print(result)


if __name__=='__main__':main()
