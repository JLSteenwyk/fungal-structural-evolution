#!/usr/bin/env python3
"""Check complete published coverage aggregates against serialized arrays and taxa."""
import argparse
from collections import Counter
from datetime import datetime,timezone
import json
import math
from pathlib import Path

import numpy as np

from audit_selected_taxon_identity_snapshot_v2 import rows,sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--summary',type=Path,required=True)
    p.add_argument('--producer',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();assert not a.output.exists()
    summary=json.loads(a.summary.read_text());producer=json.loads(a.producer.read_text())
    assert summary['status']=='complete_descriptive_full_tree_view_structural_coverage_summary'
    bindings={str(a.summary):sha(a.summary),str(a.producer):sha(a.producer),str(Path(__file__)):sha(__file__)}
    for record in [summary,producer]:
        for field in ['source_hashes','artifacts']:
            for path,digest in record.get(field,{}).items():
                assert sha(path)==digest,path
                assert path not in bindings or bindings[path]==digest,path
                bindings[path]=digest
    roots={Path(path).parent for path in summary['artifacts']};assert len(roots)==1;root=roots.pop()
    producer_roots={Path(path).parent for path in producer['artifacts']};assert len(producer_roots)==1;data=producer_roots.pop()
    axes=json.loads((data/'input_axes.json').read_text());statuses=axes['statuses'];sources=axes['sources']
    view_rows=rows(root/'all_view_source_branch_coverage.tsv');assert len(view_rows)==140
    indexed={(int(r['view_index']),r['source']):r for r in view_rows};assert len(indexed)==140
    totals=Counter()
    for vi,view in enumerate(axes['views']):
        with np.load(data/('view_%02d.npz'%vi),allow_pickle=False) as saved:values=saved['status']
        assert values.dtype==np.dtype('uint8') and values.shape==(2,125,len(view['branches']))
        for si,source in enumerate(sources):
            row=indexed[vi,source];counts=np.bincount(values[si].ravel(),minlength=5);assert len(counts)==5
            assert int(row['branch_marker_cells'])==values[si].size
            for code,status in enumerate(statuses):
                assert int(row[status])==counts[code]
                assert math.isclose(float(row[status+'_fraction']),counts[code]/values[si].size,rel_tol=1e-14,abs_tol=1e-14)
                totals[status]+=int(counts[code])
    assert dict(totals)==producer['status_counts'] and sum(totals.values())==9047500
    taxa=rows(data/'taxon_identity_source_coverage.tsv');assert len(taxa)==526
    lineage=rows(root/'all_manifest_lineage_source_coverage.tsv')
    keys={(r['study_role'],r['manifest_lineage_group'],r['source']) for r in lineage}
    assert len(keys)==len(lineage)==54
    assert keys=={(r['study_role'],r['manifest_lineage_group'],s) for r in taxa for s in sources}
    for row in lineage:
        selected=[r for r in taxa if r['study_role']==row['study_role'] and r['manifest_lineage_group']==row['manifest_lineage_group']]
        counts=[int(r[row['source']+'_usable_markers']) for r in selected]
        assert int(row['panel_entries'])==len(selected) and int(row['usable_marker_taxon_cells'])==sum(counts)
        for threshold in [1,10,50,100]:assert int(row['entries_at_least_'+str(threshold)+'_markers'])==sum(v>=threshold for v in counts)
        assert int(row['entries_with_identity_review_flags'])==sum(bool(r['identity_review_flags']) for r in selected)
    for source,expected in [('AlphaFold',31134),('ESMFold',23799)]:
        assert sum(int(r['usable_marker_taxon_cells']) for r in lineage if r['source']==source)==expected
        assert sum(int(r['panel_entries']) for r in lineage if r['source']==source)==526
    for path,digest in bindings.items():assert sha(path)==digest,path
    result=dict(status='passed_full_structural_marker_coverage_summary_readback',checked_utc=datetime.now(timezone.utc).isoformat(),
        view_source_histograms_checked=140,serialized_branch_cells_checked=9047500,lineage_source_rows_checked=54,
        full_taxon_rows_checked=526,summary=str(a.summary),summary_sha256=sha(a.summary),source_hashes=bindings,
        figure_visual_review_separate=True,scientific_eligibility=False,
        scope='All140serializedsummaryhistograms/fractions checked directly against complete NPZ status arrays, independent of case-table aggregation. All54lineage-source rows checked against all526taxa and exact qualified source totals. Full previous graph pruning remains separately required; no statistical or biological acceptance.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','scope']},indent=2))


if __name__=='__main__':main()
