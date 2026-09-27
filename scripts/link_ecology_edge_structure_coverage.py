"""Join every ever-required ecological split to source-specific structural coverage."""
from collections import Counter
import csv
import json
from pathlib import Path
import pandas as pd
from screen_duplication_domain_alignment_coverage import sha

out=Path('results/ecology/required-edge-structure-coverage-20260927-v1')
root=Path('results/ecology/bootstrap-uncertainty-summary-20260927-v1')
auditpath=Path('metadata/ecology_bootstrap_uncertainty_summary_readback_20260927.json')
audit=json.loads(auditpath.read_text());rp=root/'receipt.json';r=json.loads(rp.read_text())
assert audit['status']=='passed_full_ecology_bootstrap_summary_readback' and sha(rp)==audit['source_receipt_sha256']
pins={str(rp):sha(rp),str(auditpath):sha(auditpath)}
for name,h in r['artifacts'].items(): assert sha(root/name)==h;pins[str(root/name)]=h
splits=pd.read_csv(root/'split_uncertainty.tsv',sep='\t')
ids=set(splits.loc[splits.required_change.gt(0),'split_id'])
edges=splits[splits.split_id.isin(ids)].copy()
assert len(ids)==3 and len(edges)==12 and edges.canonical_side.nunique()==3
focal=set(edges.canonical_side);assert all(',' not in x and ';' not in x for x in focal)
summary=[];details=[]
for source in ['afdb','esmfold']:
    folder=Path(f'results/ecology/qualified-{source}-reviewed45-overlap-20260927-v1')
    receiptpath=folder/'receipt.json';receipt=json.loads(receiptpath.read_text())
    ap=Path(f'metadata/qualified_{source}_reviewed45_ecology_readback_20260927.json');a=json.loads(ap.read_text())
    assert a['status']=='passed_full_qualified_ecology_overlap_matrix_readback' and a['producer_receipt_sha256']==sha(receiptpath)
    pins[str(ap)]=sha(ap);pins[str(receiptpath)]=sha(receiptpath)
    for name,h in receipt['artifacts'].items(): assert sha(folder/name)==h;pins[str(folder/name)]=h
    taxa=pd.read_csv(folder/'taxa.tsv',sep='\t').set_index('taxon_id')
    pairs=pd.read_csv(folder/'pairs.tsv',sep='\t');markers=pd.read_csv(folder/'pair_marker_columns.tsv',sep='\t')
    for taxon in sorted(focal):
        assert taxon in taxa.index
        paired=pairs[pairs.taxon_a.eq(taxon)|pairs.taxon_b.eq(taxon)]
        assert len(paired)==44
        selected=markers[markers.taxon_a.eq(taxon)|markers.taxon_b.eq(taxon)]
        check=Counter()
        for row in selected.itertuples(index=False):
            other=row.taxon_b if row.taxon_a==taxon else row.taxon_a
            qualifies=row.shared_qualified_columns>=50
            check[other]+=int(qualifies)
            details.append(dict(source=source,focal_taxon=taxon,comparison_taxon=other,
                comparison_published_state=taxa.loc[other,'state'],marker=row.marker,
                shared_qualified_columns=int(row.shared_qualified_columns),at_least_50_columns=qualifies))
        for row in paired.itertuples(index=False):
            other=row.taxon_b if row.taxon_a==taxon else row.taxon_a
            assert check[other]==row.markers_with_at_least_50_shared_columns
        for _,edge in edges[edges.canonical_side.eq(taxon)].iterrows():
            summary.append(dict(source=source,**edge.to_dict(),species_name=taxa.loc[taxon,'species_name'],
                published_state=taxa.loc[taxon,'state'],eligible_markers=int(taxa.loc[taxon,'eligible_markers']),
                qualified_observations=int(taxa.loc[taxon,'qualified_observations']),
                reviewed_comparison_taxa=44,comparison_taxa_with_50_column_marker=sum(v>0 for v in check.values()),
                qualifying_pair_marker_links=sum(check.values()),
                distinct_markers_in_qualifying_pairs=selected.loc[selected.shared_qualified_columns.ge(50),'marker'].nunique()))
out.mkdir(parents=True,exist_ok=False)
for name,rows in [('edge_coverage.tsv',summary),('pair_marker_coverage.tsv',details)]:
    df=pd.DataFrame(rows);df.to_csv(out/name,sep='\t',index=False)
    pd.testing.assert_frame_equal(pd.read_csv(out/name,sep='\t'),df,check_dtype=False)
# Independent scalar readback of every per-source/focal denominator.
ledger=list(csv.DictReader((out/'pair_marker_coverage.tsv').open(),delimiter='\t'))
for row in csv.DictReader((out/'edge_coverage.tsv').open(),delimiter='\t'):
    selected=[x for x in ledger if x['source']==row['source'] and x['focal_taxon']==row['canonical_side'] and int(x['shared_qualified_columns'])>=50]
    assert len(selected)==int(row['qualifying_pair_marker_links'])
    assert len({x['comparison_taxon'] for x in selected})==int(row['comparison_taxa_with_50_column_marker'])
    assert len({x['marker'] for x in selected})==int(row['distinct_markers_in_qualifying_pairs'])
for p,h in pins.items(): assert sha(p)==h,p
result=dict(status='complete_ecology_required_edge_structure_coverage',edge_source_conditions=len(summary),pair_marker_rows=len(details),focal_taxa=sorted(focal),pins=pins,script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in out.iterdir()},scope='All conditions for every split required in any bootstrap reconstruction retained, including zero-required conditions. Source-specific shared-mask coverage only. Candidate partners include all 44 reviewed taxa; labels are not recoded to binary states. No independent origins, valid ecological contrasts, branch-effect estimates or power established.')
(out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
print(pd.DataFrame(summary)[['source','canonical_side','eligible_markers','comparison_taxa_with_50_column_marker','qualifying_pair_marker_links']].drop_duplicates().to_string(index=False))
