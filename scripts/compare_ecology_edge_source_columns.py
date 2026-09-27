"""Count exact original alignment columns shared by ecological pairs across predictors."""
import csv
import json
from pathlib import Path
from Bio import SeqIO
import pandas as pd
from screen_duplication_domain_alignment_coverage import sha

out=Path('results/ecology/required-edge-source-columns-20260927-v1')
edge_root=Path('results/ecology/required-edge-structure-coverage-20260927-v1')
er=json.loads((edge_root/'receipt.json').read_text())
assert er['status']=='complete_ecology_required_edge_structure_coverage'
pins={str(edge_root/'receipt.json'):sha(edge_root/'receipt.json')}
for name,h in er['artifacts'].items(): assert sha(edge_root/name)==h;pins[str(edge_root/name)]=h
focal=er['focal_taxa'];data={};taxa=None
for source in ['afdb','esmfold']:
    pp=Path(f'metadata/qualified_{source}_reviewed45_ecology_plan_20260927.json')
    plan=json.loads(pp.read_text());pins[str(pp)]=sha(pp)
    root=Path(plan['inputs']);rp=root/'receipt.json';r=json.loads(rp.read_text());ap=Path(plan['input_readback']);a=json.loads(ap.read_text())
    assert a['status']=='passed_complete_paired_inputs_from_qualified_arrays_readback' and a['source_receipt_sha256']==sha(rp)
    pins[str(rp)]=sha(rp);pins[str(ap)]=sha(ap)
    tp=Path(plan['output'])/'taxa.tsv';current=set(pd.read_csv(tp,sep='\t').taxon_id)
    assert sha(tp)==er['pins'][str(tp)]
    if taxa is None:taxa=current
    else:assert taxa==current
    masks={}
    for name,h in r['artifacts'].items():assert sha(root/name)==h;pins[str(root/name)]=h
    for name in sorted(r['artifacts']):
        if not name.endswith('/aa.faa'):continue
        marker=name.split('/')[0]
        columns=list(csv.DictReader((root/marker/'columns.tsv').open(),delimiter='\t'))
        assert [int(x['paired_column_1based']) for x in columns]==list(range(1,len(columns)+1))
        original=[int(x['alignment_column_1based']) for x in columns]
        assert len(original)==len(set(original)) and all(x['marker']==marker for x in columns)
        aa={x.id:str(x.seq) for x in SeqIO.parse(root/name,'fasta')}
        di={x.id:str(x.seq) for x in SeqIO.parse(root/marker/'3di.faa','fasta')}
        assert set(aa)==set(di)
        masks[marker]={}
        for t in sorted(taxa&set(aa)):
            assert len(aa[t])==len(di[t])==len(original)
            assert [i for i,c in enumerate(aa[t]) if c!='?']==[i for i,c in enumerate(di[t]) if c!='?']
            masks[marker][t]={col:c for col,c in zip(original,aa[t]) if c!='?'}
    data[source]=masks
markers=sorted(set(data['afdb'])|set(data['esmfold']));rows=[]
for f in sorted(focal):
    for other in sorted(taxa-{f}):
        for marker in markers:
            record=dict(focal_taxon=f,comparison_taxon=other,marker=marker)
            sets={};maps={}
            for source in ['afdb','esmfold']:
                m=data[source].get(marker,{})
                record[source+'_both_taxa_eligible']=f in m and other in m
                left,right=m.get(f,{}),m.get(other,{})
                maps[source]=(left,right);sets[source]=set(left)&set(right)
                record[source+'_shared_columns']=len(sets[source])
            common=sets['afdb']&sets['esmfold']
            matching={c for c in common if all(maps['afdb'][i][c]==maps['esmfold'][i][c] for i in [0,1])}
            record.update(four_way_observed_columns=len(common),four_way_matching_aa_columns=len(matching),
                predictor_aa_mismatch_columns=len(common-matching),at_least_50_matching_columns=len(matching)>=50)
            rows.append(record)
frame=pd.DataFrame(rows)
assert len(frame)==len(focal)*(len(taxa)-1)*len(markers)
assert not frame.duplicated(['focal_taxon','comparison_taxon','marker']).any()
# All per-source counts must recover the previously independently checked ledger.
old=pd.read_csv(edge_root/'pair_marker_coverage.tsv',sep='\t')
for source in ['afdb','esmfold']:
    keyed=old[old.source.eq(source)].set_index(['focal_taxon','comparison_taxon','marker'])
    selected=frame[frame[source+'_both_taxa_eligible']].set_index(['focal_taxon','comparison_taxon','marker'])
    assert set(keyed.index)==set(selected.index)
    assert keyed.shared_qualified_columns.sort_index().tolist()==selected[source+'_shared_columns'].sort_index().tolist()
summary=[]
for f,part in frame.groupby('focal_taxon'):
    qualified=part[part.at_least_50_matching_columns]
    summary.append(dict(focal_taxon=f,pair_marker_settings=len(part),comparison_taxa=part.comparison_taxon.nunique(),
        jointly_qualifying_comparison_taxa=qualified.comparison_taxon.nunique(),jointly_qualifying_markers=qualified.marker.nunique(),
        jointly_qualifying_pair_markers=len(qualified),four_way_matching_pair_columns=int(part.four_way_matching_aa_columns.sum()),
        predictor_aa_mismatch_pair_columns=int(part.predictor_aa_mismatch_columns.sum())))
out.mkdir(parents=True,exist_ok=False)
frame.to_csv(out/'pair_marker_columns.tsv',sep='\t',index=False)
pd.DataFrame(summary).to_csv(out/'summary.tsv',sep='\t',index=False)
pd.testing.assert_frame_equal(pd.read_csv(out/'pair_marker_columns.tsv',sep='\t'),frame)
for p,h in pins.items():assert sha(p)==h,p
result=dict(status='complete_ecological_edge_cross_predictor_columns',focal_taxa=len(focal),reviewed_taxa=len(taxa),union_markers=len(markers),pair_marker_settings=len(frame),pins=pins,script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in out.iterdir()},scope='All three focal taxa, all 44 reviewed partners and union marker grid retained, including absent observations. Original alignment positions and amino acids must match within each taxon across sources. Same marker or source coverage alone is insufficient. Pair-column counts reuse observations and are not independent residues or ecological transitions. No effect estimates.')
(out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
print(pd.DataFrame(summary).to_string(index=False))
