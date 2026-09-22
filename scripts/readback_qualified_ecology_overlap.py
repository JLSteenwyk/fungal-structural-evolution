#!/usr/bin/env python3
"""Recompute ecological overlap with dense mask matrix multiplication."""
import argparse
import csv
import itertools
import json
from pathlib import Path
import numpy as np
from readback_whole_proteome_family_coverage import sha


def rows(path):
    with Path(path).open() as handle:return list(csv.DictReader(handle,delimiter='\t'))


def fasta(path):
    records={};key=None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            key=line[1:].split()[0]
            if key in records:raise ValueError('Duplicate taxon')
            records[key]=''
        elif line.strip():records[key]+=line.strip()
    return records


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());out=Path(plan['output'])
    if args.output.exists():raise FileExistsError(args.output)
    for path,digest in plan['pins'].items():assert sha(path)==digest
    receipt=json.loads((out/'receipt.json').read_text());assert receipt['plan_sha256']==sha(args.plan)
    for name,digest in receipt['artifacts'].items():assert sha(out/name)==digest
    evidence=rows(plan['evidence']);taxa=[r['taxon_id'] for r in evidence];n=len(taxa);index={t:i for i,t in enumerate(taxa)}
    source=json.loads((Path(plan['inputs'])/'receipt.json').read_text())
    masks={};present={};dots={}
    for name,digest in source['artifacts'].items():
        assert sha(Path(plan['inputs'])/name)==digest
        if not name.endswith('/aa.faa'):continue
        marker=name.split('/')[0];aa=fasta(Path(plan['inputs'])/name);di=fasta(Path(plan['inputs'])/marker/'3di.faa')
        assert set(aa)==set(di)
        length=len(next(iter(aa.values())))
        matrix=np.zeros((n,length),dtype=np.int32);included=np.zeros(n,dtype=np.int32)
        for taxon,i in index.items():
            if taxon not in aa:continue
            assert len(aa[taxon])==len(di[taxon])==length
            left=np.array(list(aa[taxon]))!='?';right=np.array(list(di[taxon]))!='?'
            assert np.array_equal(left,right)
            matrix[i]=left;included[i]=1
        masks[marker]=matrix;present[marker]=included;dots[marker]=matrix@matrix.T
    assert len(masks)==receipt['markers']
    detail={}
    for row in rows(out/'pair_marker_columns.tsv'):
        key=(row['taxon_a'],row['taxon_b'],row['marker']);assert key not in detail
        detail[key]=int(row['shared_qualified_columns'])
    checked_pairs=0
    pair_rows=rows(out/'pairs.tsv');assert len(pair_rows)==n*(n-1)//2
    for row,(i,j) in zip(pair_rows,itertools.combinations(range(n),2)):
        a,b=evidence[i],evidence[j]
        for suffix,e in [('a',a),('b',b)]:
            assert row['taxon_'+suffix]==e['taxon_id'] and row['species_'+suffix]==e['species_name'] and row['state_'+suffix]==e['state']
        assert row['different_published_states']==str(a['state']!=b['state'])
        expected={m:int(dots[m][i,j]) for m in masks if present[m][i] and present[m][j]}
        for marker,count in expected.items():assert detail.pop((taxa[i],taxa[j],marker))==count
        selected=sorted(m for m,v in expected.items() if v>=plan['minimum_shared_columns'])
        assert int(row['both_eligible_markers'])==len(expected)
        assert int(row['markers_with_shared_columns'])==sum(v>0 for v in expected.values())
        assert int(row['markers_with_at_least_50_shared_columns'])==len(selected)
        assert row['marker_ids_at_least_50']==';'.join(selected)
        assert int(row['shared_qualified_columns'])==sum(expected.values());checked_pairs+=1
    assert not detail
    taxon_rows=rows(out/'taxa.tsv');assert len(taxon_rows)==n
    for row,e in zip(taxon_rows,evidence):
        i=index[e['taxon_id']]
        for key in ('taxon_id','species_name','state'):assert row[key]==e[key]
        assert int(row['eligible_markers'])==sum(p[i] for p in present.values())
        assert int(row['qualified_observations'])==sum(int(m[i].sum()) for m in masks.values())
    groups=rows(out/'groups.tsv');assert {r['provisional_group'] for r in groups}=={r['provisional_transition_group'] for r in evidence}
    for row in groups:
        members=[i for i,e in enumerate(evidence) if e['provisional_transition_group']==row['provisional_group']]
        assert int(row['taxa'])==len(members) and row['taxon_ids']==';'.join(taxa[i] for i in members)
        counts=[int(np.all(masks[m][members],axis=0).sum()) for m in masks if present[m][members].all()]
        assert int(row['all_taxa_eligible_markers'])==len(counts)
        assert int(row['markers_with_common_columns'])==sum(v>0 for v in counts)
        assert int(row['markers_with_at_least_50_common_columns'])==sum(v>=plan['minimum_shared_columns'] for v in counts)
        assert int(row['common_qualified_columns'])==sum(counts)
    result=dict(status='passed_full_qualified_ecology_overlap_matrix_readback',taxa=n,pairs=checked_pairs,groups=len(groups),markers=len(masks),producer_receipt_sha256=sha(out/'receipt.json'),script_sha256=sha(__file__),scope='Independent FASTA parsing and dense integer mask matrix multiplication checks every taxon, pair, pair-marker and provisional group count. Not an ecological effect or independent-transition test.')
    with args.output.open('x') as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
