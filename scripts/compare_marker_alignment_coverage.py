#!/usr/bin/env python3
"""Inventory all marker/taxon coverage changes across the two species matrices."""
import csv
import hashlib
import json
from collections import Counter,defaultdict
from pathlib import Path
from Bio import SeqIO


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    sources={};pins={};taxa=None;markers=None
    for label in ['profile','mafft']:
        root=Path('results/phylogeny')/(label+'-matrix-50-v1')
        receipt=json.loads((root/'receipt.json').read_text());pins[str(root/'receipt.json')]=sha(root/'receipt.json')
        for n,d in receipt['artifacts'].items():assert sha(root/n)==d
        matrix={r.id:str(r.seq) for r in SeqIO.parse(root/'matrix.faa','fasta')};columns=defaultdict(list)
        for r in read(root/'site_mapping.tsv'):columns[r['marker']].append(int(r['matrix_column_1based'])-1)
        assert len(matrix)==526 and len(columns)==125
        if taxa is None:taxa=set(matrix);markers=set(columns)
        assert taxa==set(matrix) and markers==set(columns)
        sources[label]=(matrix,columns)
    profileproof=Path('results/phylogeny/marker-tree-support-complete-readback-v1/receipt.json')
    pr=json.loads(profileproof.read_text());assert pr['status']=='passed_all_snapshot_input_and_graph_split_readbacks' and pr['markers']==125
    pins[str(profileproof)]=sha(profileproof)
    mafft_plan={r['marker']:r for r in read('metadata/mafft_marker_tree_input_plan_20260927.tsv')}
    rows=[];summary=[]
    for marker in sorted(markers):
        values={};sets={}
        for label,(matrix,columns) in sources.items():
            cols=columns[marker];assert cols==list(range(cols[0],cols[-1]+1))
            minimum=max(50,(3*len(cols)+9)//10)
            observed={t:sum(c in 'ACDEFGHIKLMNPQRSTVWY' for c in seq[cols[0]:cols[-1]+1]) for t,seq in matrix.items()}
            values[label]=(len(cols),minimum,observed);sets[label]={t for t,n in observed.items() if n>=minimum}
        actual={r.id for r in SeqIO.parse(Path('results/phylogeny/marker-gene-trees-v2')/marker/'input.faa','fasta')}
        assert actual==sets['profile']
        assert len(sets['mafft'])==int(mafft_plan[marker]['taxa'])
        counter=Counter()
        for taxon in sorted(taxa):
            a,b=taxon in sets['profile'],taxon in sets['mafft']
            disposition='both' if a and b else 'profile_only' if a else 'mafft_only' if b else 'neither'
            counter[disposition]+=1
            rows.append(dict(marker=marker,taxon_id=taxon,profile_columns=values['profile'][0],mafft_columns=values['mafft'][0],profile_minimum=values['profile'][1],mafft_minimum=values['mafft'][1],profile_observed=values['profile'][2][taxon],mafft_observed=values['mafft'][2][taxon],disposition=disposition))
        common=sets['profile']&sets['mafft']
        assert len(common)>=4
        summary.append(dict(marker=marker,profile_columns=values['profile'][0],mafft_columns=values['mafft'][0],profile_taxa=len(sets['profile']),mafft_taxa=len(sets['mafft']),**{k:counter[k] for k in ['both','profile_only','mafft_only','neither']},common_taxa=';'.join(sorted(common)),profile_only_taxa=';'.join(sorted(sets['profile']-common)),mafft_only_taxa=';'.join(sorted(sets['mafft']-common))))
    root=Path('results/phylogeny/marker-alignment-coverage-20260927-v1');root.mkdir(parents=True,exist_ok=False)
    for name,data in [('marker_taxa.tsv',rows),('marker_summary.tsv',summary)]:
        with (root/name).open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');w.writeheader();w.writerows(data)
    result=dict(status='complete_full_marker_alignment_coverage_comparison_pending_readback',markers=125,taxa=526,marker_taxon_rows=len(rows),dispositions=dict(Counter(r['disposition'] for r in rows)),markers_with_different_taxon_sets=sum(r['profile_only']+r['mafft_only']>0 for r in summary),common_taxa_min=min(r['both'] for r in summary),common_taxa_max=max(r['both'] for r in summary),sources=pins,script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in root.iterdir()},scope='All matrix-derived coverage cells and exact common-taxon sets. Profile sets checked against completed inputs; MAFFT sets checked against full planning grid, not unfinished inferred trees. Eligibility uses the same absolute/relative rule but column sets differ. Future tree comparisons must verify emitted MAFFT sets; pruning to common taxa is not refitting.')
    (root/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
