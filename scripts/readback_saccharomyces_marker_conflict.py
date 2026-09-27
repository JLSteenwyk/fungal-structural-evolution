#!/usr/bin/env python3
"""Reconstruct marker conflict rows directly from DendroPy gene trees."""
import csv
import hashlib
import json
from collections import defaultdict, Counter
from pathlib import Path
import dendropy


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def canon(side,taxa):return min(tuple(sorted(side)),tuple(sorted(taxa-side)),key=lambda s:(len(s),s))


def main():
    root=Path('results/phylogeny/saccharomyces-marker-conflict-20260927-v1')
    receipt=json.loads((root/'receipt.json').read_text())
    for n,d in receipt['artifacts'].items():assert sha(root/n)==d
    source=Path('results/phylogeny/marker-tree-support-complete-v1/receipt.json')
    assert sha(source)==receipt['source_support_receipt_sha256']
    inputs=json.loads(source.read_text())['inputs']
    conflict=read('results/phylogeny/pmsf-three-run-topology-sensitivity-20260927-v1/conflicts.tsv')
    first=next(r for r in conflict if r['both_meet_existing_support_label']=='True')
    targets={'profile_alignment':set(json.loads(first['split_a_taxa_json'])),'mafft_alignment':set(json.loads(first['split_b_taxa_json']))}
    rows=read(root/'marker_assessments.tsv');grouped=defaultdict(list)
    for row in rows:grouped[row['marker']].append(row)
    assert len(rows)==1000 and set(grouped)=={i['marker'] for i in inputs}
    for item in inputs:
        marker=item['marker'];path=Path('results/phylogeny/marker-gene-trees-v2')/marker/'tree.treefile'
        assert sha(path)==item['tree_sha256']
        tree=dendropy.Tree.get(path=str(path),schema='newick',rooting='force-unrooted',preserve_underscores=True)
        full={n.taxon.label for n in tree.leaf_node_iter()}
        original=[]
        for node in tree.preorder_node_iter():
            if node is tree.seed_node or node.is_leaf():continue
            side={n.taxon.label for n in node.leaf_iter()}
            original.append((side,float(node.label) if node.label else None))
        assert len(grouped[marker])==8
        for mode in ['all_available','hybrid_tips_removed']:
            taxa=full if mode=='all_available' else full-{'F27292','F332112'}
            idx={t:1<<i for i,t in enumerate(sorted(taxa))};allbits=(1<<len(taxa))-1
            projected=defaultdict(list)
            for side,value in original:
                key=canon(side&taxa,taxa)
                if len(key)>=2:projected[key].append(value)
            scores={k:min(v) if None not in v else None for k,v in projected.items()}
            bits={k:sum(idx[t] for t in k) for k in scores}
            for row in [r for r in grouped[marker] if r['mode']==mode]:
                key=canon(targets[row['alternative']]&taxa,taxa);x=sum(idx[t] for t in key)
                conflict_scores=[scores[k] for k,y in bits.items() if x&y and x&(allbits^y) and (allbits^x)&y and (allbits^x)&(allbits^y)]
                available=[v for v in conflict_scores if v is not None]
                exact=scores.get(key);maximum=max(available) if available else None;cutoff=int(row['sh_alrt_cutoff'])
                expected='uninformative_taxon_coverage' if len(key)<2 else 'supported_concordance' if exact is not None and exact>=cutoff else 'supported_conflict' if maximum is not None and maximum>=cutoff else 'unresolved'
                assert row['status']==expected
                assert row['restricted_split']==';'.join(key)
                assert row['target_taxa_present']==';'.join(sorted(targets[row['alternative']]&taxa))
                assert int(row['retained_taxa'])==len(taxa) and int(row['removed_hybrid_tips'])==len(full-taxa)
                assert (row['exact_split_present']=='True')==(key in scores)
                assert row['exact_minimum_source_sh_alrt']==('' if exact is None else str(exact))
                assert row['maximum_conflicting_minimum_source_sh_alrt']==('' if maximum is None else str(maximum))
                assert int(row['exact_source_edges'])==len(projected.get(key,[]))
                assert int(row['incompatible_splits'])==len(conflict_scores)
                assert int(row['incompatible_splits_without_complete_support'])==sum(v is None for v in conflict_scores)
    summary=read(root/'summary.tsv');assert len(summary)==8
    for row in summary:
        selected=[r for r in rows if all(r[k]==row[k] for k in ['mode','alternative','sh_alrt_cutoff'])]
        assert len(selected)==int(row['markers'])==125
        counts=Counter(r['status'] for r in selected)
        for k in ['supported_concordance','supported_conflict','unresolved','uninformative_taxon_coverage']:assert int(row[k])==counts[k]
    proof=dict(status='passed_full_saccharomyces_marker_conflict_readback',markers=125,assessment_rows=1000,summary_rows=8,source_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(__file__),scope='Independent DendroPy gene-tree parsing and bitmask incompatibility recreate every projected split/support diagnostic and count. Hybrid tips are removed only from saved topology; support is inherited conservatively, not reestimated. Does not identify the cause of discordance.')
    Path('metadata/saccharomyces_marker_conflict_readback_20260927.json').write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof,indent=2))


if __name__=='__main__':main()
