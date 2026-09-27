#!/usr/bin/env python3
"""Compare all audited marker splits with the Saccharomyces PMSF conflict."""
import csv
import hashlib
import json
from collections import defaultdict, Counter
from pathlib import Path
from Bio import Phylo


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def canonical(side,taxa):return min(tuple(sorted(side)),tuple(sorted(taxa-side)),key=lambda s:(len(s),s))


def main():
    source=Path('results/phylogeny/marker-tree-support-complete-v1')
    audit=Path('results/phylogeny/marker-tree-support-complete-readback-v1/receipt.json')
    sr=json.loads((source/'receipt.json').read_text());proof=json.loads(audit.read_text())
    assert sr['completed_markers']==sr['planned_markers']==125 and not sr['pending_markers']
    assert proof['status']=='passed_all_snapshot_input_and_graph_split_readbacks' and proof['snapshot_receipt_sha256']==sha(source/'receipt.json')
    for name,digest in sr['artifacts'].items():assert sha(source/name)==digest
    cr=Path('results/phylogeny/pmsf-three-run-topology-sensitivity-20260927-v1')
    cp=json.loads(Path('metadata/pmsf_three_run_topology_readback_20260927.json').read_text())
    assert cp['source_receipt_sha256']==sha(cr/'receipt.json')
    rr=json.loads((cr/'receipt.json').read_text());assert rr['artifacts']['conflicts.tsv']==sha(cr/'conflicts.tsv')
    conflicts=[r for r in read(cr/'conflicts.tsv') if r['both_meet_existing_support_label']=='True']
    alternatives={(r['split_a_taxa_json'],r['split_b_taxa_json']) for r in conflicts};assert len(alternatives)==1
    x,y=next(iter(alternatives));targets={'profile_alignment':set(json.loads(x)),'mafft_alignment':set(json.loads(y))}
    support=defaultdict(list)
    for r in read(source/'branch_support.tsv'):support[r['marker']].append(r)
    output=[];collapsed=[]
    for item in sr['inputs']:
        marker=item['marker'];treepath=Path('results/phylogeny/marker-gene-trees-v2')/marker/'tree.treefile'
        assert sha(treepath)==item['tree_sha256']
        taxa={n.name for n in Phylo.read(treepath,'newick').get_terminals()}
        for mode,omitted in [('all_available',set()),('hybrid_tips_removed',{'F27292','F332112'})]:
            universe=taxa-omitted
            projected=defaultdict(list)
            for edge in support[marker]:
                side=set(edge['smaller_side_taxa'].split(';'))&universe
                split=canonical(side,universe)
                if len(split)>=2:projected[split].append(None if edge['sh_alrt_percent']=='' else float(edge['sh_alrt_percent']))
            scores={s:min(v) if all(z is not None for z in v) else None for s,v in projected.items()}
            for split,vals in projected.items():
                if len(vals)>1:collapsed.append(dict(marker=marker,mode=mode,split_taxa=';'.join(split),source_edges=len(vals),minimum_source_sh_alrt=scores[split] if scores[split] is not None else '',missing_source_support=sum(v is None for v in vals)))
            for label,target in targets.items():
                split=canonical(target&universe,universe);xs=set(split)
                exact=scores.get(split);incompatible=[]
                for other,score in scores.items():
                    ys=set(other)
                    if xs&ys and xs-ys and ys-xs and universe-(xs|ys):incompatible.append((score,other))
                known=[v for v,s in incompatible if v is not None]
                maxconf=max(known) if known else None
                for cutoff in (80,95):
                    status='uninformative_taxon_coverage' if len(split)<2 else 'supported_concordance' if exact is not None and exact>=cutoff else 'supported_conflict' if maxconf is not None and maxconf>=cutoff else 'unresolved'
                    output.append(dict(marker=marker,mode=mode,alternative=label,retained_taxa=len(universe),removed_hybrid_tips=len(taxa&omitted),target_taxa_present=';'.join(sorted(target&universe)),restricted_split=';'.join(split),exact_split_present=split in projected,exact_minimum_source_sh_alrt=exact if exact is not None else '',exact_source_edges=len(projected.get(split,[])),maximum_conflicting_minimum_source_sh_alrt=maxconf if maxconf is not None else '',incompatible_splits=len(incompatible),incompatible_splits_without_complete_support=sum(v is None for v,s in incompatible),sh_alrt_cutoff=cutoff,status=status))
    summary=[]
    for mode in ['all_available','hybrid_tips_removed']:
        for label in targets:
            for cutoff in [80,95]:
                subset=[r for r in output if (r['mode'],r['alternative'],r['sh_alrt_cutoff'])==(mode,label,cutoff)]
                counts=Counter(r['status'] for r in subset);assert len(subset)==125
                summary.append(dict(mode=mode,alternative=label,sh_alrt_cutoff=cutoff,markers=125,**{k:counts[k] for k in ['supported_concordance','supported_conflict','unresolved','uninformative_taxon_coverage']}))
    root=Path('results/phylogeny/saccharomyces-marker-conflict-20260927-v1');root.mkdir(parents=True,exist_ok=False)
    for name,data in [('marker_assessments.tsv',output),('summary.tsv',summary)]:
        with (root/name).open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');w.writeheader();w.writerows(data)
    report=dict(status='complete_saccharomyces_marker_split_assessment_pending_readback',markers=125,rows=len(output),collapsed_split_records=len(collapsed),source_support_receipt_sha256=sha(source/'receipt.json'),source_support_readback_sha256=sha(audit),source_conflict_receipt_sha256=sha(cr/'receipt.json'),script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in root.iterdir()},scope='Full 125-marker saved-tree diagnostic. Hybrid-tip removal only projects existing unrooted splits; it does not refit alignments, profiles, trees or support. A collapsed path qualifies only if every source edge has support and its minimum passes the descriptive SH-aLRT cutoff. SH-aLRT is not bootstrap probability; rows are correlated and missing support remains unresolved. No cause of discordance or preferred species tree established.')
    (root/'receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
