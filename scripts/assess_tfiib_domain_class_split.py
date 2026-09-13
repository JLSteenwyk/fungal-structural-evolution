#!/usr/bin/env python3
"""Describe domain-detection classes across audited unrooted TFIIB tree splits."""
import argparse,json,math
from pathlib import Path
from Bio import Phylo
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT,sha,read_table
from prepare_paired_phylogenetic_inputs import write_table


def masks(tree,index):
    names=[t.name for t in tree.get_terminals()]
    if len(names)!=len(set(names)) or set(names)!=set(index):raise ValueError('Tree tip universe differs')
    result={}
    for c in tree.find_clades(order='postorder'):
        if c.branch_length is not None and (not math.isfinite(c.branch_length) or c.branch_length<0):raise ValueError('Invalid branch length')
        result[c]=(1<<index[c.name]) if c.is_terminal() else sum(result[x] for x in c.clades)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new assessment output')
    base=ROOT/'results/phylogeny/tfiib-domain-tree-v1';r=json.loads((base/'receipt.json').read_text());config=json.loads((base/'run_config.json').read_text())
    if sha(base/'run_config.json')!=r['config_sha256'] or sha(base/'tree.treefile')!=r['tree_sha256'] or sha(base/'tree.ufboot')!=r['bootstrap_sha256']:raise ValueError('Changed completed tree artifacts')
    annotations=ROOT/'results/phylogeny/tfiib-tip-annotations-v1';ar=checked_receipt(annotations)
    if ar['input_receipt_sha256']!=config['source_receipt_sha256']:raise ValueError('Annotation/tree input lineage differs')
    rows=read_table(annotations/'tree_tip_annotations.tsv');tips={x['entry_id']:x for x in rows};index={t:i for i,t in enumerate(sorted(tips))};full=(1<<len(index))-1
    positive=sum(1<<index[t] for t,row in tips.items() if row['BRF1_evidence_class']=='BRF1_detected');negative=full^positive
    tree=Phylo.read(base/'tree.treefile','newick');cm=masks(tree,index);branches=[];branch_masks={}
    for c,side in cm.items():
        n=side.bit_count()
        if n<2 or n>len(index)-2:continue
        label=c.name if c.name is not None else str(c.confidence)
        support=['',''] if label=='None' else [float(x) for x in label.split('/')]
        if len(support)!=2 or any(x!='' and (not math.isfinite(x) or not 0<=x<=100) for x in support):raise ValueError('Invalid dual support')
        errors=(side&negative).bit_count()+(positive&(full^side)).bit_count()
        reverse=(side&positive).bit_count()+(negative&(full^side)).bit_count()
        if reverse<errors:side=full^side;errors=reverse
        signature=hex(min(side,full^side));branch_masks[signature]=side
        branches.append({'split_id':signature,'BRF1_enriched_side_tips':side.bit_count(),'detected_on_enriched_side':(side&positive).bit_count(),'undetected_on_enriched_side':(side&negative).bit_count(),'detected_on_other_side':(positive&(full^side)).bit_count(),'undetected_on_other_side':(negative&(full^side)).bit_count(),'minimum_class_disagreements':errors,'SH_aLRT':support[0],'UFBoot':support[1],'branch_length':c.branch_length})
    if len({b['split_id'] for b in branches})!=len(branches):raise ValueError('Duplicate internal split')
    best=min(branches,key=lambda b:(b['minimum_class_disagreements'],b['split_id']));chosen=branch_masks[best['split_id']];count=0;observed=0
    for bt in Phylo.parse(base/'tree.ufboot','newick'):
        bm=masks(bt,index);splitset={min(s,full^s) for s in bm.values()};observed+=min(chosen,full^chosen) in splitset;count+=1
    if count!=1000 or count!=r['bootstrap_trees']:raise ValueError('Incomplete bootstrap output')
    tip_rows=[]
    for t,row in sorted(tips.items()):
        inside=bool(chosen&(1<<index[t]));detected=row['BRF1_evidence_class']=='BRF1_detected'
        tip_rows.append(row|{'best_split_side':'BRF1_enriched' if inside else 'BRF1_undetected_enriched','class_disagrees_with_split':inside!=detected})
    selected={side:sum(x['selected_focal_marker']=='True' and x['best_split_side']==side for x in tip_rows) for side in ['BRF1_enriched','BRF1_undetected_enriched']}
    a.output.mkdir(parents=True);write_table(a.output/'internal_split_class_counts.tsv',branches);write_table(a.output/'best_split_tip_annotations.tsv',tip_rows)
    result={'status':'complete_exploratory_tfiib_class_split_assessment','tree_receipt_sha256':sha(base/'receipt.json'),'annotation_receipt_sha256':sha(annotations/'receipt.json'),'script_sha256':sha(Path(__file__)),'tips':len(tips),'internal_splits':len(branches),'internal_splits_without_support_labels':sum(b['SH_aLRT']=='' for b in branches),'best_split':best,'bootstrap_trees_checked':count,'best_split_observed_in_bootstraps':observed,'selected_focal_marker_side_counts':selected,'interpretation':'Split chosen exploratorily by minimum disagreement with BRF1 detection classes, not a prespecified association test. Bootstrap recount measures this split support, not confidence in functional labels or independent duplication events. Unrooted topology; repeat sensitivity and species/gene reconciliation remain required.','artifacts':{p.name:sha(p) for p in a.output.iterdir() if p.is_file()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='artifacts'},indent=2))


if __name__=='__main__':main()
