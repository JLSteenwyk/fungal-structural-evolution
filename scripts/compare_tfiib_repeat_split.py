#!/usr/bin/env python3
"""Test the frozen combined-domain bipartition in a completed positional-repeat tree."""
import argparse,json,math
from pathlib import Path
from Bio import Phylo,SeqIO
from assess_tfiib_domain_class_split import masks
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT,sha,read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repeat',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable repeat comparison')
    r=json.loads((a.repeat/'receipt.json').read_text());c=json.loads((a.repeat/'run_config.json').read_text())
    if r['status']!='complete_supported_tfiib_repeat_tree' or r['config_sha256']!=sha(a.repeat/'run_config.json') or r['tree_sha256']!=sha(a.repeat/'tree.treefile') or r['bootstrap_sha256']!=sha(a.repeat/'tree.ufboot'):raise ValueError('Changed/incomplete repeat tree')
    source=ROOT/'results/phylogeny/tfiib-repeat-sensitivity-inputs-v1';checked_receipt(source)
    if c['source_receipt_sha256']!=sha(source/'receipt.json'):raise ValueError('Repeat input lineage differs')
    reference=ROOT/'results/phylogeny/tfiib-class-split-v1';ref=checked_receipt(reference)
    rows=read_table(reference/'best_split_tip_annotations.tsv');tips={x['entry_id']:x for x in rows};index={t:i for i,t in enumerate(sorted(tips))};full=(1<<len(tips))-1
    target=sum(1<<index[t] for t,row in tips.items() if row['best_split_side']=='BRF1_enriched');target_id=min(target,full^target)
    if hex(target_id)!=ref['best_split']['split_id']:raise ValueError('Reference split mask differs')
    positive=sum(1<<index[t] for t,row in tips.items() if row['BRF1_evidence_class']=='BRF1_detected');negative=full^positive
    tree=Phylo.read(a.repeat/'tree.treefile','newick');cm=masks(tree,index);branches=[]
    for clade,side in cm.items():
        n=side.bit_count()
        if not 2<=n<=len(tips)-2:continue
        label=clade.name if clade.name is not None else str(clade.confidence)
        support=['',''] if label=='None' else [float(x) for x in label.split('/')]
        if len(support)!=2 or any(x!='' and (not math.isfinite(x) or not 0<=x<=100) for x in support):raise ValueError('Malformed support')
        errors=(side&negative).bit_count()+(positive&(full^side)).bit_count();reverse=(side&positive).bit_count()+(negative&(full^side)).bit_count()
        if reverse<errors:side=full^side;errors=reverse
        branches.append({'split_id':hex(min(side,full^side)),'BRF1_enriched_side_tips':side.bit_count(),'detected_on_enriched_side':(side&positive).bit_count(),'undetected_on_enriched_side':(side&negative).bit_count(),'detected_on_other_side':(positive&(full^side)).bit_count(),'undetected_on_other_side':(negative&(full^side)).bit_count(),'minimum_class_disagreements':errors,'SH_aLRT':support[0],'UFBoot':support[1],'branch_length':clade.branch_length})
    if len({b['split_id'] for b in branches})!=len(branches):raise ValueError('Duplicate internal split')
    best=min(branches,key=lambda b:(b['minimum_class_disagreements'],b['split_id']));best_id=int(best['split_id'],16);found=[b for b in branches if int(b['split_id'],16)==target_id];count=0;reference_count=0;best_count=0
    for bt in Phylo.parse(a.repeat/'tree.ufboot','newick'):
        bm=masks(bt,index);splits={min(s,full^s) for s in bm.values()};reference_count+=target_id in splits;best_count+=best_id in splits;count+=1
    if count!=1000 or count!=r['bootstrap_trees']:raise ValueError('Bootstrap count differs')
    side_by_tip=[]
    for t,row in sorted(tips.items()):
        side_by_tip.append({**row,'repeat_best_split_canonical_side':bool(best_id&(1<<index[t]))})
    a.output.mkdir(parents=True);write_table(a.output/'internal_split_class_counts.tsv',branches);write_table(a.output/'tip_split_membership.tsv',side_by_tip)
    result={'status':'complete_tfiib_repeat_reference_split_assessment','repeat_order':r['repeat_order'],'repeat_receipt_sha256':sha(a.repeat/'receipt.json'),'reference_assessment_receipt_sha256':sha(reference/'receipt.json'),'script_sha256':sha(Path(__file__)),'mask_helper_sha256':sha(Path(__file__).with_name('assess_tfiib_domain_class_split.py')),'tips':len(tips),'bootstrap_trees_checked':count,'combined_reference_split_present_in_ML_tree':bool(found),'combined_reference_split_branch':found[0] if found else None,'combined_reference_split_bootstrap_count':reference_count,'repeat_exploratory_best_class_split':best,'repeat_best_class_split_bootstrap_count':best_count,'best_class_split_equals_combined_reference':best_id==target_id,'interpretation':'Exact unrooted bipartition comparison on identical gene-copy tips. Reference split fixed from the combined-domain analysis; repeat best class split is separately exploratory. Repeat data overlap the combined alignment, so agreement is not independent replication. Split support does not establish biochemical labels, root, duplication timing or orthology. Other repeat and species/gene reconciliation remain required.','artifacts':{f.name:sha(f) for f in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='artifacts'},indent=2))

if __name__=='__main__':main()
