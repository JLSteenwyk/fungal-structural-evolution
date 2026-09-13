#!/usr/bin/env python3
"""Attach curated orthology caveats to frozen paired-marker review tables."""
import argparse,json
from collections import Counter
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT,sha,read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use an immutable review snapshot')
    cp=ROOT/'config/marker_orthology_review.json';caveats={r['marker']:r for r in json.loads(cp.read_text())['records']}
    evidence=ROOT/'results/phylogeny/tfiib-class-split-v1';er=checked_receipt(evidence)
    tiprows=read_table(evidence/'best_split_tip_annotations.tsv');selected={}
    for r in tiprows:
        if r['selected_focal_marker']=='True':
            key=(r['taxon_id'],r['protein_id'])
            if key in selected:raise ValueError('Ambiguous selected protein evidence')
            selected[key]=r
    outputs=[];tips=[];pins={}
    for source,input_folder,map_folder in [('AlphaFold','paired-inputs-gdm-expanded-v1','gdm-expanded-v1'),('ESMFold','paired-inputs-esmfold-partial-v1','esmfold-partial-v1')]:
        base=ROOT/'results/phylogeny'/input_folder;mapping=ROOT/'results/structural_markers'/map_folder;receipt=checked_receipt(base);checked_receipt(mapping)
        if receipt['source_receipts']['snapshot']['sha256']!=sha(mapping/'receipt.json'):raise ValueError('Source mapping differs')
        pins[source]=sha(base/'receipt.json');markers=read_table(base/'marker_summary.tsv');ready={r['marker'] for r in markers if r['status']=='ready_for_inference'}
        for r in markers:outputs.append({'source':source,**r,'curated_orthology_caveat':caveats.get(r['marker'],{}).get('status','none_recorded'),'confirmatory_interpretation_status':'withhold_pending_copy_reconciliation' if r['marker'] in caveats else 'general_validation_still_required'})
        eligible={r['taxon_id'] for r in read_table(base/'taxon_coverage.tsv') if r['marker']=='4986044at2759' and r['taxon_eligible']=='True'}
        for r in read_table(mapping/'marker_structure_links.tsv'):
            if r['marker']!='4986044at2759' or r['taxon_id'] not in eligible:continue
            annotation=selected.get((r['taxon_id'],r['protein_id']))
            if annotation and annotation['sequence_id']!='S'+r['sequence_sha256']:raise ValueError('Selected sequence identity differs')
            tips.append({'source':source,'marker':r['marker'],'taxon_id':r['taxon_id'],'protein_id':r['protein_id'],'sequence_sha256':r['sequence_sha256'],'source_marker_ready':r['marker'] in ready,'family_tree_entry':annotation['entry_id'] if annotation else 'unresolved','family_tree_side':annotation['best_split_side'] if annotation else 'outside_eligible_family_tree','BRF1_evidence_class':annotation['BRF1_evidence_class'] if annotation else 'unresolved'})
    a.output.mkdir(parents=True);write_table(a.output/'paired_marker_review.tsv',outputs);write_table(a.output/'flagged_marker_copy_coverage.tsv',tips)
    summary={source:dict(Counter(r['family_tree_side'] for r in tips if r['source']==source)) for source in pins}
    result={'status':'complete_marker_caveat_review_overlay','caveat_config_sha256':sha(cp),'family_split_receipt_sha256':sha(evidence/'receipt.json'),'paired_input_receipt_sha256':pins,'script_sha256':sha(Path(__file__)),'marker_review_rows':len(outputs),'flagged_marker_copy_side_counts':summary,'interpretation':'Read-only review overlay; frozen fits and input masks retained. Flagged marker must not receive confirmatory single-ortholog interpretation without copy reconciliation. Unflagged markers are not certified valid. Exact taxon/protein/sequence joins distinguish eligible family-tree copies from unresolved candidates.','artifacts':{p.name:sha(p) for p in a.output.iterdir() if p.is_file()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
