#!/usr/bin/env python3
"""Independently check all 315 source identities, figure cells and PDF labels."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import re

from PIL import Image
from pypdf import PdfReader

from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def inspect(completion,source_rows,panels,placements,pdf=None):
    cohorts=['full_primary','exclude_sparse_boundary_taxon','exclude_below10_both_alignments',
             'exclude_curated_hybrids','exclude_hybrids_and_uncertain_labels']
    candidates=['coal:profile:uncontracted','coal:profile:sh_alrt10','coal:profile:sh_alrt80',
                'coal:mafft:uncontracted','coal:mafft:sh_alrt10','coal:mafft:sh_alrt80']
    references=['concat:profile_profile:ml','concat:profile_profile:consensus',
                'concat:profile_mafft:ml','concat:profile_mafft:consensus',
                'concat:mafft_profile:ml','concat:mafft_profile:consensus',
                'concat:mafft_mafft:ml','concat:mafft_mafft:consensus']
    summaries={r['cohort']:r for r in completion['cohort_summaries']}
    assert set(summaries)==set(cohorts) and len(panels)==5
    index={}
    for r in source_rows:
        key=(r['cohort'],frozenset([r['view_a'],r['view_b']]))
        assert len(key[1])==2 and key not in index
        index[key]=r
    expected={(c,frozenset([a,b])) for c in cohorts for a in candidates for b in references}
    expected|={(c,frozenset([a,b])) for c in cohorts for i,a in enumerate(candidates) for b in candidates[:i]}
    assert set(index)==expected and len(index)==315 and len(placements)==315
    placed={}
    for r in placements:
        key=(r['cohort'],frozenset([r['view_a'],r['view_b']]))
        assert key not in placed and key in index
        assert set(r)==set(index[key])|{'panel','row','column'}
        assert all(r[field]==value for field,value in index[key].items())
        placed[key]=r
    assert set(placed)==expected
    label_counts=[Counter(),Counter()];cross_count=internal_count=0
    for ci,cohort in enumerate(cohorts):
        p,s=panels[ci],summaries[cohort]
        assert p['cohort']==cohort and p['taxa']==s['taxa'] and p['roles']==s['roles']
        assert sum(p['roles'].values())==p['taxa']
        assert p['coalescent_views']==candidates and p['reference_views']==references
        for field in ['shared_all_coalescent','shared_all_references','shared_all_views']:
            assert p[field]==s[field] and 0<=p[field]<=p['taxa']-3
        for field,width in [('cross_rf',8),('internal_rf_lower_triangle',6)]:
            assert len(p[field])==6 and all(len(row)==width for row in p[field])
        for i,a in enumerate(candidates):
            for panel,other,matrix in [('coalescent_reference',references,p['cross_rf']),
                                       ('coalescent_internal',candidates,p['internal_rf_lower_triangle'])]:
                for j,b in enumerate(other):
                    if panel=='coalescent_internal' and j>=i:
                        assert matrix[i][j] is None
                        continue
                    key=(cohort,frozenset([a,b]));r=index[key];q=placed[key]
                    assert q['panel']==panel and int(q['row'])==i and int(q['column'])==j
                    assert r['comparison_scope']==('coalescent_vs_concatenated' if panel=='coalescent_reference'
                                                   else 'coalescent_alignment_support_sensitivity')
                    assert int(r['shared_splits'])+int(r['unique_splits_a'])==p['taxa']-3
                    assert int(r['shared_splits'])+int(r['unique_splits_b'])==p['taxa']-3
                    rf=int(r['unique_splits_a'])+int(r['unique_splits_b'])
                    assert rf==int(r['rf_distance'])
                    value=rf/(2*(p['taxa']-3))
                    assert 0<=value<=1 and abs(value-float(r['normalized_rf']))<1e-12
                    assert matrix[i][j] is not None and abs(matrix[i][j]-value)<1e-12
                    label_counts[panel=='coalescent_internal'][f'{value:.3f}']+=1
                    if panel=='coalescent_reference':cross_count+=1
                    else:internal_count+=1
    assert cross_count==240 and internal_count==75
    if pdf is not None:
        reader=PdfReader(pdf);assert len(reader.pages)==2
        for page,counts in zip(reader.pages,label_counts):
            # Colorbar ticks use different formatting; these are the 3-decimal cell labels.
            text=page.extract_text();actual=Counter(re.findall(r'(?<![\d.])\d\.\d{3}(?!\d)',text))
            assert actual==counts,(actual-counts,counts-actual)
            assert all(str(summaries[c]['taxa'])+' entries' in text for c in cohorts)
            assert 'No preferred tree' in text and 'SH-aLRT contraction' in text
    return dict(comparisons=315,cross_reference_cells=cross_count,internal_candidate_cells=internal_count,
                panels=10,pdf_pages=2,tree_views=70,cohorts=5)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();plan=json.loads(args.plan.read_text());root=Path(plan['output'])
    bindings=dict(plan['pins']);bind(bindings,args.plan)
    completion=json.loads(Path(plan['completion']).read_text())
    assert completion['status']=='complete_verified_full_coalescent_reference_tree_comparisons'
    assert completion['comparison_rows']==315 and completion['tree_views']==70 and completion['cohorts']==5
    assert completion['exact_process_journals_checked']==2
    bind(bindings,plan['completion']);bind(bindings,completion['full_hash_archive'],completion['full_hash_archive_sha256'])
    archive=json.loads(Path(completion['full_hash_archive']).read_text());assert len(archive['services'])==2
    for path,digest in archive['source_hashes'].items():bind(bindings,path,digest)
    bind(bindings,completion['producer_receipt'],completion['producer_receipt_sha256'])
    bind(bindings,completion['independent_readback'],completion['independent_readback_sha256'])
    native=json.loads(Path(completion['producer_receipt']).read_text())
    prior_reader=json.loads(Path(completion['independent_readback']).read_text())
    assert prior_reader['producer_receipt_sha256']==sha(completion['producer_receipt'])
    assert native['cohort_summaries']==prior_reader['cohort_summaries']==completion['cohort_summaries']
    for name,digest in native['artifacts'].items():bind(bindings,Path(completion['producer_receipt']).parent/name,digest)
    receipt_path=root/'receipt.json';receipt=json.loads(receipt_path.read_text());bind(bindings,receipt_path)
    assert receipt['status']=='complete_full_coalescent_reference_comparison_figure_pending_cell_readback'
    assert receipt['plan_sha256']==sha(args.plan) and receipt['comparison_completion_sha256']==sha(plan['completion'])
    assert set(receipt['artifacts'])=={'panel_data.json','figure_pair_cells.tsv','coalescent_reference_rf.png',
                                      'coalescent_internal_rf.png','coalescent_reference_sensitivity.pdf'}
    for path,digest in receipt['source_hashes'].items():bind(bindings,path,digest)
    for name,digest in receipt['artifacts'].items():bind(bindings,root/name,digest)
    verify(bindings)
    with (Path(completion['producer_receipt']).parent/'comparisons.tsv').open() as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    with (root/'figure_pair_cells.tsv').open() as f:placements=list(csv.DictReader(f,delimiter='\t'))
    summary=inspect(completion,rows,json.loads((root/'panel_data.json').read_text()),placements,
                    root/'coalescent_reference_sensitivity.pdf')
    assert all(receipt[k]==value for k,value in summary.items())
    for name in ['coalescent_reference_rf.png','coalescent_internal_rf.png']:
        with Image.open(root/name) as im:
            assert im.format=='PNG' and im.size==(2196,2376)
            im.verify()
    verify(bindings)
    result=dict(status='passed_full_315_coalescent_comparison_figure_cell_and_pdf_readback',**summary,
        plan_sha256=sha(args.plan),producer_receipt_sha256=sha(receipt_path),source_hashes=bindings,
        visual_inspection_required=True,scientific_eligibility=False,scope=plan['scope'])
    with args.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if isinstance(v,(str,int,bool))}),flush=True)


if __name__=='__main__':main()
