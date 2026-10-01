"""Closed projection source I/O for complete same-cohort tree comparisons."""
import csv
import json
from pathlib import Path
from retained_taxon_projection_sources import load as load_projection
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha

PAIR_FIELDS=['policy','view_a','view_b','comparison_scope','support_rule',
             'shared_internal_splits','unique_internal_splits_a','unique_internal_splits_b',
             'rf_distance','normalized_rf','incompatible_split_pairs','both_UFB95_incompatible_pairs']
PRESENCE_FIELDS=['policy','view','split_taxa_json','present','projected_branch_length_sum',
                 'empirical_projected_ufboot_percent','sh_alrt_percent']
CONFLICT_FIELDS=['policy','view_a','view_b','split_a_taxa_json','split_b_taxa_json',
                 'empirical_projected_ufboot_percent_a','empirical_projected_ufboot_percent_b',
                 'both_UFB95','witness_quartet']
SUMMARY_FIELDS=['policies','tree_views','comparison_rows','split_presence_rows','incompatible_pairs']
RULE='empirical_projected_UFB95_SH_aLRT_unavailable'


def load(plan,plan_path):
    bindings=dict(plan['pins']);bind(bindings,plan_path)
    pp=Path(plan['projection_plan']);projection=json.loads(pp.read_text())
    sources,universe,policies,prior=load_projection(projection,pp)
    for p,d in prior.items():bind(bindings,p,d)
    closed_path=Path(plan['projection_completion']);closed=json.loads(closed_path.read_text())
    assert closed['status']=='complete_verified_full_retained_taxon_baseline_projections' and closed['exact_process_journals_checked']==2
    archive=Path(closed['full_hash_archive']);bind(bindings,archive,closed['full_hash_archive_sha256'])
    proof=json.loads(archive.read_text());assert proof['status']=='complete_verified_retained_taxon_projection_archive' and len(proof['services'])==2
    assert closed['bound_source_hashes']==len(proof['source_hashes'])
    for p,d in proof['source_hashes'].items():bind(bindings,p,d)
    bind(bindings,closed_path)
    root=Path(projection['output']);rp=root/'receipt.json';r=json.loads(rp.read_text())
    assert closed['producer_receipt_sha256']==sha(rp) and closed['independent_readback_sha256']==sha(closed['independent_readback'])
    assert r['status']=='complete_full_retained_taxon_baseline_projections_pending_independent_readback' and r['plan_sha256']==sha(pp)
    assert r['tree_views']==32 and r['edge_rows']==33088 and r['projected_bootstrap_tree_states']==16000
    for name,digest in r['artifacts'].items():bind(bindings,root/name,digest)
    groups={p:dict(**spec,views={}) for p,spec in policies.items()}
    for policy,group in groups.items():
        for label,source in sources.items():
            for kind in ['ml','consensus']:
                group['views'][label+':'+kind]=dict(run=label,kind=kind,splits={},source=source['spec'])
    with (root/'projected_tree_edges.tsv').open() as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            group=groups[row['policy']];view=group['views'][row['baseline_run']+':'+row['tree_type']]
            key=tuple(json.loads(row['split_taxa_json']))
            if len(key)==1:continue
            assert len(key)==len(set(key)) and set(key)<=group['taxa'] and row['sh_alrt_percent']==''
            assert key==min(tuple(sorted(key)),tuple(sorted(group['taxa'].difference(key))),key=lambda x:(len(x),x))
            assert key not in view['splits'];view['splits'][key]=row
    assert len(groups)==4
    for group in groups.values():
        assert len(group['views'])==8 and all(len(v['splits'])==len(group['taxa'])-3 for v in group['views'].values())
    verify(bindings)
    return groups,universe,bindings
