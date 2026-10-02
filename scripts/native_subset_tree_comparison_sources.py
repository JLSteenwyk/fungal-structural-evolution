"""Closed actual subset fits and all matching coalescent/projected sources."""
import csv
import itertools
import json
from pathlib import Path

from coalescent_tree_comparison_sources import load as load_candidates,canonical,METRICS
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha

PAIR_FIELDS=['cohort','view_a','view_b','comparison_scope','support_rule_a','support_rule_b',
    'splits_a','splits_b','possible_internal_splits','shared_splits','unique_splits_a','unique_splits_b',
    'rf_distance','normalized_rf','normalized_rf_by_observed_splits','incompatible_pairs','both_support_rules_met_pairs']
PRESENCE_FIELDS=['cohort','view','family','support_rule','split_taxa_json','present',*METRICS,'support_rule_met']
CONFLICT_FIELDS=['cohort','view_a','view_b','split_a_taxa_json','split_b_taxa_json','support_rule_a','support_rule_b',
    'support_rule_met_a','support_rule_met_b','both_support_rules_met','witness_quartet']
BOUNDARY_FIELDS=['cohort','view','ingroup','outgroup','boundary_taxa_json','boundary_present','support_rule',*METRICS,'support_rule_met']
SUMMARY_FIELDS=['cohorts','tree_views','comparison_rows','split_presence_rows','incompatible_pairs',
    'role_boundary_rows','boundary_views_with_role_split','cohort_summaries']


def supported(row,rule):
    if rule=='local_PP95':return row['local_posterior']>=.95
    if rule=='native_SH80_and_empirical_UFB95':return row['sh_alrt_percent']>=80 and row['empirical_ufb_percent']>=95
    assert rule in ['native_consensus_empirical_UFB95_SH_unavailable','projected_empirical_UFB95_SH_unavailable']
    assert row['sh_alrt_percent'] is None
    return row['empirical_ufb_percent']>=95


def pair_grid(views):
    return [(a,b) for a,b in itertools.combinations(sorted(views),2)
            if views[a]['family']=='native_pmsf' or views[b]['family']=='native_pmsf']


def comparison_scope(a,b):
    families={a['family'],b['family']};assert 'native_pmsf' in families
    if families=={'native_pmsf'}:return 'native_alignment_guide_view_sensitivity'
    return 'native_vs_coalescent' if 'coalescent' in families else 'native_vs_projected_reference'


def load(plan,plan_path):
    bindings=dict(plan['pins']);bind(bindings,plan_path)
    cp=Path(plan['candidate_plan']);candidate_plan=json.loads(cp.read_text())
    groups,universe,prior=load_candidates(candidate_plan,cp)
    for p,d in prior.items():bind(bindings,p,d)
    complete=json.loads(Path(plan['candidate_completion']).read_text())
    assert complete['status']=='complete_verified_full_coalescent_reference_tree_comparisons'
    assert complete['comparison_rows']==315 and complete['tree_views']==70 and complete['exact_process_journals_checked']==2
    bind(bindings,plan['candidate_completion']);bind(bindings,complete['full_hash_archive'],complete['full_hash_archive_sha256'])
    archive=json.loads(Path(complete['full_hash_archive']).read_text());assert len(archive['services'])==2
    for p,d in archive['source_hashes'].items():bind(bindings,p,d)
    groups.pop('full_primary');assert len(groups)==4
    for g in groups.values():assert len(g['views'])==14 and all(v['prune'] for v in g['views'].values() if v['family']=='concatenated')
    native=json.loads(Path(plan['native_completion']).read_text())
    assert native['status']=='complete_verified_full_native_taxon_pmsf_collection'
    assert native['run_count']==16 and native['tree_views']==32 and native['raw_bootstrap_trees']==16000
    assert native['site_profiles']==902216 and native['exact_process_journals_checked']==2
    bind(bindings,plan['native_completion']);bind(bindings,native['full_hash_archive'],native['full_hash_archive_sha256'])
    proof=json.loads(Path(native['full_hash_archive']).read_text());assert len(proof['services'])==2
    for p,d in proof['source_hashes'].items():bind(bindings,p,d)
    bind(bindings,native['native_batch_receipt'],native['native_batch_receipt_sha256'])
    bind(bindings,native['independent_readback'],native['independent_readback_sha256'])
    raw=json.loads(Path(native['native_batch_receipt']).read_text());reader=json.loads(Path(native['independent_readback']).read_text())
    assert reader['native_batch_receipt_sha256']==sha(native['native_batch_receipt'])
    assert len(raw['runs'])==len(reader['runs'])==16
    for name,digest in reader['artifacts'].items():bind(bindings,Path(native['independent_readback']).parent/name,digest)
    found=set()
    for source,checked in zip(raw['runs'],reader['runs']):
        for field in ['label','policy','alignment','guide_alignment']:assert source[field]==checked[field]
        cohort=source['policy'];g=groups[cohort];assert source['taxa']==checked['taxa']==len(g['taxa'])
        assert checked['bootstrap_trees']==1000
        table=Path(native['independent_readback']).parent/(source['label']+'.support.tsv')
        with table.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        for kind in ['ml','consensus']:
            identity=(cohort,source['alignment'],source['guide_alignment'],kind);assert identity not in found;found.add(identity)
            metrics={}
            for row in rows:
                if row['tree']!=kind:continue
                key=canonical(json.loads(row['split_taxa_json']),g['taxa']);assert key not in metrics and 2<=len(key)<=len(g['taxa'])//2
                metrics[key]=dict(branch_length=float(row['branch_length']),branch_length_unit='AA_substitutions_per_site',
                    local_posterior=None,effective_genes=None,empirical_ufb_percent=float(row['empirical_ufboot_percent']),
                    sh_alrt_percent=float(row['sh_alrt_percent']) if row['sh_alrt_percent'] else None)
            assert len(metrics)<=len(g['taxa'])-3
            if kind=='ml':assert len(metrics)==len(g['taxa'])-3
            name='native:'+source['alignment']+'_'+source['guide_alignment']+':'+kind
            assert name not in g['views']
            g['views'][name]=dict(family='native_pmsf',kind=kind,prune=False,
                source_tree=str(Path(source['run'])/('pmsf.treefile' if kind=='ml' else 'pmsf.contree')),
                support_rule='native_SH80_and_empirical_UFB95' if kind=='ml' else 'native_consensus_empirical_UFB95_SH_unavailable',splits=metrics)
    assert len(found)==32 and all(len(g['views'])==22 and len(pair_grid(g['views']))==140 for g in groups.values())
    verify(bindings)
    return groups,universe,bindings
