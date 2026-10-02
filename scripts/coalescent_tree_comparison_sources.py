"""Closed numerical and same-taxon reference sources for all coalescent comparisons."""
import itertools
import json
from pathlib import Path

from four_run_pmsf_sources import load_sources
from retained_tree_comparison_sources import load as load_retained
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha

METRICS=['branch_length','branch_length_unit','local_posterior','effective_genes','empirical_ufb_percent','sh_alrt_percent']
PAIR_FIELDS=['cohort','view_a','view_b','comparison_scope','support_rule_a','support_rule_b',
    'shared_splits','unique_splits_a','unique_splits_b','rf_distance','normalized_rf',
    'incompatible_pairs','both_support_rules_met_pairs']
PRESENCE_FIELDS=['cohort','view','family','support_rule','split_taxa_json','present',*METRICS,'support_rule_met']
CONFLICT_FIELDS=['cohort','view_a','view_b','split_a_taxa_json','split_b_taxa_json',
    'support_rule_a','support_rule_b','support_rule_met_a','support_rule_met_b',
    'both_support_rules_met','witness_quartet']
BOUNDARY_FIELDS=['cohort','view','ingroup','outgroup','boundary_taxa_json','boundary_present',
    'support_rule',*METRICS,'support_rule_met']
SUMMARY_FIELDS=['cohorts','tree_views','comparison_rows','split_presence_rows','incompatible_pairs',
                'role_boundary_rows','boundary_views_with_role_split','cohort_summaries']


def canonical(side,taxa):
    return min(tuple(sorted(side)),tuple(sorted(set(taxa)-set(side))),key=lambda s:(len(s),s))


def rule_met(row,rule):
    if rule=='local_PP95':return row['local_posterior']>=.95
    if rule=='original_SH80_and_empirical_UFB95':
        return row['sh_alrt_percent']>=80 and row['empirical_ufb_percent']>=95
    assert rule in ['original_consensus_empirical_UFB95_SH_unavailable','projected_empirical_UFB95_SH_unavailable']
    assert row['sh_alrt_percent'] is None
    return row['empirical_ufb_percent']>=95


def pair_grid(views):
    return [(a,b) for a,b in itertools.combinations(sorted(views),2)
            if views[a]['family']=='coalescent' or views[b]['family']=='coalescent']


def comparison_scope(a,b):
    return 'coalescent_vs_concatenated' if a['family']!=b['family'] else 'coalescent_alignment_support_sensitivity'


def load(plan,plan_path):
    bindings=dict(plan['pins']);bind(bindings,plan_path)
    rp=Path(plan['retained_reference_plan']);retained_plan=json.loads(rp.read_text())
    retained,universe,prior=load_retained(retained_plan,rp)
    for p,d in prior.items():bind(bindings,p,d)
    bp=Path(plan['primary_reference_plan']);baseline_plan=json.loads(bp.read_text())
    primary,manifest,prior=load_sources(baseline_plan,bp)
    for p,d in prior.items():bind(bindings,p,d)
    primary_closed=json.loads(Path(plan['primary_reference_completion']).read_text())
    assert primary_closed['status']=='complete_verified_full_four_run_pmsf_ML_consensus_sensitivity'
    assert len(primary_closed['services'])==2
    bind(bindings,plan['primary_reference_completion'])
    for p,d in primary_closed['source_hashes'].items():bind(bindings,p,d)
    groups={'full_primary':dict(taxa=universe,outgroups={r['taxon_id'] for r in manifest if r['study_role']=='outgroup'},
        roles=dict(ingroup=501,outgroup=25),views={})}
    def reference_view(run,kind,rows,projected,cohort):
        metrics={}
        for row in rows:
            if not projected and row['tree']!=kind:continue
            side=tuple(json.loads(row['split_taxa_json']))
            key=canonical(side,groups[cohort]['taxa'])
            assert len(key)>=2 and key not in metrics
            metrics[key]=dict(branch_length=float(row['projected_branch_length_sum'] if projected else row['branch_length']),
                branch_length_unit='original_path_sum_AA_substitutions_per_site' if projected else 'AA_substitutions_per_site',
                local_posterior=None,effective_genes=None,
                empirical_ufb_percent=float(row['empirical_projected_ufboot_percent'] if projected else row['empirical_ufboot_percent']),
                sh_alrt_percent=None if projected or kind=='consensus' else float(row['sh_alrt_percent']))
        rule=('projected_empirical_UFB95_SH_unavailable' if projected else
              'original_SH80_and_empirical_UFB95' if kind=='ml' else 'original_consensus_empirical_UFB95_SH_unavailable')
        return dict(family='concatenated',kind=kind,run=run,source_tree=str(Path(primary[run]['spec']['run'])/('pmsf.treefile' if kind=='ml' else 'pmsf.contree')),
                    prune=projected,support_rule=rule,splits=metrics)
    for label,source in primary.items():
        for kind in ['ml','consensus']:
            groups['full_primary']['views']['concat:'+label+':'+kind]=reference_view(label,kind,source['rows'],False,'full_primary')
    for cohort,group in retained.items():
        groups[cohort]=dict(taxa=group['taxa'],outgroups=group['outgroups'],roles=group['roles'],views={})
        for label,view in group['views'].items():
            groups[cohort]['views']['concat:'+label]=reference_view(view['run'],view['kind'],list(view['splits'].values()),True,cohort)
    cp=Path(plan['coalescent_completion']);completion=json.loads(cp.read_text())
    if plan['mode']=='full_batch':
        assert completion['status']=='complete_verified_full_native_coalescent_quartet_numerics'
        assert completion['cases']==30 and completion['branch_gene_states']==1938750
        assert completion['exact_process_journals_checked']==2
        readback_path=completion['independent_readback']
    else:
        assert plan['mode']=='complete_named_case_preflight'
        assert completion['status']=='complete_verified_named_native_coalescent_global_local_quartet_preflight'
        assert completion['cases']==1 and completion['exact_process_journals_checked']==1
        readback_path=completion['independent_readback']
    bind(bindings,cp);bind(bindings,completion['full_hash_archive'],completion['full_hash_archive_sha256'])
    closed=json.loads(Path(completion['full_hash_archive']).read_text())
    assert len(closed['services'])==completion['exact_process_journals_checked']
    for p,d in closed['source_hashes'].items():bind(bindings,p,d)
    readback=json.loads(Path(readback_path).read_text());bind(bindings,readback_path,completion['independent_readback_sha256'])
    native_plan=json.loads(Path(plan['native_plan']).read_text());bind(bindings,plan['native_plan'])
    inputs=json.loads(Path(json.loads(Path(native_plan['input_plan']).read_text())['output']).joinpath('receipt.json').read_text())
    cases={r['case']:r for r in inputs['summaries']}
    chosen_cohorts=set()
    for summary in readback['summaries']:
        spec=cases[summary['case']];cohort=spec['cohort'];chosen_cohorts.add(cohort)
        assert groups[cohort]['roles']==spec['roles'] and len(groups[cohort]['taxa'])==spec['expected_taxa']
        local_root=Path(summary['local_readback_receipt']).parent
        branches=json.loads((local_root/'branches.json').read_text())
        bind(bindings,local_root/'branches.json')
        metrics={}
        for row in branches:
            key=canonical(row['canonical_side'],groups[cohort]['taxa']);assert key not in metrics
            metrics[key]=dict(branch_length=row['native_map_length'],branch_length_unit='MAP_coalescent_units',
                local_posterior=row['native']['pp1'],effective_genes=row['native']['EN'],empirical_ufb_percent=None,sh_alrt_percent=None)
        groups[cohort]['views']['coal:'+spec['alignment']+':'+spec['support_policy']]=dict(family='coalescent',kind='estimated',
            source_tree=str(Path(native_plan['output'])/spec['case']/'species.tree'),prune=False,support_rule='local_PP95',
            case=spec['case'],alignment=spec['alignment'],support_policy=spec['support_policy'],splits=metrics)
    groups={c:g for c,g in groups.items() if c in chosen_cohorts}
    for group in groups.values():
        assert all(len(v['splits'])==len(group['taxa'])-3 for v in group['views'].values())
        assert len([v for v in group['views'].values() if v['family']=='concatenated'])==8
        if plan['mode']=='full_batch':assert len(group['views'])==14 and len(pair_grid(group['views']))==63
    if plan['mode']=='full_batch':assert len(groups)==5 and sum(map(lambda g:len(pair_grid(g['views'])),groups.values()))==315
    else:assert len(groups)==1 and sum(len(g['views']) for g in groups.values())==9
    verify(bindings)
    return groups,universe,bindings
