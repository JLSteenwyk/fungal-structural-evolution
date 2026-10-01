"""Full audited crossed PMSF source/provenance I/O, without split algorithms."""
import csv
import json
from pathlib import Path
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha

PAIR_FIELDS = ['view_a', 'view_b', 'comparison_scope', 'support_rule_a', 'support_rule_b',
               'shared_internal_splits', 'unique_internal_splits_a', 'unique_internal_splits_b',
               'rf_distance', 'normalized_rf', 'incompatible_split_pairs', 'both_supported_incompatible_pairs']
PRESENCE_FIELDS = ['view', 'run', 'tree_type', 'split_taxa_json', 'present', 'branch_length',
                   'sh_alrt', 'empirical_ufb']
CONFLICT_FIELDS = ['view_a', 'view_b', 'split_a_taxa_json', 'split_b_taxa_json',
                   'support_rule_a', 'support_rule_b', 'sh_alrt_a', 'sh_alrt_b',
                   'empirical_ufb_a', 'empirical_ufb_b', 'both_support_criteria_met', 'witness_quartet']
BOUNDARY_FIELDS = ['view', 'ingroup_count', 'outgroup_count', 'boundary_split_present',
                  'sh_alrt', 'empirical_ufb', 'boundary_taxa_json']
SUMMARY_FIELDS = ['taxa', 'runs', 'tree_views', 'internal_splits_per_view',
                  'comparison_rows', 'split_presence_rows', 'incompatible_pairs',
                  'shared_all_ml', 'shared_all_consensus', 'shared_all_eight',
                  'boundary_views_with_role_split']


def load_sources(plan, plan_path):
    bindings = dict(plan['pins']); bind(bindings, plan_path)
    c = json.loads(Path(plan['fourth_completion']).read_text())
    assert c['status'] == 'complete_verified_fourth_crossed_pmsf_profile_tree_and_bootstrap_readback' and len(c['services']) == 2
    for path, digest in c['source_hashes'].items(): bind(bindings, path, digest)
    assert len(plan['runs']) == 4
    sources = {}
    for label, spec in plan['runs'].items():
        run = Path(spec['run']); audit = Path(spec['audit']); rp = run / 'receipt.json'; ap = audit / 'receipt.json'
        r = json.loads(rp.read_text()); a = json.loads(ap.read_text()); config = json.loads((run / 'config.json').read_text())
        assert r['status'] == 'complete_pmsf_execution_pending_full_audit' and r['returncode'] == 0
        assert a['status'] == 'passed_pmsf_profile_tree_and_bootstrap_readback' and a['source_receipt_sha256'] == sha(rp)
        assert a['taxa'] == plan['expected']['taxa'] and a['bootstrap_trees'] == 1000 and a['internal_support_rows'] == 1046
        assert r['config_sha256'] == sha(run / 'config.json') and config['resource_plan']['model'] == 'LG+C20+F+G4'
        assert config['command'][config['command'].index('-m') + 1] == 'LG+C20+F+G4'
        for root, record in [(run, r), (audit, a)]:
            bind(bindings, root / 'receipt.json')
            for name, digest in record['artifacts'].items(): bind(bindings, root / name, digest)
        for path, digest in config['pinned_files'].items(): bind(bindings, path, digest)
        bind(bindings, 'scripts/audit_species_pmsf.py', a['script_sha256'])
        with (audit / 'branch_support.tsv').open() as f: rows = list(csv.DictReader(f, delimiter='\t'))
        assert len(rows) == 1046 and {r['tree'] for r in rows} == {'ml', 'consensus'}
        sources[label] = dict(spec=spec, rows=rows, audit=a)
    manifest = list(csv.DictReader(Path(plan['manifest']).open(), delimiter='\t'))
    assert len(manifest) == len({r['taxon_id'] for r in manifest}) == plan['expected']['taxa']
    assert {r['study_role'] for r in manifest} == {'ingroup', 'outgroup'}
    assert sum(r['study_role'] == 'ingroup' for r in manifest) == 501
    assert sum(r['study_role'] == 'outgroup' for r in manifest) == 25
    bind(bindings, plan['manifest']); bind(bindings, plan['fourth_completion']); verify(bindings)
    return sources, manifest, bindings
