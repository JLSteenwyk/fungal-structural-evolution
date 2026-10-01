"""Full closed source and retained-manifest I/O for baseline tree projections."""
import csv
import json
from pathlib import Path
from four_run_pmsf_sources import load_sources
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha

EDGE_FIELDS=['baseline_run','policy','tree_type','split_taxa_json','projected_branch_length_sum',
             'original_edge_components','empirical_projected_ufboot_percent','sh_alrt_percent']
BOOT_FIELDS=['baseline_run','policy','split_taxa_json','replicates_with_projected_split',
             'empirical_projected_ufboot_percent']
BOUNDARY_FIELDS=['baseline_run','policy','tree_type','retained_taxa','retained_ingroup','retained_outgroup',
                 'boundary_taxa_json','boundary_present','empirical_projected_ufboot_percent','sh_alrt_percent']


def load(plan, plan_path):
    bindings=dict(plan['pins']);bind(bindings,plan_path)
    sp=Path(plan['baseline_plan']);source_plan=json.loads(sp.read_text())
    closed=json.loads(Path(plan['baseline_completion']).read_text())
    assert closed['status']=='complete_verified_full_four_run_pmsf_ML_consensus_sensitivity' and len(closed['services'])==2
    assert closed['source_hashes'][str(sp)]==sha(sp)
    sources,manifest,source_bindings=load_sources(source_plan,sp)
    for record in [closed,dict(source_hashes=source_bindings)]:
        for p,d in record['source_hashes'].items():bind(bindings,p,d)
    bind(bindings,plan['baseline_completion'])
    ic=json.loads(Path(plan['input_completion']).read_text())
    assert ic['status']=='complete_verified_native_species_taxon_refit_inputs' and len(ic['services'])==2 and ic['summary']['matrices']==8
    for p,d in ic['source_hashes'].items():bind(bindings,p,d)
    bind(bindings,plan['input_completion'])
    ip=Path(plan['inputs'])/'receipt.json';bind(bindings,ip);inputs=json.loads(ip.read_text());assert len(inputs['matrices'])==8
    universe={r['taxon_id'] for r in manifest};manifest_by_taxon={r['taxon_id']:r for r in manifest};policies={}
    for spec in inputs['matrices']:
        assert spec['policy'] in ic['summary']['policies']
        root=Path(spec['path']);rp=root/'receipt.json';tp=root/'taxa.tsv'
        bind(bindings,rp,spec['matrix_receipt_sha256']);bind(bindings,root/'matrix.faa',spec['matrix_sha256']);bind(bindings,tp)
        with tp.open() as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
        retained={r['taxon_id'] for r in rows};assert len(rows)==len(retained)==spec['taxa'] and retained<universe
        assert all(r==manifest_by_taxon[r['taxon_id']] for r in rows)
        roles=ic['summary']['policy_roles'][spec['policy']];assert spec['roles']==roles
        assert all(sum(r['study_role']==role for r in rows)==n for role,n in roles.items())
        result=dict(taxa=retained,outgroups={r['taxon_id'] for r in rows if r['study_role']=='outgroup'},roles=roles)
        if spec['policy'] in policies:assert policies[spec['policy']]==result
        else:policies[spec['policy']]=result
    assert len(sources)==len(policies)==4 and len(universe)==526
    verify(bindings)
    return sources,universe,policies,bindings
