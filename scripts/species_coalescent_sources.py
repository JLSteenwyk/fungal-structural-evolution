"""Closed marker, support and taxon-cohort I/O for full ASTRAL sensitivities."""
import csv
from decimal import Decimal
import json
from pathlib import Path

from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha
from species_gcf_sources import load as load_markers
from retained_taxon_projection_sources import load as load_policies


def load(plan, path):
    marker_plan = json.loads(Path(plan['marker_source_plan']).read_text())
    _, universe, markers, bindings = load_markers(marker_plan, plan['marker_source_plan'])
    cohort_plan = json.loads(Path(plan['cohort_source_plan']).read_text())
    _, other_universe, policies, prior = load_policies(cohort_plan, plan['cohort_source_plan'])
    assert universe == other_universe
    for p, d in prior.items(): bind(bindings, p, d)
    for p, d in plan['pins'].items(): bind(bindings, p, d)
    bind(bindings, path)
    closure = json.loads(Path(plan['marker_completion']).read_text())
    assert closure['status'] == 'complete_verified_full_native_species_gene_concordance'
    assert closure['markers_per_alignment'] == 125 and closure['exact_process_journals_checked'] == 2
    bind(bindings, plan['marker_completion'])
    bind(bindings, closure['full_hash_archive'], closure['full_hash_archive_sha256'])
    archive = json.loads(Path(closure['full_hash_archive']).read_text())
    assert len(archive['services']) == 2
    for p, d in archive['source_hashes'].items(): bind(bindings, p, d)
    species = json.loads(Path(marker_plan['species_plan']).read_text())
    with Path(species['manifest']).open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    manifest = {row['taxon_id']: row for row in rows}
    assert len(manifest) == len(rows) == 526
    cohorts = dict(full_primary=dict(taxa=universe,
                   outgroups={r['taxon_id'] for r in rows if r['study_role'] == 'outgroup'},
                   roles=dict(ingroup=501, outgroup=25)), **policies)
    assert len(cohorts) == 5
    support = {}
    for alignment, items in markers.items():
        support_path = Path(marker_plan['markers'][alignment]['support']) / 'branch_support.tsv'
        bind(bindings, support_path)
        support[alignment] = {item['marker']: [] for item in items}
        with support_path.open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                assert row['marker'] in support[alignment]
                assert row['support_status'] == ('reported' if row['sh_alrt_percent'] else 'not_reported')
                side = frozenset(row['smaller_side_taxa'].split(';'))
                assert len(side) == int(row['smaller_side_size']) >= 2
                value = Decimal(row['sh_alrt_percent']) if row['sh_alrt_percent'] else None
                assert value is None or 0 <= value <= 100
                support[alignment][row['marker']].append(dict(side=side, support=value))
        assert all(len(support[alignment][item['marker']]) == item['taxa'] - 3 for item in items)
    assert plan['support_policies'] == {'uncontracted': None, 'sh_alrt10': 10, 'sh_alrt80': 80}
    verify(bindings)
    return universe, markers, cohorts, support, bindings
