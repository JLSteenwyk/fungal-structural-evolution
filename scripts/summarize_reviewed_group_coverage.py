#!/usr/bin/env python3
"""Describe source, taxon and marker coverage before/after representative-edge filtering."""
import argparse
import json
from collections import defaultdict
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha, read_table
from prepare_paired_phylogenetic_inputs import write_table


def counts(rows):
    before = {r['cluster_input_id'] for r in rows}
    after = {r['cluster_input_id'] for r in rows if r['reviewed_group_id']}
    passed = {r['cluster_input_id'] for r in rows if r['disposition'] == 'bidirectional_edge_pass'}
    representatives = {r['cluster_input_id'] for r in rows if r['disposition'] == 'representative_no_self_edge_test'}
    return {'input_models': len(before), 'retained_models': len(after),
            'deferred_models': len(before - after),
            'retained_representatives': len(representatives),
            'retained_nonself_members': len(passed),
            'retained_fraction': len(after) / len(before) if before else '',
            'nonself_retained_fraction': len(passed) / (len(before) - len(representatives)) if before - representatives else '',
            'input_markers': len({r['marker'] for r in rows}),
            'retained_markers': len({r['marker'] for r in rows if r['reviewed_group_id']}),
            'input_taxa': len({r['taxon_id'] for r in rows}),
            'retained_taxa': len({r['taxon_id'] for r in rows if r['reviewed_group_id']})}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['groups', 'manifest', 'output']:
        p.add_argument('--' + name, required=True, type=Path)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use immutable output')
    receipt = checked_receipt(a.groups)
    manifest = read_table(a.manifest)
    taxa = {r['taxon_id']: r for r in manifest}
    assert len(taxa) == len(manifest)
    links = read_table(a.groups / 'model_taxon_marker_links.tsv')
    models = read_table(a.groups / 'retained_membership.tsv') + read_table(a.groups / 'deferred_members.tsv')
    by_model = {r['cluster_input_id']: r for r in models}
    assert len(by_model) == len(models) == receipt['input_models']
    assert len(links) == receipt['model_taxon_marker_links']
    assert {r['cluster_input_id'] for r in links} == set(by_model)
    by_taxon, by_marker, by_lineage, by_source, by_stratum = [defaultdict(list) for _ in range(5)]
    for row in links:
        taxon = taxa[row['taxon_id']]
        model = by_model[row['cluster_input_id']]
        assert row['species_name'] == taxon['species_name'] and row['study_role'] == taxon['study_role']
        for field in ['source', 'model_id', 'reviewed_group_id', 'disposition']:
            assert row[field] == model[field]
        by_taxon[row['taxon_id']].append(row)
        by_marker[row['marker']].append(row)
        by_lineage[row['study_role'], row['major_lineage']].append(row)
        by_source[row['source']].append(row)
        length = int(model['length'])
        length_bin = '1-128' if length <= 128 else '129-256' if length <= 256 else '257-512' if length <= 512 else '513-1024' if length <= 1024 else '>1024'
        confidence = float(model['mean_ca_plddt'])
        confidence_bin = '<70' if confidence < 70 else '70-<90' if confidence < 90 else '>=90'
        by_stratum[row['source'], length_bin, confidence_bin].append(row)
    taxon_rows = [dict(taxon_id=t, species_name=r['species_name'], study_role=r['study_role'], lineage=r['lineage'], **counts(by_taxon.get(t, []))) for t, r in taxa.items()]
    marker_rows = [dict(marker=m, **counts(rows)) for m, rows in sorted(by_marker.items())]
    lineage_rows = [dict(study_role=role, major_lineage=lineage, **counts(rows)) for (role, lineage), rows in sorted(by_lineage.items())]
    source_rows = [dict(source=source, **counts(rows)) for source, rows in sorted(by_source.items())]
    stratum_rows = [dict(source=source, length_bin=length, mean_plddt_bin=confidence, **counts(rows)) for (source, length, confidence), rows in sorted(by_stratum.items())]
    a.output.mkdir(parents=True)
    for name, rows in [('taxon_coverage.tsv', taxon_rows), ('marker_coverage.tsv', marker_rows), ('lineage_coverage.tsv', lineage_rows), ('source_coverage.tsv', source_rows), ('source_length_confidence.tsv', stratum_rows)]:
        write_table(a.output / name, rows)
    result = {'status': 'complete_descriptive_reviewed_group_coverage',
              'group_receipt_sha256': sha(a.groups / 'receipt.json'), 'manifest_sha256': sha(a.manifest),
              'script_sha256': sha(Path(__file__)), 'manifest_taxa': len(taxa),
              'taxa_with_input_models': sum(r['input_models'] > 0 for r in taxon_rows),
              'taxa_with_retained_models': sum(r['retained_models'] > 0 for r in taxon_rows),
              'taxa_losing_all_available_models': [r['taxon_id'] for r in taxon_rows if r['input_models'] and not r['retained_models']],
              'source_summary': source_rows,
              'interpretation': 'Descriptive filtering coverage in frozen marker model collection; no phylogenetic inference. All manifest taxa retained in taxon table, including absent input coverage. Counts deduplicate models within each stratum; identical models linked to multiple taxa/markers mean taxon, marker and lineage totals are not additive. Representatives retained by definition are reported separately from tested nonself members. Mean pLDDT is not a residue-level mask; source/length/confidence strata are confounded and do not establish prediction accuracy or evolutionary effects.',
              'artifacts': {f.name: sha(f) for f in a.output.iterdir()}}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
