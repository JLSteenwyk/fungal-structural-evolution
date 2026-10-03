#!/usr/bin/env python3
"""Audit every selected entry against frozen assembly and NCBI taxonomy sources.

This writes a new evidence overlay, never a replacement analysis manifest.
NCBI species-rank records and distinct accessions are not species delimitation.
"""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import tarfile


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(2**20), b''):
            h.update(chunk)
    return h.hexdigest()


def rows(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def dump_fields(line):
    parts = [p.strip() for p in line.decode().rstrip('\n').split('|')]
    assert parts[-1] == ''
    return parts[:-1]


def table(path, records):
    with path.open('x') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter='\t',
                                lineterminator='\n')
        writer.writeheader()
        writer.writerows(records)
    assert rows(path) == [{k: str(v) for k, v in r.items()} for r in records]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists() and not args.receipt.exists()
    manifest_path = Path('metadata/analysis_manifest.tsv')
    sources_path = Path('metadata/source_receipts.json')
    review_path = Path('metadata/taxon_label_review.tsv')
    hybrid_path = Path('config/curated_hybrid_evidence.json')
    override_path = Path('config/taxonomy_overrides.json')
    selected_path = Path('metadata/species_assembly_candidates.tsv')
    complex_path = Path('metadata/neocallimastix_identity_inputs_completed_20260927.json')
    cautions_path = Path('metadata/serendipita_oliveonia_taxonomic_review_20260927.json')
    pinned_paths = [Path(__file__), manifest_path, sources_path, review_path,
                   hybrid_path, override_path, selected_path, complex_path, cautions_path]
    source_receipts = json.loads(sources_path.read_text())
    archive = Path(source_receipts['ncbi_taxonomy']['path'])
    assert sha(archive) == source_receipts['ncbi_taxonomy']['sha256']
    pinned_paths.append(archive)
    manifest = rows(manifest_path)
    assert len(manifest) == 526
    assert Counter(r['study_role'] for r in manifest) == {'ingroup': 501, 'outgroup': 25}
    assert len({r['taxon_id'] for r in manifest}) == len(manifest)
    selected = {r['assembly_accession']: r for r in rows(selected_path)}
    flags = {r['taxon_id']: r for r in rows(review_path)}
    hybrids = {r['taxon_id']: r for r in json.loads(hybrid_path.read_text())['records']}
    overrides = json.loads(override_path.read_text())
    complex_ids = set(json.loads(complex_path.read_text())['group_taxa'])
    cautions = {r['taxon_id']: r for r in json.loads(cautions_path.read_text())['records']}
    for override in overrides.values():
        evidence = Path(override['evidence_path'])
        assert sha(evidence) == override['evidence_sha256']
        pinned_paths.append(evidence)
    wanted_assemblies = {r['assembly_accession'] for r in manifest}
    catalogue = defaultdict(list)
    for key, source in source_receipts.items():
        if not key.endswith('_fungi') and not key.endswith('_protozoa') and not key.endswith('_invertebrate'):
            continue
        path = Path(source['path'])
        assert sha(path) == source['sha256']
        pinned_paths.append(path)
        with path.open() as stream:
            header = None
            for line in stream:
                if line.startswith('#'):
                    candidate = line.lstrip('#').strip().split('\t')
                    if candidate[0] == 'assembly_accession':
                        header = candidate
                    continue
                assert header is not None
                fields = line.rstrip('\n').split('\t')
                assert len(fields) == len(header)
                if fields[0] in wanted_assemblies:
                    catalogue[fields[0]].append((key, dict(zip(header, fields))))
    bindings = {str(p): sha(p) for p in pinned_paths}
    nodes, merged, deleted = {}, {}, set()
    with tarfile.open(archive) as tar:
        with tar.extractfile('merged.dmp') as stream:
            for line in stream:
                old, new = dump_fields(line)
                assert old not in merged
                merged[old] = new
        with tar.extractfile('delnodes.dmp') as stream:
            deleted = {dump_fields(line)[0] for line in stream}
        with tar.extractfile('nodes.dmp') as stream:
            for line in stream:
                fields = dump_fields(line)
                taxid = fields[0]
                assert taxid not in nodes
                nodes[taxid] = (fields[1], fields[2])

        def canonical(original):
            current, seen = original, set()
            while current in merged:
                assert current not in seen, ('merged ID cycle', original)
                seen.add(current)
                current = merged[current]
            return current

        def lineage(original):
            current, path, seen = canonical(original), [], set()
            while current:
                if current not in nodes:
                    return path, 'deleted' if current in deleted else 'missing'
                assert current not in seen, ('node cycle', original)
                seen.add(current)
                parent, rank = nodes[current]
                path.append((current, rank))
                if parent == current:
                    assert current == '1', ('unexpected self-parent', current)
                    return path, 'resolved'
                current = parent
            return path, 'missing_identifier'

        prepared = []
        needed_names = set()
        for record in manifest:
            accession = record['assembly_accession']
            matches = catalogue.get(accession, [])
            assert len(matches) <= 1, ('ambiguous catalogue accession', accession)
            cat_source, cat = matches[0] if matches else ('', {})
            if record['study_role'] == 'ingroup':
                assert cat and accession in selected
                # Independent raw catalogue join must reproduce curated identifiers.
                for field in ['taxid', 'species_taxid', 'organism_name', 'infraspecific_name',
                              'biosample', 'bioproject', 'pubmed_id']:
                    assert cat[field] == selected[accession][field], (accession, field)
            assert cat or accession.startswith('figshare:'), ('uncatalogued assembly', accession)
            override = overrides.get(record['taxon_id'], {})
            effective_id = override.get('species_taxid', record['species_taxid'])
            tax_path, tax_status = lineage(effective_id)
            assembly_path, assembly_status = lineage(cat.get('taxid', ''))
            nearest = next((n for n, rank in tax_path if rank == 'species'), '')
            assembly_species = next((n for n, rank in assembly_path if rank == 'species'), '')
            needed_names.update(n for n, _ in tax_path + assembly_path)
            prepared.append((record, cat_source, cat, override, effective_id, tax_path,
                             tax_status, assembly_path, assembly_status, nearest, assembly_species))
        scientific, name_records, ranked = {}, defaultdict(list), {}
        with tar.extractfile('names.dmp') as stream:
            for line in stream:
                taxid, name, unique, kind = dump_fields(line)
                if taxid in needed_names:
                    name_records[taxid].append((name, kind))
                    if kind == 'scientific name':
                        assert taxid not in scientific
                        scientific[taxid] = name
        with tar.extractfile('rankedlineage.dmp') as stream:
            for line in stream:
                fields = dump_fields(line)
                if fields[0] in needed_names:
                    assert len(fields) == 10 and fields[0] not in ranked
                    ranked[fields[0]] = fields
    results, checks = [], []
    groups = defaultdict(list)
    for record, source, cat, override, effective, path, status, apath, astatus, nearest, anearest in prepared:
        resolved = canonical(effective)
        rank = nodes.get(resolved, ('', ''))[1]
        declared = canonical(cat.get('species_taxid', ''))
        expected_species_name = scientific.get(nearest, '')
        ranked_row = ranked.get(resolved)
        ranked_species_name = (ranked_row[1] if rank == 'species' else ranked_row[2]) if ranked_row else ''
        agrees = bool(nearest and ranked_row and expected_species_name == ranked_species_name)
        literal_name = scientific.get(resolved, '')
        aliases = name_records.get(resolved, [])
        raw_flag = flags.get(record['taxon_id'], {}).get('flags', '')
        risk = []
        if raw_flag == 'incompletely_identified_species_label': risk.append('uncertain_species_label')
        if record['taxon_id'] in hybrids: risk.append('curated_hybrid')
        if record['taxon_id'] in complex_ids: risk.append('neocallimastix_species_complex')
        if record['taxon_id'] in cautions: risk.append('strain_linked_taxonomic_caution')
        if rank != 'species': risk.append('effective_identifier_not_species_rank')
        if status != 'resolved': risk.append('unresolved_taxonomy_path')
        if not agrees: risk.append('rankedlineage_crosscheck_unresolved')
        if cat and anearest != declared: risk.append('assembly_species_ancestor_disagrees_with_catalogue')
        if cat and declared != resolved: risk.append('catalogue_species_identifier_disagrees_with_effective_identifier')
        if nearest: groups[(record['study_role'], nearest)].append(record['taxon_id'])
        results.append(dict(taxon_id=record['taxon_id'], study_role=record['study_role'],
            frozen_species_name=record['species_name'], frozen_species_taxid=record['species_taxid'],
            effective_species_taxid=effective, canonical_effective_taxid=resolved,
            effective_identifier_rank=rank, snapshot_scientific_name=literal_name,
            frozen_name_matches_scientific=str(record['species_name'] == literal_name).lower(),
            frozen_name_record_classes=';'.join(kind for name,kind in aliases if name == record['species_name']),
            taxonomy_path_status=status, nearest_species_rank_taxid=nearest,
            nearest_species_rank_name=expected_species_name,
            node_path_taxids=';'.join(n for n,_ in path),
            node_path_ranks=';'.join(rank for _,rank in path),
            rankedlineage_species_name=ranked_species_name,
            independent_rankedlineage_name_agreement=str(agrees).lower(),
            snapshot_kingdom=next((scientific[n] for n,r in path if r=='kingdom'), ''),
            snapshot_phylum=next((scientific[n] for n,r in path if r=='phylum'), ''),
            snapshot_genus=next((scientific[n] for n,r in path if r=='genus'), ''),
            assembly_accession=record['assembly_accession'], catalogue_source=source,
            catalogue_taxid=cat.get('taxid',''), catalogue_species_taxid=cat.get('species_taxid',''),
            catalogue_organism_name=cat.get('organism_name',''),
            catalogue_strain=cat.get('infraspecific_name',''), catalogue_isolate=cat.get('isolate',''),
            catalogue_biosample=cat.get('biosample',''), catalogue_bioproject=cat.get('bioproject',''),
            catalogue_pubmed_ids=cat.get('pubmed_id',''),
            assembly_taxonomy_path_status=astatus, assembly_nearest_species_taxid=anearest,
            assembly_type_representation=cat.get('assembly_type',''),
            catalogue_paired_accession=cat.get('gbrs_paired_asm',''),
            catalogue_paired_comparison=cat.get('paired_asm_comp',''),
            identifier_override_evidence=override.get('evidence_path',''),
            identity_review_flags=';'.join(risk), biological_species_delimitation='unverified'))
        checks.append(dict(taxon_id=record['taxon_id'], path_reaches_root=str(bool(path and path[-1][0]=='1')).lower(),
                           rankedlineage_name_agrees=str(agrees).lower(),
                           raw_catalogue_join='matched' if cat else 'not_applicable_figshare_deposition'))
    collisions = [dict(study_role=role, species_rank_taxid=taxid, species_rank_name=scientific[taxid],
                       taxon_ids=';'.join(ids), entry_count=len(ids))
                  for (role,taxid),ids in sorted(groups.items()) if len(ids)>1]
    collision_ids = {tid for r in collisions for tid in r['taxon_ids'].split(';')}
    for row in results:
        if row['taxon_id'] in collision_ids:
            row['identity_review_flags'] += (';' if row['identity_review_flags'] else '')+'shared_species_rank_identifier'
    exact_accessions = Counter(r['assembly_accession'] for r in results)
    paired_selected = [dict(taxon_id=r['taxon_id'], accession=r['assembly_accession'],
                           paired_accession=r['catalogue_paired_accession'],
                           paired_comparison=r['catalogue_paired_comparison']) for r in results
                       if r['catalogue_paired_accession'] in exact_accessions]
    args.output.mkdir(parents=True, exist_ok=False)
    table(args.output/'taxon_identity.tsv', results)
    table(args.output/'source_crosschecks.tsv', checks)
    collision_path = args.output/'species_rank_collisions.json'
    with collision_path.open('x') as stream:
        stream.write(json.dumps(dict(collisions=collisions, paired_selected_accessions=paired_selected,
            repeated_exact_accessions={k:v for k,v in exact_accessions.items() if v>1}),indent=2)+'\n')
    assert json.loads(collision_path.read_text())['collisions'] == collisions
    for p,h in bindings.items(): assert sha(p) == h, ('source changed during audit',p)
    report = dict(status='complete_selected_taxon_frozen_snapshot_identity_audit',
        checked_utc=datetime.now(timezone.utc).isoformat(), entries=len(results),
        role_counts=dict(Counter(r['study_role'] for r in results)),
        raw_catalogue_matches=sum(bool(r['catalogue_source']) for r in results),
        figshare_depositions=sum(r['assembly_accession'].startswith('figshare:') for r in results),
        identifier_overrides_used=[r['taxon_id'] for r in results if r['identifier_override_evidence']],
        rank_counts=dict(Counter(r['effective_identifier_rank'] for r in results)),
        resolved_paths=sum(r['taxonomy_path_status']=='resolved' for r in results),
        independent_rankedlineage_name_agreements=sum(r['independent_rankedlineage_name_agreement']=='true' for r in results),
        raw_catalogue_species_ancestor_disagreements=[r['taxon_id'] for r in results if 'assembly_species_ancestor_disagrees_with_catalogue' in r['identity_review_flags']],
        frozen_scientific_name_differences=[r['taxon_id'] for r in results if r['frozen_name_matches_scientific']=='false'],
        entry_review_flags=dict(Counter(flag for r in results for flag in r['identity_review_flags'].split(';') if flag)),
        species_rank_collisions=collisions,
        distinct_species_rank_identifiers_by_role={role:len({r['nearest_species_rank_taxid'] for r in results if r['study_role']==role and r['nearest_species_rank_taxid']}) for role in ['ingroup','outgroup']},
        biological_unique_species_count_verified=False, curated_hybrid_screen_comprehensive=False,
        source_hashes=bindings, artifacts={str(p):sha(p) for p in args.output.iterdir()},
        resource_plan=dict(cpus=1,address_space_gib=4,cpu_seconds=300,wall_seconds=600,
            output_allowance_mib=16,planning_runtime_seconds=[10,180],gpu=False,charges=False),
        sampling_inputs_changed=False, scientific_eligibility=False,
        scope='All526entries joined by exact accession to frozen raw catalogues or retained as five distinct Figshare depositions. Full nodes/merged/deleted traversal independently compared with names/rankedlineage. Prior verified identifier overlay retained separately. Species-rank identifiers are database records, not accepted species boundaries; uncertain names, curated hybrids and prior species-complex cautions remain explicit. No current-taxonomy claim, relabeling, pruning, new inference or live job mutation.')
    with args.receipt.open('x') as stream: stream.write(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['source_hashes','artifacts']},indent=2))


if __name__ == '__main__':
    main()
