#!/usr/bin/env python3
"""Join exact-sequence structures to marker alignments and retained matrix sites."""
import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from Bio import AlignIO, SeqIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PROVIDER = 'GDM'
DEFAULT_TOOL = 'AlphaFold Monomer v2.0 pipeline'

def inventory_key(row):
    key = row.get('record_id') or row.get('uniprot_accession')
    if not isinstance(key, str) or not key.strip():
        raise ValueError('Inventory record needs a local record ID or actual UniProt accession')
    return key


def source_candidates(candidates, provider, tool):
    return [m for m in candidates if m.get('provider') == provider and m.get('tool') == tool]



def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def residue_positions(aligned_sequence):
    number = 0
    positions = []
    for residue in aligned_sequence:
        if residue in '.-':
            positions.append(None)
        else:
            number += 1
            positions.append(number)
    return positions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--provider', default=DEFAULT_PROVIDER)
    parser.add_argument('--tool', default=DEFAULT_TOOL)
    parser.add_argument('--inventory', type=Path, default=ROOT / 'data/raw/afdb_models.jsonl')
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use an immutable new structure-mapping snapshot')
    inventory_raw = args.inventory.read_bytes()
    inventory_lines = inventory_raw.splitlines(keepends=True)
    if inventory_lines and not inventory_lines[-1].endswith(b'\n'):
        inventory_lines.pop()  # Exclude only an unfinished concurrent append.
    frozen_inventory = b''.join(inventory_lines)
    latest = {}
    for line in frozen_inventory.decode().splitlines():
        try:
            row = json.loads(line)
        except json.JSONDecodeError as error:
            raise ValueError('Malformed complete inventory record') from error
        latest[inventory_key(row)] = row
    models = defaultdict(list)
    for row in latest.values():
        if row['status'] == 'verified':
            for model in row['models']:
                origin = {'source_record_id': inventory_key(row)}
                if row.get('uniprot_accession'):
                    origin['uniprot_accession'] = row['uniprot_accession']
                models[model['sequence_sha256']].append(dict(model, **origin))
    with (ROOT / 'results/phylogeny/markers-full-v1/protein_mapping.tsv').open() as handle:
        markers = list(csv.DictReader(handle, delimiter='\t'))
    source = defaultdict(dict)
    for row in markers:
        source[row['marker']][row['taxon_id']] = row
    matrix = ROOT / 'results/phylogeny/profile-matrix-50-v1'
    matrix_receipt = json.loads((matrix / 'receipt.json').read_text())
    if digest(matrix / 'site_mapping.tsv') != matrix_receipt['artifacts']['site_mapping.tsv']:
        raise ValueError('Changed species-matrix site map')
    sites = defaultdict(dict)
    with (matrix / 'site_mapping.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            sites[row['marker']][int(row['alignment_column_1based'])] = int(row['matrix_column_1based'])
    args.output.mkdir(parents=True)
    (args.output / 'source_inventory.jsonl').write_bytes(frozen_inventory)
    links, summaries, verified_models, confidence_cache = [], [], {}, {}
    source_audit = []
    residue_path = args.output / 'matrix_to_structure_residues.tsv.gz'
    with gzip.open(residue_path, 'wt') as out:
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        writer.writerow(['marker', 'taxon_id', 'protein_id', 'model_id', 'model_version',
                         'matrix_column_1based', 'protein_residue_1based', 'ca_plddt'])
        mapped_residues = 0
        for marker, taxa in sorted(source.items()):
            raw = ROOT / f'results/phylogeny/markers-full-v1/unaligned/{marker}.faa'
            with raw.open() as handle:
                sequences = {r.id: str(r.seq) for r in SeqIO.parse(handle, 'fasta')}
            profile = ROOT / f'results/phylogeny/profile-alignments-full-v1/{marker}'
            receipt = json.loads(profile.with_suffix('.receipt.json').read_text())
            if digest(profile.with_suffix('.sto')) != receipt['stockholm_sha256']:
                raise ValueError('Changed full profile alignment')
            alignment = {r.id: str(r.seq) for r in AlignIO.read(profile.with_suffix('.sto'), 'stockholm')}
            linked = 0
            for taxon, row in sorted(taxa.items()):
                sequence = sequences[taxon]
                if hashlib.sha256(sequence.encode()).hexdigest() != row['sequence_sha256']:
                    raise ValueError('Marker mapping sequence differs from source')
                all_candidates = models.get(row['sequence_sha256'], [])
                candidates = source_candidates(all_candidates, args.provider, args.tool)
                source_audit.append({'marker': marker, 'taxon_id': taxon, 'sequence_sha256': row['sequence_sha256'], 'all_candidate_models': len(all_candidates), 'source_eligible_models': len(candidates), 'status': 'eligible_source_available' if candidates else ('only_other_sources_available' if all_candidates else 'no_model_in_snapshot')})
                if not candidates:
                    continue
                # Rank only within the explicitly selected pipeline; never compare confidence across pipelines.
                model = max(candidates, key=lambda m: (m['mean_ca_plddt'], m['version'], m['model_id']))
                key = model['path']
                if key not in confidence_cache:
                    path = ROOT / key
                    if digest(path) != model['sha256']:
                        raise ValueError('Changed structural coordinates')
                    cif = MMCIF2Dict(str(path))
                    polymer = [''.join(s.split()) for s in cif['_entity_poly.pdbx_seq_one_letter_code_can']]
                    if polymer != [sequence]:
                        raise ValueError('CIF sequence differs from marker')
                    ca = {int(position): float(value) for atom, position, value in zip(
                        cif['_atom_site.label_atom_id'], cif['_atom_site.label_seq_id'], cif['_atom_site.B_iso_or_equiv']) if atom == 'CA'}
                    if set(ca) != set(range(1, len(sequence) + 1)):
                        raise ValueError('CIF lacks complete marker residue coverage')
                    confidence_cache[key] = ca
                    verified_models[key] = model
                aligned = alignment[taxon]
                if aligned.upper().replace('-', '').replace('.', '') != sequence:
                    raise ValueError('Stockholm sequence differs from source')
                positions = residue_positions(aligned)
                observed = []
                for profile_column, matrix_column in sorted(sites[marker].items()):
                    stockholm_column = receipt['retained_stockholm_columns_1based'][profile_column - 1]
                    position = positions[stockholm_column - 1]
                    if position is None:
                        continue
                    confidence = confidence_cache[key][position]
                    observed.append(confidence)
                    writer.writerow([marker, taxon, row['protein_id'], model['model_id'], model['version'], matrix_column, position, confidence])
                    mapped_residues += 1
                links.append({'marker': marker, 'taxon_id': taxon, 'protein_id': row['protein_id'],
                              'sequence_sha256': row['sequence_sha256'], 'model_id': model['model_id'],
                              'model_version': model['version'], 'model_path': key, 'model_sha256': model['sha256'],
                              'model_provider': model.get('provider'), 'model_tool': model.get('tool'),
                              'all_source_candidate_models': len(all_candidates),
                              'exact_sequence_candidate_models': len(candidates), 'retained_marker_residues': len(observed),
                              'mean_retained_ca_plddt': sum(observed) / len(observed) if observed else '',
                              'fraction_retained_ca_plddt_ge70': sum(p >= 70 for p in observed) / len(observed) if observed else ''})
                linked += 1
            summaries.append({'marker': marker, 'single_copy_taxa': len(taxa), 'taxa_with_exact_sequence_structure': linked})
    for name, rows in [('marker_structure_links.tsv', links), ('marker_coverage.tsv', summaries), ('source_selection_audit.tsv', source_audit)]:
        with (args.output / name).open('w') as handle:
            if rows:
                writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
                writer.writeheader()
                writer.writerows(rows)
    (args.output / 'model_provenance.json').write_text(json.dumps(list(verified_models.values()), indent=2) + '\n')
    result = {'source_policy': {'provider': args.provider, 'tool': args.tool, 'selection': 'Highest mean CA pLDDT, version and model ID only within this exact provider/tool; no cross-pipeline fallback. Confidence ranking does not establish accuracy.'},
              'source_inventory_path': str(args.inventory),
              'source_inventory_sha256': hashlib.sha256(frozen_inventory).hexdigest(),
              'script_sha256': digest(Path(__file__)),
              'marker_proteins_screened': len(markers), 'marker_proteins_linked': len(links),
              'distinct_taxa_linked': len({r['taxon_id'] for r in links}), 'distinct_models': len(verified_models),
              'matrix_residue_links': mapped_residues, 'markers_with_at_least_4_linked_taxa': sum(r['taxa_with_exact_sequence_structure'] >= 4 for r in summaries),
              'matrix_receipt_sha256': digest(matrix / 'receipt.json'),
              'assignment_basis': 'Exact complete amino-acid sequence, including reuse across taxa with identical sequences. Taxon identities and source accessions are retained.',
              'caveat': 'Partial acquisition snapshot. Confidence is descriptive; branch inference, PAE/domain assessment and source-sensitivity checks remain pending.',
              'artifacts': {p.name: digest(p) for p in args.output.iterdir()}}
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
