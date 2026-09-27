#!/usr/bin/env python3
"""Check every case-homolog CA grid against raw mmCIF and independent metadata."""
import csv
import gzip
import json
import math
import subprocess
import time
from collections import Counter, defaultdict
from pathlib import Path

import psutil
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from screen_duplication_domain_alignment_coverage import sha

BASE = Path('results/experimental_structures')
MAPPING = BASE / 'whole-domain-case-ca-mapping-20260927-v1'
COORDINATES = BASE / 'whole-domain-case-coordinates-20260927-v1'
METADATA = BASE / 'whole-domain-case-metadata-20260927-v1'
OUTPUT = BASE / 'whole-domain-case-ca-readback-20260927-v1'
LAUNCH = Path('metadata/case_experimental_ca_mapping_launch_20260927.json')
AA = {k.upper(): v for k, v in protein_letters_3to1.items()}


def status_for(records, expected):
    if len(records) == 0:
        return 'no_CA_observation'
    if len(records) > 1:
        return 'multiple_CA_records'
    atom = records[0]
    if AA.get(atom['label_comp_id']) != expected:
        return 'nonstandard_or_mismatching_monomer'
    try:
        xyz = tuple(float(atom[name]) for name in ('Cartn_x', 'Cartn_y', 'Cartn_z'))
        occupancy = float(atom['occupancy'])
    except ValueError:
        return 'invalid_coordinates_or_occupancy'
    if not all(math.isfinite(x) for x in (*xyz, occupancy)) or not 0 < occupancy <= 1:
        return 'invalid_coordinates_or_occupancy'
    if occupancy < 1 or atom['label_alt_id'] not in {'.', '?'}:
        return 'alternate_or_partial_occupancy'
    return 'unambiguous_full_occupancy_CA'


def main():
    launch = json.loads(LAUNCH.read_text())
    sources = {str(LAUNCH): sha(LAUNCH)}
    while True:
        try:
            process = psutil.Process(launch['pid'])
            if process.create_time() != launch['created'] or process.status() == psutil.STATUS_ZOMBIE:
                break
            assert process.cmdline() == launch['cmdline']
        except psutil.NoSuchProcess:
            break
        print('Waiting for exact CA producer', launch['pid'], flush=True)
        time.sleep(30)
    state = dict(line.split('=', 1) for line in subprocess.check_output(
        ['systemctl', '--user', 'show', launch['unit'], '-p', 'ActiveState',
         '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
    assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0')
    assert sha(launch['cmdline'][1]) == launch['script_sha256']
    receipt = json.loads((MAPPING / 'receipt.json').read_text())
    config = json.loads((MAPPING / 'config.json').read_text())
    assert receipt['status'] == 'complete_experimental_CA_correspondence_screen'
    assert receipt['config_sha256'] == sha(MAPPING / 'config.json')
    assert config['metadata_receipt_sha256'] == sha(METADATA / 'receipt.json')
    assert config['coordinate_receipt_sha256'] == sha(COORDINATES / 'receipt.json')
    archive_audit = BASE / 'whole-domain-case-coordinate-readback-20260927-v1/receipt.json'
    assert config['archive_audit_receipt_sha256'] == sha(archive_audit)
    for path in [MAPPING / 'receipt.json', MAPPING / 'config.json', METADATA / 'receipt.json', COORDINATES / 'receipt.json', archive_audit]:
        sources[str(path)] = sha(path)
    archives = {x['entry_id']: x for x in json.loads((COORDINATES / 'receipt.json').read_text())['results']}
    targets = defaultdict(dict)
    for item in json.loads((METADATA / 'receipt.json').read_text())['responses']:
        if item['kind'] != 'polymer_entity':
            continue
        path = METADATA / 'polymer_entity' / (item['identifier'] + '.json')
        assert sha(path) == item['response_sha256']
        sources[str(path)] = sha(path)
        entry, entity = item['identifier'].rsplit('_', 1)
        targets[entry][entity] = ''.join(json.loads(path.read_text())['entity_poly']['pdbx_seq_one_letter_code_can'].split())
    assert set(targets) - set(config['deferred_entries']) == set(archives) == set(receipt['entry_receipt_sha256'])
    totals, summaries, raw_atoms = Counter(), [], 0
    fields = ['label_atom_id','label_entity_id','label_asym_id','label_seq_id','label_comp_id','label_alt_id','pdbx_PDB_model_num','Cartn_x','Cartn_y','Cartn_z','occupancy','B_iso_or_equiv','auth_asym_id','auth_seq_id','pdbx_PDB_ins_code']
    for entry, rh in receipt['entry_receipt_sha256'].items():
        rp = MAPPING / (entry + '.receipt.json')
        assert sha(rp) == rh
        er = json.loads(rp.read_text())
        table = MAPPING / (entry + '.residues.tsv.gz')
        assert er['config_sha256'] == sha(MAPPING / 'config.json')
        assert sha(table) == er['table_sha256']
        path = COORDINATES / (entry + '.cif.gz')
        assert sha(path) == archives[entry]['gzip_sha256']
        sources[str(rp)] = sha(rp)
        sources[str(table)] = sha(table)
        with gzip.open(path, 'rt') as handle:
            data = MMCIF2Dict(handle)
        sequences = dict(zip(data['_entity_poly.entity_id'], data['_entity_poly.pdbx_seq_one_letter_code_can']))
        for entity, seq in targets[entry].items():
            assert ''.join(sequences[entity].split()) == seq
        chain_entity = {chain: entity for chain, entity in zip(data['_struct_asym.id'], data['_struct_asym.entity_id']) if entity in targets[entry]}
        assert set(chain_entity.values()) == set(targets[entry])
        models = set(data['_atom_site.pdbx_PDB_model_num'])
        raw = defaultdict(list)
        cols = {field: data['_atom_site.' + field] for field in fields}
        assert len({len(x) for x in cols.values()}) == 1
        for i, atom_name in enumerate(cols['label_atom_id']):
            if atom_name != 'CA' or cols['label_entity_id'][i] not in targets[entry]:
                continue
            atom = {field: cols[field][i] for field in fields}
            assert chain_entity[atom['label_asym_id']] == atom['label_entity_id']
            pos = int(atom['label_seq_id'])
            assert 1 <= pos <= len(targets[entry][atom['label_entity_id']])
            raw[(atom['pdbx_PDB_model_num'], atom['label_asym_id'], pos)].append(atom)
        expected_grids = {(model, chain) for model in models for chain in chain_entity}
        seen = defaultdict(set)
        grid_counts = defaultdict(Counter)
        counts = Counter()
        with gzip.open(table, 'rt') as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                assert row['entry_id'] == entry
                model, chain, pos = row['model_number'], row['label_asym_id'], int(row['label_seq_id'])
                key = (model, chain)
                assert key in expected_grids and pos not in seen[key]
                seen[key].add(pos)
                entity = chain_entity[chain]
                assert row['entity_id'] == entity
                assert 1 <= pos <= len(targets[entry][entity])
                assert row['target_aa'] == targets[entry][entity][pos - 1]
                atoms = raw.pop((model, chain, pos), [])
                assert json.loads(row['atom_records_json']) == atoms
                assert int(row['CA_record_count']) == len(atoms)
                status = status_for(atoms, row['target_aa'])
                assert row['CA_status'] == status
                counts[status] += 1
                grid_counts[key][status] += 1
                raw_atoms += len(atoms)
        assert not raw and set(seen) == expected_grids
        for (model, chain), positions in seen.items():
            entity = chain_entity[chain]
            length = len(targets[entry][entity])
            assert positions == set(range(1, length + 1))
            summaries.append(dict(entry_id=entry, entity_id=entry+'_'+entity,
                                  model_number=model, label_asym_id=chain, canonical_length=length,
                                  unambiguous_CA=grid_counts[model, chain]['unambiguous_full_occupancy_CA'],
                                  status_counts_json=json.dumps(dict(grid_counts[model, chain]), sort_keys=True)))
        assert len(models) == er['deposited_models'] and len(chain_entity) == er['selected_chains']
        assert dict(counts) == er['status_counts'] and sum(counts.values()) == er['position_rows']
        totals.update(counts)
        print(entry, sum(counts.values()), 'raw atom and grid checks passed', flush=True)
    assert len(archives) == receipt['entries']
    assert dict(totals) == receipt['status_counts'] and sum(totals.values()) == receipt['position_rows']
    for path, digest in sources.items():
        assert sha(path) == digest
    OUTPUT.mkdir(parents=True, exist_ok=False)
    table = OUTPUT / 'chain_model_coverage.tsv'
    with table.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(summaries)
    with table.open() as handle:
        assert list(csv.DictReader(handle, delimiter='\t')) == [{k: str(v) for k, v in row.items()} for row in summaries]
    result = dict(status='complete_full_case_CA_raw_atom_and_grid_readback', terminal_state=state,
                  source_hashes=sources, checker_sha256=sha(__file__), entries=len(archives),
                  chain_models=len(summaries), positions=sum(totals.values()), CA_atoms=raw_atoms,
                  status_counts=dict(totals), artifacts={table.name: sha(table)},
                  scope='Every selected entity, chain, model, canonical sequence position, exported CA atom field and status checked against raw mmCIF and metadata. Producer and checker share the Biopython mmCIF parser; mapping and status logic are separate. Entity coverage is not fungal query/domain coverage. No geometry, experimental quality or training-independence validation.')
    (OUTPUT / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
