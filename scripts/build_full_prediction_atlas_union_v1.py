#!/usr/bin/env python3
"""Integrate full exact-sequence source availability without choosing predictors."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import time

from Bio.SeqIO.FastaIO import SimpleFastaParser

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args(); assert not args.receipt.exists()
    plan = json.loads(args.plan.read_text()); pins = dict(plan['pins']); verify(pins)
    root = Path(plan['output']); root.mkdir(exist_ok=False)
    started = time.monotonic(); counts = Counter(); model_counts = Counter()
    def state(stage):
        p = root / 'state.tmp'
        p.write_text(json.dumps(dict(stage=stage, counts=dict(counts), model_counts=dict(model_counts),
                                    elapsed_seconds=time.monotonic()-started))+'\n')
        p.replace(root / 'state.json')
    closure = json.loads(Path(plan['afdb_closure']).read_text())
    assert closure['status'] == 'complete_verified_full_completed_retrieval_catalog_refresh'
    assert closure['taxa'] == plan['expected_taxa'] == 526
    assert closure['proteins_screened'] == plan['expected_proteins'] == 5815847
    afdb_root = Path(plan['afdb_catalog']); raw = json.loads((afdb_root/'receipt.json').read_text())
    assert raw['proteins_linked'] == plan['expected_afdb_links'] == closure['proteins_linked']
    assert raw['unique_models'] == plan['expected_afdb_models'] == closure['unique_models']
    for name, digest in raw['artifacts'].items(): bind(pins, afdb_root/name, digest)
    er = json.loads(Path(plan['esmfold_readback']).read_text())
    assert er['status'] == 'passed_complete_disjoint_prediction_inventory_readback'
    assert er['models'] == plan['expected_esmfold_models'] == 25322
    assert er['producer_receipt_sha256'] == sha(Path(plan['esmfold_inventory']).parent/'receipt.json')
    for path, digest in er['input_hashes'].items(): bind(pins, path, digest)
    verify(pins)
    representatives = json.loads(Path(plan['representatives']).read_text())
    assert representatives['status'] == 'complete' and representatives['processed_taxa'] == 526
    with Path(plan['manifest']).open() as handle:
        manifest = {r['taxon_id']: r for r in csv.DictReader(handle, delimiter='\t')}
    assert len(manifest) == 526 and set(manifest) == {r['taxon_id'] for r in representatives['taxa']}
    assert Counter(r['study_role'] for r in manifest.values()) == dict(ingroup=501, outgroup=25)
    database = root/'prediction_atlas.sqlite'; db = sqlite3.connect(database)
    db.execute('PRAGMA foreign_keys=ON')
    db.executescript('''CREATE TABLE models(source TEXT,sequence_sha256 TEXT,model_id TEXT,version INTEGER,path TEXT,length INTEGER,mean_ca_plddt REAL,fraction_ca_plddt_below50 REAL,coordinate_sha256 TEXT,prediction_config_sha256 TEXT,PRIMARY KEY(source,sequence_sha256),UNIQUE(source,model_id,version),UNIQUE(source,path));
CREATE TABLE afdb_links(taxon_id TEXT,protein_id TEXT,sequence_sha256 TEXT,model_id TEXT,version INTEGER,path TEXT,PRIMARY KEY(taxon_id,protein_id));
CREATE TABLE proteins(taxon_id TEXT,protein_id TEXT,sequence_sha256 TEXT,length INTEGER,has_afdb INTEGER CHECK(has_afdb IN(0,1)),has_esmfold INTEGER CHECK(has_esmfold IN(0,1)),PRIMARY KEY(taxon_id,protein_id));''')
    esm = {}; esm_bytes = 0; esm_pae_bytes = 0; esm_receipts = 0
    state('importing_complete_source_models')
    with (root/'model_metadata.jsonl').open('x') as out:
        for source, path in [('AFDB', afdb_root/'models.jsonl'), ('ESMFold', Path(plan['esmfold_inventory']))]:
            with path.open() as handle:
                for line in handle:
                    item = json.loads(line)
                    if source == 'ESMFold':
                        assert item['status'] == 'verified' and len(item['models']) == 1
                        m = item['models'][0]
                        assert item['record_id'] == 'S'+m['sequence_sha256']
                        assert m['provider'] == 'local' and m['tool'] == 'ESMFold v1'
                        assert m['sequence_sha256'] not in esm
                        esm[m['sequence_sha256']] = m
                        assert sha(m['path']) == m['sha256']
                        assert sha(m['prediction_receipt_path']) == m['prediction_receipt_sha256']
                        assert sha(m['local_pae_npz_path']) == m['local_pae_npz_sha256']
                        esm_bytes += Path(m['path']).stat().st_size
                        esm_pae_bytes += Path(m['local_pae_npz_path']).stat().st_size; esm_receipts += 1
                    else:
                        m = item
                        assert m['provider'] == 'GDM' and m['tool'] == 'AlphaFold Monomer v2.0 pipeline'
                    assert isinstance(m['length'], int) and not isinstance(m['length'], bool) and m['length'] > 0
                    assert isinstance(m['version'], int) and not isinstance(m['version'], bool) and m['version'] > 0
                    for key in ['sequence_sha256', 'sha256']:
                        assert len(m[key]) == 64 and set(m[key]) <= set('0123456789abcdef')
                    for key, maximum in [('mean_ca_plddt', 100), ('fraction_ca_plddt_below50', 1)]:
                        value = m[key]
                        assert isinstance(value, (int, float)) and not isinstance(value, bool)
                        assert math.isfinite(value) and 0 <= value <= maximum
                    db.execute('INSERT INTO models VALUES (?,?,?,?,?,?,?,?,?,?)',
                        (source, m['sequence_sha256'], m['model_id'], m['version'], m['path'], m['length'],
                         m['mean_ca_plddt'], m['fraction_ca_plddt_below50'], m['sha256'], m.get('prediction_config_sha256', '')))
                    out.write(json.dumps(dict(source=source, model=m), separators=(',', ':'), allow_nan=False)+'\n')
                    model_counts[source] += 1
                    if model_counts[source] % 10000 == 0: db.commit(); state('importing_complete_source_models')
    assert model_counts == dict(AFDB=plan['expected_afdb_models'], ESMFold=25322)
    state('importing_all_afdb_protein_links')
    with (afdb_root/'protein_model_links.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            actual = db.execute("SELECT model_id,version,path FROM models WHERE source='AFDB' AND sequence_sha256=?", (row['sequence_sha256'],)).fetchone()
            assert actual == (row['model_id'], int(row['version']), row['model_path'])
            db.execute('INSERT INTO afdb_links VALUES (?,?,?,?,?,?)',
                (row['taxon_id'], row['protein_id'], row['sequence_sha256'], row['model_id'], int(row['version']), row['model_path']))
            counts['afdb_source_links'] += 1
            if counts['afdb_source_links'] % 10000 == 0: db.commit(); state('importing_all_afdb_protein_links')
    assert counts['afdb_source_links'] == plan['expected_afdb_links']
    db.commit(); summaries = []; esm_used = set(); afdb_used = set()
    fields = ['taxon_id','protein_id','sequence_sha256','length','availability','afdb_model_id','afdb_version','afdb_path','esmfold_model_id','esmfold_version','esmfold_path','esmfold_prediction_config_sha256']
    state('matching_every_representative_protein')
    with (root/'protein_source_links.tsv').open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n'); writer.writeheader()
        for entry in representatives['taxa']:
            fasta = Path(entry['path']); bind(pins, fasta, entry['sha256']); assert sha(fasta) == entry['sha256']
            taxon = manifest[entry['taxon_id']]; local = Counter(); seen = set()
            with fasta.open() as inp:
                for title, sequence in SimpleFastaParser(inp):
                    protein = title.split()[0]; assert protein not in seen; seen.add(protein)
                    digest = hashlib.sha256(sequence.encode()).hexdigest(); length = len(sequence)
                    af = db.execute('SELECT sequence_sha256,model_id,version,path FROM afdb_links WHERE taxon_id=? AND protein_id=?', (entry['taxon_id'], protein)).fetchone()
                    em = esm.get(digest)
                    if af:
                        assert af[0] == digest
                        assert db.execute("SELECT length FROM models WHERE source='AFDB' AND sequence_sha256=?", (digest,)).fetchone() == (length,)
                        afdb_used.add(digest)
                    if em: assert em['length'] == length; esm_used.add(digest)
                    availability = 'both' if af and em else 'afdb_only' if af else 'esmfold_only' if em else 'neither'
                    row = dict(taxon_id=entry['taxon_id'],protein_id=protein,sequence_sha256=digest,length=length,
                        availability=availability,afdb_model_id=af[1] if af else '',afdb_version=af[2] if af else '',
                        afdb_path=af[3] if af else '',esmfold_model_id=em['model_id'] if em else '',
                        esmfold_version=em['version'] if em else '',esmfold_path=em['path'] if em else '',
                        esmfold_prediction_config_sha256=em['prediction_config_sha256'] if em else '')
                    writer.writerow(row)
                    db.execute('INSERT INTO proteins VALUES (?,?,?,?,?,?)', (entry['taxon_id'],protein,digest,length,int(af is not None),int(em is not None)))
                    local.update(proteins=1, **{availability:1}); counts.update(proteins=1, **{availability:1})
                    if counts['proteins'] % 10000 == 0: db.commit(); state('matching_every_representative_protein')
            assert local['proteins'] == entry['selected_proteins']
            assert sum(local[k] for k in ['both','afdb_only','esmfold_only','neither']) == local['proteins']
            summaries.append(dict(taxon_id=entry['taxon_id'],species_name=taxon['species_name'],study_role=taxon['study_role'],
                representative_proteins=local['proteins'],afdb_only=local['afdb_only'],esmfold_only=local['esmfold_only'],
                both=local['both'],neither=local['neither']))
    db.commit(); state('checking_full_union_integrity')
    assert counts['proteins'] == 5815847 and len(summaries) == 526
    assert counts['afdb_only'] + counts['both'] == plan['expected_afdb_links']
    assert len(afdb_used) == plan['expected_afdb_models']
    assert db.execute('SELECT count(*) FROM proteins').fetchone() == (5815847,)
    assert db.execute('SELECT count(*) FROM models').fetchone() == (plan['expected_afdb_models']+25322,)
    assert db.execute('PRAGMA integrity_check').fetchone() == ('ok',)
    assert not db.execute('PRAGMA foreign_key_check').fetchall()
    db.close()
    with (root/'taxon_source_coverage.tsv').open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(summaries)
    verify(pins)
    for path in [args.plan, Path(__file__), *root.iterdir()]: bind(pins, path)
    result = dict(status='completed_full_prediction_atlas_union_pending_independent_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(),taxa=526,representative_proteins=5815847,
        model_counts=dict(model_counts),availability={k:counts[k] for k in ['afdb_only','esmfold_only','both','neither']},
        esmfold_models_with_representative_links=len(esm_used),esmfold_models_without_representative_links=25322-len(esm_used),
        esmfold_coordinate_bytes_rehashed=esm_bytes,esmfold_pae_bytes_rehashed=esm_pae_bytes,
        esmfold_prediction_receipts_rehashed=esm_receipts,afdb_coordinate_rehash_repeated=False,
        elapsed_seconds=time.monotonic()-started,source_hashes=pins,scientific_eligibility=False,gpu=False,new_predictions=0,
        scope='All 526taxa and5,815,847representatives matched to both full source inventories by exact sequence '
              'and length. Predictor identities and all original model metadata/configurations retained; '
              'both and neither explicit, no cross-predictor ranking. Every existing ESMFold coordinate/PAE '
              'and prediction receipt freshly hashed; AFDB coordinates rely on the bound completed full '
              'catalog hash run. Whole independent source/FASTA/link/model/coverage replay and original '
              'execution closure required. No residue-confidence qualification, inferred homology, '
              'prediction accuracy, missing-at-random assumption or completed evolutionary aim.')
    with args.receipt.open('x') as handle: json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__ == '__main__': main()
