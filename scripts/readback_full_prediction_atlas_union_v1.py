#!/usr/bin/env python3
"""Independently reconstruct every full atlas source and representative link."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sqlite3

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def fasta_records(path):
    title = None; parts = []
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith('>'):
                if title is not None: yield title.split()[0], ''.join(parts).replace(' ', '').replace('\r', '')
                title = line[1:].rstrip(); parts = []
            else:
                assert title is not None or not line.strip(), 'Sequence before FASTA header'
                if title is not None: parts.append(line.rstrip())
    if title is not None: yield title.split()[0], ''.join(parts).replace(' ', '').replace('\r', '')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args(); assert not args.receipt.exists()
    reader_plan = json.loads(args.plan.read_text()); pins = dict(reader_plan['pins']); verify(pins)
    plan_path = Path(reader_plan['producer_plan']); plan = json.loads(plan_path.read_text())
    producer_path = Path(reader_plan['producer_receipt']); transport_path = Path(reader_plan['producer_transport'])
    producer, transport = [json.loads(p.read_text()) for p in [producer_path, transport_path]]
    assert producer['status'] == 'completed_full_prediction_atlas_union_pending_independent_readback'
    assert transport['validation_sha256'] == sha(producer_path)
    assert transport['original_tool_terminal_exit_code'] == 0
    assert transport['whole_wrapper_initial_and_terminal_payloads_matched']
    assert transport['manager_start_records'] == transport['manager_completion_records'] == 1
    for mapping in [producer['source_hashes'], transport['source_hashes'], plan['pins']]:
        for path, digest in mapping.items(): bind(pins, path, digest)
    for path in [args.plan, plan_path, producer_path, transport_path, Path(__file__)]: bind(pins, path)
    verify(pins)
    root = Path(plan['output']); db = sqlite3.connect('file:'+str((root/'prediction_atlas.sqlite').resolve())+'?mode=ro', uri=True)
    afdb_lengths = {}; esm = {}; model_counts = Counter(); models_seen = set()
    with (root/'model_metadata.jsonl').open() as combined:
        emitted = (json.loads(line) for line in combined)
        for source, path in [('AFDB', Path(plan['afdb_catalog'])/'models.jsonl'), ('ESMFold', Path(plan['esmfold_inventory']))]:
            with path.open() as inp:
                for line in inp:
                    raw = json.loads(line); m = raw if source == 'AFDB' else raw['models'][0]
                    if source == 'ESMFold':
                        assert raw['status'] == 'verified' and len(raw['models']) == 1
                        assert raw['record_id'] == 'S'+m['sequence_sha256']
                    observed = next(emitted)
                    assert observed == dict(source=source, model=m)
                    identity = (source, m['sequence_sha256']); assert identity not in models_seen; models_seen.add(identity)
                    expected = (source,m['sequence_sha256'],m['model_id'],m['version'],m['path'],m['length'],
                        m['mean_ca_plddt'],m['fraction_ca_plddt_below50'],m['sha256'],m.get('prediction_config_sha256',''))
                    actual = db.execute('SELECT * FROM models WHERE source=? AND sequence_sha256=?',identity).fetchone()
                    assert actual == expected
                    if source == 'AFDB': afdb_lengths[m['sequence_sha256']] = m['length']
                    else: esm[m['sequence_sha256']] = m
                    model_counts[source] += 1
        assert next(emitted, None) is None
    assert model_counts == dict(AFDB=2935733, ESMFold=25322)
    assert db.execute('SELECT count(*) FROM models').fetchone() == (sum(model_counts.values()),)
    source_links = {}
    with (Path(plan['afdb_catalog'])/'protein_model_links.tsv').open() as inp:
        for row in csv.DictReader(inp, delimiter='\t'):
            key = (row['taxon_id'], row['protein_id']); assert key not in source_links
            expected = (row['sequence_sha256'], row['model_id'], int(row['version']), row['model_path'])
            assert db.execute('SELECT sequence_sha256,model_id,version,path FROM afdb_links WHERE taxon_id=? AND protein_id=?',key).fetchone() == expected
            source_links[key] = expected
    assert len(source_links) == producer['availability']['afdb_only']+producer['availability']['both'] == 2994868
    assert db.execute('SELECT count(*) FROM afdb_links').fetchone() == (len(source_links),)
    representatives = json.loads(Path(plan['representatives']).read_text())
    with Path(plan['manifest']).open() as inp: manifest = {r['taxon_id']: r for r in csv.DictReader(inp, delimiter='\t')}
    assert len(manifest) == 526
    coverage = {}; totals = Counter(); esm_used = set(); afdb_used = set()
    with (root/'protein_source_links.tsv').open() as inp:
        rows = iter(csv.DictReader(inp, delimiter='\t'))
        for entry in representatives['taxa']:
            taxon = entry['taxon_id']; local = Counter(); seen = set()
            assert sha(entry['path']) == entry['sha256']
            for protein, sequence in fasta_records(entry['path']):
                assert protein not in seen; seen.add(protein)
                digest = hashlib.sha256(sequence.encode()).hexdigest(); length = len(sequence)
                af = source_links.pop((taxon,protein),None); em = esm.get(digest)
                flags = (af is not None, em is not None)
                status = {(False,False):'neither',(True,False):'afdb_only',(False,True):'esmfold_only',(True,True):'both'}[flags]
                if af: assert af[0] == digest and afdb_lengths[digest] == length; afdb_used.add(digest)
                if em: assert em['length'] == length; esm_used.add(digest)
                expected = dict(taxon_id=taxon,protein_id=protein,sequence_sha256=digest,length=str(length),availability=status,
                    afdb_model_id=af[1] if af else '',afdb_version=str(af[2]) if af else '',afdb_path=af[3] if af else '',
                    esmfold_model_id=em['model_id'] if em else '',esmfold_version=str(em['version']) if em else '',
                    esmfold_path=em['path'] if em else '',esmfold_prediction_config_sha256=em['prediction_config_sha256'] if em else '')
                assert next(rows) == expected
                assert db.execute('SELECT * FROM proteins WHERE taxon_id=? AND protein_id=?',(taxon,protein)).fetchone() == (taxon,protein,digest,length,int(flags[0]),int(flags[1]))
                local.update(proteins=1, **{status:1}); totals.update(proteins=1, **{status:1})
            assert local['proteins'] == entry['selected_proteins']
            assert taxon not in coverage
            coverage[taxon] = dict(taxon_id=taxon,species_name=manifest[taxon]['species_name'],study_role=manifest[taxon]['study_role'],
                representative_proteins=str(local['proteins']),afdb_only=str(local['afdb_only']),esmfold_only=str(local['esmfold_only']),both=str(local['both']),neither=str(local['neither']))
        assert next(rows,None) is None and not source_links
    seen = set()
    with (root/'taxon_source_coverage.tsv').open() as inp:
        for row in csv.DictReader(inp,delimiter='\t'):
            taxon = row['taxon_id']; assert taxon not in seen; seen.add(taxon); assert row == coverage[taxon]
    assert seen == set(manifest) == set(coverage) and len(seen) == 526
    assert totals['proteins'] == 5815847
    availability = {k:totals[k] for k in ['afdb_only','esmfold_only','both','neither']}
    assert availability == producer['availability']
    assert len(afdb_used) == 2935733
    assert len(esm_used) == producer['esmfold_models_with_representative_links']
    assert 25322-len(esm_used) == producer['esmfold_models_without_representative_links']
    assert db.execute('SELECT count(*) FROM proteins').fetchone() == (5815847,)
    assert db.execute('PRAGMA integrity_check').fetchone() == ('ok',)
    assert not db.execute('PRAGMA foreign_key_check').fetchall()
    db.close(); verify(pins)
    result = dict(status='passed_full_prediction_atlas_union_independent_source_replay',checked_utc=datetime.now(timezone.utc).isoformat(),
        producer_receipt_sha256=sha(producer_path),taxa=526,representative_proteins=5815847,model_counts=dict(model_counts),
        availability=availability,esmfold_models_with_representative_links=len(esm_used),
        esmfold_models_without_representative_links=25322-len(esm_used),source_hashes=pins,
        scientific_eligibility=False,gpu=False,new_predictions=0,
        scope='Every original full source model and SQLite model field, every AFDB source link, every '
              'representative FASTA sequence/hash/length/source disposition and emitted TSV/SQLite row '
              'independently reconstructed without producer imports. All526taxon rows/zeros and source '
              'unlinked ESMFold models checked. Shares SQLite and SHA helpers, uses separate FASTA parser. '
              'Does not repeat producer coordinate/PAE byte hashes, confidence qualification or inference.')
    with args.receipt.open('x') as handle: json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__ == '__main__': main()
