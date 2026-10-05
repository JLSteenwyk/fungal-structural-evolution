#!/usr/bin/env python3
"""Replay every original annotation line and protein against the full coordinate index."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import sqlite3
import time
from urllib.parse import unquote_to_bytes

from Bio.SeqIO.FastaIO import SimpleFastaParser
from reference_measurement_union_sources import bind, verify


def decode(text, mode):
    result = {}
    if mode == 'creolimax_gtf':
        pairs = re.findall(r'([A-Za-z0-9_]+)\s+"(.*?)"', text)
        for key, value in pairs:
            result.setdefault(key, []).append(value)
        if not pairs and text.strip() not in ('', '.'):
            result['_unparsed_bare_gtf_label'] = [text]
    else:
        for part in text.split(';'):
            if not part or part == '.':
                continue
            key, separator, values = part.partition('=')
            assert separator
            for value in values.split(','):
                result.setdefault(key, []).append(unquote_to_bytes(value).decode('utf-8'))
    return result


def candidate_references(kind, attrs, mode):
    if kind != 'CDS':
        return []
    rules = {
        'ncbi_protein_gff': [('protein_id', 'explicit_protein_id')],
        'transcript_gff': [('Parent', 'direct_transcript_parent_candidate')],
        'creolimax_gtf': [('transcript_id', 'explicit_gtf_transcript_candidate'),
                         ('gene_id', 'gene_level_candidate_requires_review')]}
    return sorted({(value, basis) for key, basis in rules[mode] for value in attrs.get(key, [])})


def check_taxon(entry, report):
    verify(entry['pins'])
    assert report['taxon_id'] == entry['taxon_id'] and report['mapping_mode'] == entry['mapping_mode']
    verify(report['source_hashes'])
    path = Path(report['database'])
    assert path.stat().st_size == report['database_bytes']
    connection = sqlite3.connect(path.resolve().as_uri() + '?mode=ro', uri=True)
    assert connection.execute('PRAGMA integrity_check').fetchone() == ('ok',)
    assert not connection.execute('PRAGMA foreign_key_check').fetchall()
    records = iter(connection.execute('SELECT * FROM products ORDER BY source_ordinal'))
    selected = dict(connection.execute('SELECT protein_id,sequence_sha256 FROM products WHERE selected_representative=1'))
    identifiers = set()
    decisions = Counter()
    with open(entry['mapping_path']) as map_handle, open(entry['decisions_path']) as decision_handle, open(entry['proteome_path']) as fasta:
        maps = iter(csv.DictReader(map_handle, delimiter='\t'))
        decs = iter(csv.DictReader(decision_handle, delimiter='\t'))
        for ordinal, (title, sequence) in enumerate(SimpleFastaParser(fasta), 1):
            ident = title.split()[0]
            assert ident not in identifiers
            identifiers.add(ident)
            mapping, decision, saved = next(maps), next(decs), next(records)
            assert mapping['protein_id'] == ident and mapping['taxon_id'] == entry['taxon_id']
            assert int(mapping['protein_length']) == len(sequence)
            assert all(decision[key] == value for key, value in mapping.items())
            included = decision['decision'] != 'alternative_product_retained_in_source'
            digest = hashlib.sha256(sequence.encode('ascii')).hexdigest()
            assert saved[:5] == (ident, ordinal, title, len(sequence), digest)
            assert json.loads(saved[5]) == mapping and json.loads(saved[6]) == decision and saved[7] == int(included)
            decisions[decision['decision']] += 1
        assert next(records, None) is None and next(maps, None) is None and next(decs, None) is None
    with open(entry['representative_path']) as handle:
        for title, sequence in SimpleFastaParser(handle):
            assert selected.pop(title.split()[0]) == hashlib.sha256(sequence.encode('ascii')).hexdigest()
    assert not selected and len(identifiers) == entry['source_products'] == report['source_products']
    chosen = len(identifiers) - decisions['alternative_product_retained_in_source']
    assert chosen == report['selected_representatives'] == entry['selected_representatives']
    assert dict(decisions) == report['decisions']
    counts, flags, link_counts = Counter(), Counter(), Counter()
    features = iter(connection.execute('SELECT * FROM features ORDER BY row_id'))
    other_lines = iter(connection.execute('SELECT * FROM nonfeature_lines ORDER BY line_number'))
    parents = iter(connection.execute('SELECT feature_row_id,parent_id FROM feature_parents ORDER BY rowid'))
    references = iter(connection.execute('SELECT * FROM cds_product_refs ORDER BY rowid'))
    provisional = iter(connection.execute('SELECT * FROM provisional_orfs ORDER BY source_ordinal'))
    number = 0
    if entry['mapping_mode'] == 'verified_orf_coordinates':
        with open(entry['annotation_path']) as handle:
            for ordinal, row in enumerate(csv.DictReader(handle, delimiter='\t'), 1):
                saved = next(provisional)
                assert saved[:7] == (ordinal, row['protein_id'], row['contig'], int(row['start']),
                                     int(row['end']), row['strand'], row['status'])
                assert json.loads(saved[7]) == row and row['protein_id'] in identifiers
                counts['provisional_orf_row'] += 1
                flags['provisional_orf:' + row['status']] += 1
    else:
        opener = gzip.open if entry['annotation_path'].endswith('.gz') else open
        tail = False
        with opener(entry['annotation_path'], 'rt', newline='') as handle:
            for lineno, raw in enumerate(handle, 1):
                line = raw.rstrip('\r\n')
                tail = tail or line.startswith('##FASTA')
                if tail or line.startswith('#') or not line.strip():
                    assert next(other_lines) == (lineno, line)
                    continue
                fields = line.split('\t')
                assert len(fields) == 9
                attrs = decode(fields[-1], entry['mapping_mode'])
                number += 1
                saved = next(features)
                ident = attrs['ID'][0] if len(attrs.get('ID', [])) == 1 else None
                assert saved[:11] == (number, lineno, fields[0], fields[1], fields[2],
                    int(fields[3]), int(fields[4]), fields[5], fields[6], fields[7], ident)
                assert json.loads(saved[11]) == attrs and saved[12] == line
                for parent in attrs.get('Parent', []):
                    assert next(parents) == (number, parent)
                kind = fields[2]
                counts[kind] += 1
                for key in ('exception', 'transl_except', 'partial', 'start_range', 'end_range',
                            'Is_circular', 'part', 'transl_table', '_unparsed_bare_gtf_label'):
                    if key in attrs:
                        flags[kind + ':' + key] += 1
                if kind == 'CDS':
                    flags['CDS:phase:' + fields[7]] += 1
                    if ident is None:
                        flags['CDS:nonunique_or_absent_ID'] += 1
                    refs = candidate_references(kind, attrs, entry['mapping_mode'])
                    if not refs:
                        flags['CDS:no_direct_product_reference'] += 1
                    for product, basis in refs:
                        included = product in identifiers
                        assert next(references) == (number, product, basis, int(included))
                        link_counts[basis + (':source_product' if included else ':outside_source_proteome')] += 1
        if 'expected_feature_counts' in entry:
            assert dict(counts) == entry['expected_feature_counts']
    for cursor in (features, other_lines, parents, references, provisional):
        assert next(cursor, None) is None
    assert number == report['feature_rows'] and dict(counts) == report['feature_counts']
    assert dict(flags) == report['flags'] and dict(link_counts) == report['candidate_link_counts']
    without = None if entry['mapping_mode'] == 'verified_orf_coordinates' else connection.execute('''
      SELECT count(*) FROM products p WHERE NOT EXISTS
      (SELECT 1 FROM cds_product_refs c WHERE c.product_id=p.protein_id)''').fetchone()[0]
    assert without == report['products_without_direct_cds_candidate']
    connection.close()
    verify(entry['pins'])
    return dict(taxon_id=entry['taxon_id'], source_products=len(identifiers),
                selected_representatives=chosen, feature_rows=number, flags=dict(flags))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    verify(plan['pins'])
    original = json.loads(Path(plan['producer_plan']).read_text())
    receipt = json.loads(Path(plan['producer_receipt']).read_text())
    transport = json.loads(Path(plan['producer_transport']).read_text())
    assert receipt['status'] == 'complete_full_annotation_coordinate_registry_pending_independent_readback'
    assert transport['status'] == 'verified_original_software_wait_and_whole_wrapper_payloads'
    assert transport['original_tool_terminal_exit_code'] == 0
    verify(transport['source_hashes'])
    reports = {r['taxon_id']: r for r in receipt['taxa_reports']}
    assert len(reports) == len(original['entries']) == 526
    output = Path(plan['output'])
    output.mkdir(parents=True, exist_ok=False)
    assert not args.receipt.exists()
    start = time.monotonic()
    checks = []
    with ProcessPoolExecutor(max_workers=plan['cpu']) as pool:
        futures = [pool.submit(check_taxon, entry, reports[entry['taxon_id']]) for entry in original['entries']]
        for future in as_completed(futures):
            checks.append(future.result())
            state = dict(stage='independently_replaying_full_annotations', checked_taxa=len(checks), expected_taxa=526,
                         feature_rows=sum(r['feature_rows'] for r in checks), elapsed_seconds=time.monotonic()-start)
            (output / 'state.json').write_text(json.dumps(state, indent=2)+'\n')
            print(json.dumps(state), flush=True)
    checks.sort(key=lambda r: r['taxon_id'])
    pins = dict(plan['pins'])
    for path, digest in receipt['source_hashes'].items():
        bind(pins, path, digest)
    for path in (args.plan, plan['producer_receipt'], plan['producer_transport']):
        bind(pins, path)
    verify(pins)
    result = dict(status='passed_full_independent_annotation_coordinate_registry_replay',
        checked_utc=datetime.now(timezone.utc).isoformat(), taxa=526, taxa_checks=checks,
        feature_rows=sum(r['feature_rows'] for r in checks), source_products=sum(r['source_products'] for r in checks),
        selected_representatives=sum(r['selected_representatives'] for r in checks),
        elapsed_seconds=time.monotonic()-start, source_hashes=pins, scientific_eligibility=False,
        genome_cds_reconstruction_complete=False, gpu=False, new_predictions=0,
        scope='Every indexed original feature/nonfeature line, parent/candidate reference and source/representative protein '
              'replayed with separate decoding and FASTA parsing; shared Python/SQLite/Biopython dependencies remain. '
              'Does not establish genome-CDS agreement, annotation correctness, gene duplication or evolutionary events.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_hashes', 'taxa_checks')}), flush=True)


if __name__ == '__main__':
    main()
