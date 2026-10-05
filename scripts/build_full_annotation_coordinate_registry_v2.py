#!/usr/bin/env python3
"""Index every source annotation and protein product without changing the gene baseline."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import sqlite3
import time
from urllib.parse import unquote

from Bio import SeqIO
from reference_measurement_union_sources import bind, verify


def attributes(text, mode):
    """Retain repeated keys, encoded commas and nonstandard bare GTF labels."""
    result = defaultdict(list)
    if mode == 'creolimax_gtf':
        fields = re.findall(r'(\w+)\s+"([^\"]*)"\s*;?', text)
        for key, value in fields:
            result[key].append(value)
        if not fields and text.strip() not in ('', '.'):
            result['_unparsed_bare_gtf_label'].append(text)
    else:
        for field in text.split(';'):
            if field in ('', '.'):
                continue
            key, value = field.split('=', 1)
            result[key].extend(unquote(v) for v in value.split(','))
    return dict(result)


def initialize(connection):
    connection.executescript('''
    PRAGMA journal_mode=DELETE;
    PRAGMA synchronous=FULL;
    PRAGMA temp_store=FILE;
    PRAGMA cache_size=-32768;
    CREATE TABLE features (
      row_id INTEGER PRIMARY KEY, line_number INTEGER NOT NULL UNIQUE,
      seqid TEXT NOT NULL, source TEXT NOT NULL, feature_type TEXT NOT NULL,
      start INTEGER NOT NULL, end INTEGER NOT NULL, score TEXT NOT NULL,
      strand TEXT NOT NULL, phase TEXT NOT NULL, feature_id TEXT,
      attributes_json TEXT NOT NULL, raw_line TEXT NOT NULL);
    CREATE TABLE nonfeature_lines (line_number INTEGER PRIMARY KEY, raw_line TEXT NOT NULL);
    CREATE TABLE feature_parents (feature_row_id INTEGER NOT NULL, parent_id TEXT NOT NULL,
      FOREIGN KEY(feature_row_id) REFERENCES features(row_id));
    CREATE TABLE cds_product_refs (feature_row_id INTEGER NOT NULL, product_id TEXT NOT NULL,
      link_basis TEXT NOT NULL, in_source_proteome INTEGER NOT NULL,
      FOREIGN KEY(feature_row_id) REFERENCES features(row_id));
    CREATE TABLE products (protein_id TEXT PRIMARY KEY, source_ordinal INTEGER NOT NULL UNIQUE,
      description TEXT NOT NULL, protein_length INTEGER NOT NULL, sequence_sha256 TEXT NOT NULL,
      mapping_json TEXT NOT NULL, decision_json TEXT NOT NULL,
      selected_representative INTEGER NOT NULL);
    CREATE TABLE provisional_orfs (source_ordinal INTEGER PRIMARY KEY, protein_id TEXT NOT NULL,
      contig TEXT NOT NULL, start INTEGER NOT NULL, end INTEGER NOT NULL, strand TEXT NOT NULL,
      status TEXT NOT NULL, source_row_json TEXT NOT NULL);
    ''')


def load_products(entry, connection):
    selected = {}
    for record in SeqIO.parse(entry['representative_path'], 'fasta'):
        assert record.id not in selected
        selected[record.id] = hashlib.sha256(str(record.seq).encode('ascii')).hexdigest()
    count = 0
    decisions = Counter()
    identifiers = set()
    with open(entry['mapping_path']) as map_handle, open(entry['decisions_path']) as decision_handle:
        maps = iter(csv.DictReader(map_handle, delimiter='\t'))
        decs = iter(csv.DictReader(decision_handle, delimiter='\t'))
        for record in SeqIO.parse(entry['proteome_path'], 'fasta'):
            original, decision = next(maps), next(decs)
            assert all(decision[k] == v for k, v in original.items())
            assert original['taxon_id'] == entry['taxon_id']
            assert original['protein_id'] == record.id and int(original['protein_length']) == len(record.seq)
            assert record.id not in identifiers
            sequence_hash = hashlib.sha256(str(record.seq).encode('ascii')).hexdigest()
            included = record.id in selected
            assert included == (decision['decision'] != 'alternative_product_retained_in_source')
            if included:
                assert selected.pop(record.id) == sequence_hash
            count += 1
            identifiers.add(record.id)
            decisions[decision['decision']] += 1
            connection.execute('INSERT INTO products VALUES (?,?,?,?,?,?,?,?)',
                (record.id, count, record.description, len(record.seq), sequence_hash,
                 json.dumps(original, separators=(',', ':')), json.dumps(decision, separators=(',', ':')), int(included)))
        assert next(maps, None) is None and next(decs, None) is None and not selected
    assert count == entry['source_products']
    chosen = count - decisions['alternative_product_retained_in_source']
    assert chosen == entry['selected_representatives']
    return identifiers, dict(source_products=count, selected_representatives=chosen, decisions=dict(decisions))


def guard(root, budget):
    if shutil.disk_usage(root).free < budget['emergency_free_disk_gib'] * 2**30:
        raise RuntimeError('Annotation registry emergency disk reserve reached')
    total = sum(p.stat().st_size for p in root.glob('*.sqlite'))
    if total > budget['output_allowance_gib'] * 2**30:
        raise RuntimeError('Annotation registry aggregate output allowance reached')


def build_taxon(entry, output, budget):
    verify(entry['pins'])
    root = Path(output)
    path = root / (entry['taxon_id'] + '.sqlite')
    assert not path.exists()
    connection = sqlite3.connect(path)
    initialize(connection)
    protein_ids, product_counts = load_products(entry, connection)
    counts, flags, link_counts = Counter(), Counter(), Counter()
    row_id = 0
    if entry['mapping_mode'] == 'verified_orf_coordinates':
        with open(entry['annotation_path']) as handle:
            for ordinal, row in enumerate(csv.DictReader(handle, delimiter='\t'), 1):
                assert row['protein_id'] in protein_ids
                connection.execute('INSERT INTO provisional_orfs VALUES (?,?,?,?,?,?,?,?)',
                    (ordinal, row['protein_id'], row['contig'], int(row['start']), int(row['end']),
                     row['strand'], row['status'], json.dumps(row, separators=(',', ':'))))
                counts['provisional_orf_row'] += 1
                flags['provisional_orf:' + row['status']] += 1
        assert counts['provisional_orf_row'] == entry['source_products']
    else:
        opener = gzip.open if entry['annotation_path'].endswith('.gz') else open
        fasta_tail = False
        with opener(entry['annotation_path'], 'rt', newline='') as handle:
            for line_number, raw in enumerate(handle, 1):
                # Source bytes remain authoritative. Line endings alone are omitted in this index.
                line = raw.rstrip('\r\n')
                if line.startswith('##FASTA'):
                    fasta_tail = True
                if fasta_tail or line.startswith('#') or not line.strip():
                    connection.execute('INSERT INTO nonfeature_lines VALUES (?,?)', (line_number, line))
                    continue
                fields = line.split('\t')
                assert len(fields) == 9
                seqid, source, kind, start, end, score, strand, phase, text = fields
                start, end = int(start), int(end)
                assert start > 0 and end >= start
                attrs = attributes(text, entry['mapping_mode'])
                ids = attrs.get('ID', [])
                ident = ids[0] if len(ids) == 1 else None
                row_id += 1
                counts[kind] += 1
                for key in ('exception', 'transl_except', 'partial', 'start_range', 'end_range',
                            'Is_circular', 'part', 'transl_table', '_unparsed_bare_gtf_label'):
                    if key in attrs:
                        flags[kind + ':' + key] += 1
                if kind == 'CDS':
                    flags['CDS:phase:' + phase] += 1
                    if ident is None:
                        flags['CDS:nonunique_or_absent_ID'] += 1
                connection.execute('INSERT INTO features VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)',
                    (row_id, line_number, seqid, source, kind, start, end, score, strand, phase,
                     ident, json.dumps(attrs, separators=(',', ':')), line))
                for parent in attrs.get('Parent', []):
                    connection.execute('INSERT INTO feature_parents VALUES (?,?)', (row_id, parent))
                refs = set()
                if kind == 'CDS':
                    if entry['mapping_mode'] == 'ncbi_protein_gff':
                        refs.update((p, 'explicit_protein_id') for p in attrs.get('protein_id', []))
                    elif entry['mapping_mode'] == 'transcript_gff':
                        refs.update((p, 'direct_transcript_parent_candidate') for p in attrs.get('Parent', []))
                    elif entry['mapping_mode'] == 'creolimax_gtf':
                        refs.update((p, 'explicit_gtf_transcript_candidate') for p in attrs.get('transcript_id', []))
                        refs.update((p, 'gene_level_candidate_requires_review') for p in attrs.get('gene_id', []))
                    else:
                        raise ValueError('Unknown mapping mode')
                    if not refs:
                        flags['CDS:no_direct_product_reference'] += 1
                    for product, basis in sorted(refs):
                        connection.execute('INSERT INTO cds_product_refs VALUES (?,?,?,?)',
                            (row_id, product, basis, int(product in protein_ids)))
                        link_counts[basis + (':source_product' if product in protein_ids else ':outside_source_proteome')] += 1
                if row_id % 100000 == 0:
                    connection.commit()
                    guard(root, budget)
        if 'expected_feature_counts' in entry:
            assert dict(counts) == entry['expected_feature_counts'], entry['taxon_id']
    connection.commit()
    connection.executescript('''
    CREATE INDEX feature_type_index ON features(feature_type);
    CREATE INDEX feature_id_index ON features(feature_id);
    CREATE INDEX feature_coordinates_index ON features(seqid,start,end);
    CREATE INDEX parent_id_index ON feature_parents(parent_id);
    CREATE INDEX parent_feature_index ON feature_parents(feature_row_id);
    CREATE INDEX cds_product_index ON cds_product_refs(product_id);
    CREATE INDEX cds_feature_index ON cds_product_refs(feature_row_id);
    ''')
    assert connection.execute('PRAGMA integrity_check').fetchone() == ('ok',)
    assert not connection.execute('PRAGMA foreign_key_check').fetchall()
    product_without_cds = connection.execute('''SELECT count(*) FROM products p WHERE NOT EXISTS
      (SELECT 1 FROM cds_product_refs c WHERE c.product_id=p.protein_id)''').fetchone()[0]
    if entry['mapping_mode'] == 'verified_orf_coordinates':
        product_without_cds = None  # ORF coordinates are separately retained, never promoted to gene/CDS features.
    connection.close()
    guard(root, budget)
    verify(entry['pins'])
    pins = dict(entry['pins'])
    bind(pins, path)
    result = dict(taxon_id=entry['taxon_id'], mapping_mode=entry['mapping_mode'],
        annotation_path=entry['annotation_path'], database=str(path), database_bytes=path.stat().st_size,
        feature_rows=row_id, feature_counts=dict(counts), flags=dict(flags), candidate_link_counts=dict(link_counts),
        products_without_direct_cds_candidate=product_without_cds, **product_counts,
        source_hashes=pins, scientific_eligibility=False,
        scope='All original annotation feature rows, Parent references, source products and unchanged representative decisions. '
              'Candidate links do not establish CDS reconstruction, expression, gene copies, contamination or evolutionary events.')
    with (root / (entry['taxon_id'] + '.receipt.json')).open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    verify(plan['pins'])
    assert len(plan['entries']) == plan['expected_taxa'] == 526
    root = Path(plan['output'])
    root.mkdir(parents=True, exist_ok=False)
    assert not args.receipt.exists()
    start = time.monotonic()
    reports = []
    with ProcessPoolExecutor(max_workers=plan['cpu']) as pool:
        futures = [pool.submit(build_taxon, entry, str(root), plan['budget']) for entry in plan['entries']]
        for future in as_completed(futures):
            reports.append(future.result())
            state = dict(stage='indexing_full_source_annotations', completed_taxa=len(reports),
                expected_taxa=526, feature_rows=sum(r['feature_rows'] for r in reports),
                source_products=sum(r['source_products'] for r in reports),
                database_bytes=sum(r['database_bytes'] for r in reports), elapsed_seconds=time.monotonic()-start)
            temporary = root / 'state.partial'
            temporary.write_text(json.dumps(state, indent=2)+'\n')
            temporary.replace(root / 'state.json')
            print(json.dumps(state), flush=True)
    reports.sort(key=lambda r: r['taxon_id'])
    assert sum(r['source_products'] for r in reports) == plan['expected_source_products']
    assert sum(r['selected_representatives'] for r in reports) == plan['expected_selected_representatives']
    pins = dict(plan['pins'])
    for report in reports:
        for path, digest in report['source_hashes'].items():
            bind(pins, path, digest)
        bind(pins, root / (report['taxon_id'] + '.receipt.json'))
    bind(pins, args.plan)
    verify(pins)
    result = dict(status='complete_full_annotation_coordinate_registry_pending_independent_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(), taxa=526,
        source_products=plan['expected_source_products'], selected_representatives=plan['expected_selected_representatives'],
        feature_rows=sum(r['feature_rows'] for r in reports),
        database_bytes=sum(r['database_bytes'] for r in reports), taxa_reports=reports,
        elapsed_seconds=time.monotonic()-start, source_hashes=pins,
        scientific_eligibility=False, genome_cds_reconstruction_complete=False, gpu=False, new_predictions=0,
        scope='Full selected source annotation and product inventory, not independently replayed or genome-matched. '
              'All alternatives, unresolved mappings and provisional ORF exceptions retained; all eight evolutionary aims unfinished.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_hashes', 'taxa_reports')}), flush=True)


if __name__ == '__main__':
    main()
