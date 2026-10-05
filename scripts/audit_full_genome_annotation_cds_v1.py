#!/usr/bin/env python3
"""Compare every annotation coordinate and available deposited CDS with original assembly DNA."""
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

from Bio import SeqIO
from genomic_cds_join_v1 import coordinate_status, reconstruct
from reference_measurement_union_sources import bind, verify


def load_genome(entry):
    report = entry['genome_report']
    if report['status'] != 'verified_publisher_bound_genomic_dna':
        assert report['status'] == 'genome_acquisition_or_format_error'
        return {}, {}
    opener = gzip.open if report['path'].endswith('.gz') else open
    genomes, metadata = {}, {}
    with opener(report['path'], 'rt') as handle, opener(report['path'], 'rt') as header_handle:
        headers = (line[1:].rstrip('\r\n') for line in header_handle if line.startswith('>'))
        for record in SeqIO.parse(handle, 'fasta'):
            assert record.id not in genomes
            sequence = str(record.seq)
            genomes[record.id] = sequence
            header = next(headers)
            assert header.split()[0] == record.id
            metadata[record.id] = dict(sequence_id=record.id, description=header,
                length=len(sequence), sequence_sha256=hashlib.sha256(sequence.encode('ascii')).hexdigest(),
                uppercase_sequence_sha256=hashlib.sha256(sequence.upper().encode('ascii')).hexdigest())
        assert next(headers, None) is None
    with open(report['contig_metadata']) as handle:
        observed = {r['sequence_id']: r for r in map(json.loads, handle)}
    assert observed == metadata
    return genomes, {key: len(value) for key, value in genomes.items()}


def load_targets(entry):
    if entry['target_cds'] is None:
        return []
    source = entry['target_cds']
    path = source['path']
    opener = gzip.open if path.endswith('.gz') else open
    rows = []
    with opener(path, 'rt') as handle:
        for ordinal, record in enumerate(SeqIO.parse(handle, 'fasta'), 1):
            if source['identifier_mode'] == 'ncbi_protein_id_header':
                labels = re.findall(r'\[protein_id=([^\]]+)\]', record.description)
                product = labels[0] if len(labels) == 1 else None
            else:
                assert source['identifier_mode'] == 'exact_fasta_id'
                product = record.id
            sequence = str(record.seq).upper()
            rows.append(dict(ordinal=ordinal, cds_id=record.id, description=record.description,
                product_id=product, sequence=sequence, sequence_length=len(sequence),
                sequence_sha256=hashlib.sha256(sequence.encode('ascii')).hexdigest()))
    if 'expected_records' in source:
        assert len(rows) == source['expected_records']
    return rows


def candidate_groups(connection, mode):
    groups = defaultdict(lambda: defaultdict(list))
    if mode == 'verified_orf_coordinates':
        for ordinal, ident, seqid, start, end, strand, status, raw in connection.execute('SELECT * FROM provisional_orfs ORDER BY source_ordinal'):
            groups[ident][('provisional_orf', ordinal)].append(dict(row_id=ordinal, seqid=seqid,
                start=start, end=end, strand=strand, phase='unknown', attrs={},
                original_orf_status=status, link_bases=['original_provisional_orf_coordinate']))
        return groups
    query = '''SELECT c.product_id,c.link_basis,f.row_id,f.seqid,f.start,f.end,
      f.strand,f.phase,f.feature_id,f.attributes_json FROM cds_product_refs c
      JOIN features f ON f.row_id=c.feature_row_id ORDER BY c.product_id,f.row_id,c.link_basis'''
    for product, basis, row_id, seqid, start, end, strand, phase, feature_id, attributes in connection.execute(query):
        attrs = json.loads(attributes)
        if mode == 'ncbi_protein_gff':
            key = ('ncbi_feature_id', feature_id) if feature_id is not None else ('missing_unique_feature_id', row_id)
        elif mode == 'transcript_gff':
            key = ('published_transcript_parent_candidate', product)
        else:
            assert mode == 'creolimax_gtf'
            transcripts = attrs.get('transcript_id', [])
            key = ('gtf_transcript_candidate', transcripts[0]) if len(transcripts) == 1 else ('unresolved_gtf_transcript', row_id)
        existing = groups[product][key]
        if existing and existing[-1]['row_id'] == row_id:
            existing[-1]['link_bases'].append(basis)
        else:
            existing.append(dict(row_id=row_id, seqid=seqid, start=start, end=end, strand=strand,
                                 phase=phase, attrs=attrs, link_bases=[basis]))
    return groups


def compare_target(target, candidates, multiplicity):
    matches = [i for i, candidate in enumerate(candidates) if candidate['sequence'] == target['sequence']]
    if target['product_id'] is None:
        status = 'missing_or_ambiguous_target_product_id'
    elif not candidates:
        status = 'no_genomic_cds_candidate'
    elif not any(candidate['sequence'] is not None for candidate in candidates):
        status = 'all_genomic_candidates_require_review'
    elif not matches:
        status = 'no_exact_unmodified_genomic_cds_match'
    elif len(matches) > 1:
        status = 'multiple_exact_genomic_candidates_require_review'
    elif multiplicity != 1:
        status = 'exact_genomic_candidate_multiple_target_records_require_review'
    elif any('gene_level_candidate_requires_review' in p['link_bases'] for p in (candidates[matches[0]],)):
        status = 'exact_sequence_gene_level_candidate_requires_review'
    elif 'annotation_exception_retained_not_corrected' in candidates[matches[0]]['flags']:
        status = 'exact_genomic_candidate_with_annotation_exception_requires_review'
    else:
        status = 'one_exact_unmodified_genomic_cds_candidate'
    return dict(**{k:v for k,v in target.items() if k != 'sequence'}, status=status,
                candidate_count=len(candidates), matching_candidate_indices=matches,
                exact_sequence_match=bool(matches), candidates=[{k:v for k,v in c.items() if k != 'sequence'} for c in candidates])


def audit_taxon(entry, output, budget):
    verify(entry['pins'])
    root = Path(output) / entry['taxon_id']
    root.mkdir(exist_ok=False)
    connection = sqlite3.connect(Path(entry['registry_report']['database']).resolve().as_uri() + '?mode=ro', uri=True)
    genomes, lengths = load_genome(entry)
    circular = {seqid for seqid, attrs in connection.execute("SELECT seqid,attributes_json FROM features WHERE feature_type='region'")
                if json.loads(attrs).get('Is_circular') == ['true']}
    feature_counts, coordinate_counts = Counter(), Counter()
    coordinate_path = root / 'feature_coordinates.tsv.gz'
    with gzip.open(coordinate_path, 'wt') as handle:
        columns = ['source_kind', 'row_id', 'seqid', 'feature_type', 'start', 'end', 'strand',
                   'original_sequence_length', 'circular_annotation', 'coordinate_status']
        writer = csv.DictWriter(handle, columns, delimiter='\t', lineterminator='\n'); writer.writeheader()
        for row_id, seqid, kind, start, end, strand in connection.execute('SELECT row_id,seqid,feature_type,start,end,strand FROM features ORDER BY row_id'):
            status = coordinate_status(seqid, start, end, lengths, circular) if genomes else 'original_genome_unavailable'
            feature_counts[kind] += 1; coordinate_counts[kind + ':' + status] += 1
            writer.writerow(dict(source_kind='original_annotation_feature', row_id=row_id, seqid=seqid,
                feature_type=kind, start=start, end=end, strand=strand, original_sequence_length=lengths.get(seqid, ''),
                circular_annotation=seqid in circular, coordinate_status=status))
        for ordinal, ident, seqid, start, end, strand, status, raw in connection.execute('SELECT * FROM provisional_orfs ORDER BY source_ordinal'):
            observed = coordinate_status(seqid, start, end, lengths, set()) if genomes else 'original_genome_unavailable'
            coordinate_counts['provisional_orf:' + observed] += 1
            feature_counts['provisional_orf_row'] += 1
            writer.writerow(dict(source_kind='original_provisional_orf', row_id=ordinal, seqid=seqid,
                feature_type='provisional_orf', start=start, end=end, strand=strand,
                original_sequence_length=lengths.get(seqid, ''), circular_annotation=False, coordinate_status=observed))
    assert dict(feature_counts) == entry['registry_report']['feature_counts']
    groups = candidate_groups(connection, entry['mapping_mode'])
    candidates = {}
    candidate_counts = Counter()
    candidate_path = root / 'genomic_candidates.jsonl.gz'
    with gzip.open(candidate_path, 'wt') as handle:
        for product in sorted(groups):
            values = []
            for key, parts in sorted(groups[product].items(), key=lambda x: str(x[0])):
                if not genomes:
                    report = dict(status='original_genome_unavailable', flags=[], sequence_sha256=None,
                                  sequence_length=None, feature_row_ids=[p['row_id'] for p in parts])
                    sequence = None
                else:
                    report, sequence = reconstruct(parts, genomes, circular, lengths)
                if key[0] in ('missing_unique_feature_id', 'unresolved_gtf_transcript'):
                    report.update(status='unresolved_cds_feature_group_requires_review', sequence_sha256=None, sequence_length=None)
                    sequence = None
                if any(p.get('original_orf_status', 'exact_translation') != 'exact_translation' for p in parts):
                    report['flags'].append('original_provisional_orf_review_status_retained')
                report.update(group_key=list(key), link_bases=sorted({basis for p in parts for basis in p['link_bases']}))
                candidate_counts[report['status']] += 1
                handle.write(json.dumps(dict(product_id=product, **report), separators=(',', ':'))+'\n')
                values.append(dict(report, sequence=sequence))
            candidates[product] = values
    targets = load_targets(entry)
    multiplicity = Counter(row['product_id'] for row in targets)
    target_counts, linked = Counter(), defaultdict(list)
    target_path = root / 'cds_record_comparisons.jsonl.gz'
    with gzip.open(target_path, 'wt') as handle:
        for target in targets:
            result = compare_target(target, candidates.get(target['product_id'], []), multiplicity[target['product_id']])
            target_counts[result['status']] += 1
            linked[target['product_id']].append(dict(cds_id=target['cds_id'], ordinal=target['ordinal'], status=result['status']))
            handle.write(json.dumps(result, separators=(',', ':'))+'\n')
    product_path = root / 'source_product_dispositions.jsonl.gz'
    product_count = selected_count = 0
    with gzip.open(product_path, 'wt') as handle:
        for ident, ordinal, description, length, sequence_hash, mapping, decision, selected in connection.execute('SELECT * FROM products ORDER BY source_ordinal'):
            product_count += 1; selected_count += selected
            row = dict(protein_id=ident, source_ordinal=ordinal, sequence_sha256=sequence_hash,
                protein_length=length, original_mapping=json.loads(mapping), original_decision=json.loads(decision),
                selected_representative=bool(selected), genomic_candidate_count=len(candidates.get(ident, [])),
                original_cds_targets=linked.get(ident, []),
                missing_target_reason=entry['missing_target_reason'] if entry['target_cds'] is None else None)
            handle.write(json.dumps(row, separators=(',', ':'))+'\n')
    assert product_count == entry['registry_report']['source_products']
    assert selected_count == entry['registry_report']['selected_representatives']
    connection.close()
    if shutil.disk_usage(root).free < budget['emergency_free_disk_gib'] * 2**30:
        raise RuntimeError('Genomic CDS audit disk reserve reached')
    pins = dict(entry['pins'])
    for path in (coordinate_path, candidate_path, target_path, product_path):
        bind(pins, path)
    verify(pins)
    result = dict(taxon_id=entry['taxon_id'], status='all_source_coordinates_and_available_cds_compared_pending_independent_readback',
        mapping_mode=entry['mapping_mode'], source_products=product_count, selected_representatives=selected_count,
        feature_counts=dict(feature_counts), coordinate_counts=dict(coordinate_counts),
        candidate_counts=dict(candidate_counts), target_records=len(targets), target_status_counts=dict(target_counts),
        target_cds_source=entry['target_cds'], missing_target_reason=entry['missing_target_reason'],
        artifact_paths=[str(p) for p in (coordinate_path, candidate_path, target_path, product_path)],
        source_hashes=pins, scientific_eligibility=False,
        scope='Every source feature/ORF and protein retained; available original CDS records compared with unmodified genomic candidates. '
              'No phase trimming, recoding, frameshift repair, gene-copy or selection admission. Creolimax transcript sequence is not a CDS target.')
    with (root / 'receipt.json').open('x') as handle:
        json.dump(result, handle, indent=2); handle.write('\n')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text()); verify(plan['pins'])
    assert plan['expected_taxa'] == len(plan['entries']) == 526
    for name, expected in [('assembly_readback', 'passed_full_selected_genome_disposition_and_original_fasta_reconstruction'),
                           ('annotation_readback', 'passed_full_independent_annotation_coordinate_registry_replay')]:
        receipt = json.loads(Path(plan[name]).read_text())
        assert receipt['status'] == expected and receipt['taxa'] == 526
        verify(receipt['source_hashes'])
        transport = json.loads(Path(plan[name+'_transport']).read_text())
        assert transport['status'] == 'verified_original_software_wait_and_whole_wrapper_payloads'
        assert transport['original_tool_terminal_exit_code'] == 0
        verify(transport['source_hashes'])
    root = Path(plan['output']); root.mkdir(parents=True, exist_ok=False)
    assert not args.receipt.exists()
    start = time.monotonic(); reports = []
    with ProcessPoolExecutor(max_workers=plan['cpu']) as pool:
        futures = [pool.submit(audit_taxon, entry, str(root), plan['budget']) for entry in plan['entries']]
        for future in as_completed(futures):
            reports.append(future.result())
            state = dict(stage='full_genome_annotation_cds_comparison', completed_taxa=len(reports), expected_taxa=526,
                source_products=sum(r['source_products'] for r in reports), target_records=sum(r['target_records'] for r in reports),
                elapsed_seconds=time.monotonic()-start)
            (root / 'state.json').write_text(json.dumps(state, indent=2)+'\n')
            print(json.dumps(state), flush=True)
    assert sum(r['source_products'] for r in reports) == 5927745
    assert sum(r['selected_representatives'] for r in reports) == 5815847
    pins = dict(plan['pins'])
    for report in reports:
        for path, digest in report['source_hashes'].items():
            bind(pins, path, digest)
        bind(pins, root / report['taxon_id'] / 'receipt.json')
    bind(pins, args.plan); verify(pins)
    result = dict(status='complete_full_genome_annotation_cds_comparison_pending_independent_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(), taxa=526, source_products=5927745,
        selected_representatives=5815847, taxa_reports=sorted(reports, key=lambda r:r['taxon_id']),
        source_hashes=pins, elapsed_seconds=time.monotonic()-start, scientific_eligibility=False,
        gpu=False, new_predictions=0, all_eight_aims_incomplete=True,
        scope='All selected sources and products, every coordinate and available CDS target; errors/ambiguities/alternatives retained. '
              'Independent full genomic/sequence/order/disposition replay remains required. Not contamination, haplotig or evolutionary acceptance.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2); handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_hashes','taxa_reports')}, indent=2))


if __name__ == '__main__':
    main()
