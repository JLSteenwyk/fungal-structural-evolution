#!/usr/bin/env python3
"""Independently reconstruct every genomic/CDS candidate and all full-source dispositions."""
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
import sqlite3
import time

from reference_measurement_union_sources import bind, verify


COMPLEMENT = bytes.maketrans(b'ACGTRYSWKMBDHVNacgtryswkmbdhvn',
                             b'TGCAYRSWMKVHDBNtgcayrswmkvhdbn')


def fasta_bytes(path):
    """Own streaming FASTA decoder; retain original header whitespace and sequence case."""
    opener = gzip.open if str(path).endswith('.gz') else open
    title = None
    chunks = []
    with opener(path, 'rb') as handle:
        for raw in handle:
            line = raw.rstrip(b'\r\n')
            if not line.strip():
                continue
            if line.startswith(b'>'):
                if title is not None:
                    assert chunks
                    yield title, b''.join(chunks)
                title = line[1:].decode('utf-8')
                assert title.split()
                chunks = []
            else:
                assert title is not None
                chunks.append(b''.join(line.split()))
        if title is not None:
            assert chunks
            yield title, b''.join(chunks)


def bounds(seqid, start, end, sizes, circular):
    size = sizes.get(seqid)
    if size is None:
        return 'missing_genomic_sequence_id'
    if not (0 < start <= end):
        return 'invalid_coordinate_interval'
    if end <= size:
        return 'within_original_sequence_bounds'
    if seqid not in circular:
        return 'out_of_bounds_without_circular_annotation'
    if end > size * 2 or end - start >= size:
        return 'circular_span_requires_review'
    return 'documented_circular_virtual_coordinates'


def genomic_sequences(entry):
    report = entry['genome_report']
    if report['status'] != 'verified_publisher_bound_genomic_dna':
        assert report['status'] == 'genome_acquisition_or_format_error'
        return {}, {}
    genome = {}
    with open(report['contig_metadata']) as handle:
        rows = iter(map(json.loads, handle))
        for title, sequence in fasta_bytes(report['path']):
            ident = title.split()[0]
            assert ident not in genome and sequence
            assert next(rows) == dict(sequence_id=ident, description=title, length=len(sequence),
                sequence_sha256=hashlib.sha256(sequence).hexdigest(),
                uppercase_sequence_sha256=hashlib.sha256(sequence.upper()).hexdigest())
            genome[ident] = sequence
        assert next(rows, None) is None
    assert genome
    return genome, {ident:len(sequence) for ident, sequence in genome.items()}


def expected_candidate(parts, genome, sizes, circular):
    checks = [bounds(p['seqid'], p['start'], p['end'], sizes, circular) for p in parts]
    row = dict(feature_row_ids=[p['row_id'] for p in parts], coordinate_statuses=checks,
        original_phases=[p['phase'] for p in parts], sequence_sha256=None, sequence_length=None,
        annotation_exceptions=[dict(row_id=p['row_id'], attributes={key:value for key,value in p['attrs'].items()
            if key in ['exception','transl_except','partial','start_range','end_range','pseudo','pseudogene','part','transl_table']}) for p in parts],
        order_source=None, ordered_feature_row_ids=[], status=None, flags=[])
    if set(checks) - {'within_original_sequence_bounds','documented_circular_virtual_coordinates'}:
        row['status'] = 'coordinate_review_required'
        return row, None
    if set(p['strand'] for p in parts) - {'+','-'}:
        row['status'] = 'unsupported_or_unknown_cds_strand'
        return row, None
    if any('part' in p['attrs'] and p['attrs']['part'] for p in parts):
        ranks = {}
        denominators = set()
        valid = True
        for p in parts:
            values = p['attrs'].get('part', [])
            if len(values) != 1 or not re.fullmatch('[1-9][0-9]*/[1-9][0-9]*', values[0]):
                valid = False
                break
            numerator, denominator = map(int, values[0].split('/'))
            ranks.setdefault(numerator, []).append(p)
            denominators.add(denominator)
        if not valid:
            row['status'] = 'incomplete_or_invalid_explicit_part_order'
            return row, None
        if denominators != {len(parts)} or set(ranks) != set(range(1,len(parts)+1)) or any(len(v)!=1 for v in ranks.values()):
            row['status'] = 'conflicting_explicit_part_order'
            return row, None
        ordered = [ranks[n][0] for n in range(1,len(parts)+1)]
        row['order_source'] = 'complete_explicit_part_X_over_Y'
    else:
        seqids = {p['seqid'] for p in parts}
        strands = {p['strand'] for p in parts}
        if len(seqids) != 1 or len(strands) != 1:
            row['status'] = 'multiple_sequences_or_strands_without_explicit_order'
            return row, None
        if len(parts)>1 and seqids & circular:
            row['status'] = 'circular_multipart_order_requires_review'
            return row, None
        ordered = sorted(parts,key=lambda p:(p['start'],p['end'],p['row_id']))
        if next(iter(strands)) == '-':
            ordered.reverse()
        row['order_source'] = 'single_sequence_strand_coordinate_order'
    intervals = sorted(parts,key=lambda p:(p['seqid'],p['start'],p['end']))
    if any(left['seqid']==right['seqid'] and right['start']<=left['end'] for left,right in zip(intervals,intervals[1:])):
        row['flags'].append('overlapping_parts_concatenated_without_repair')
    if any(p['phase']!='0' for p in ordered):
        row['flags'].append('original_nonzero_or_unknown_phase_not_trimmed')
    if any(p['attrs'].get('exception') or p['attrs'].get('transl_except') for p in ordered):
        row['flags'].append('annotation_exception_retained_not_corrected')
    chunks = []
    for p in ordered:
        data = genome[p['seqid']]
        start = (p['start']-1) % len(data)
        span = p['end']-p['start']+1
        # Build at most two slices, without duplicating an entire chromosome.
        first = data[start:min(start+span,len(data))]
        fragment = first + data[:span-len(first)]
        if p['strand']=='-':
            fragment = fragment.translate(COMPLEMENT)[::-1]
        chunks.append(fragment.upper())
    sequence = b''.join(chunks)
    row.update(status='unmodified_genomic_cds_candidate', ordered_feature_row_ids=[p['row_id'] for p in ordered],
               sequence_sha256=hashlib.sha256(sequence).hexdigest(), sequence_length=len(sequence))
    return row, sequence


def source_groups(connection, mode):
    grouped = defaultdict(dict)
    if mode=='verified_orf_coordinates':
        for ordinal,protein,seqid,start,end,strand,status,raw in connection.execute('SELECT * FROM provisional_orfs ORDER BY source_ordinal'):
            grouped[protein][('provisional_orf',ordinal)] = [dict(row_id=ordinal,seqid=seqid,start=start,end=end,
                strand=strand,phase='unknown',attrs={},original_orf_status=status,
                link_bases=['original_provisional_orf_coordinate'])]
        return grouped
    # Build features once and collect reference bases by physical feature row, not by reference occurrence.
    query='''SELECT product_id,feature_row_id,link_basis FROM cds_product_refs ORDER BY product_id,feature_row_id,link_basis'''
    references=defaultdict(lambda:defaultdict(list))
    for product,row_id,basis in connection.execute(query):
        references[product][row_id].append(basis)
    features={}
    for row_id,seqid,start,end,strand,phase,ident,raw in connection.execute("SELECT row_id,seqid,start,end,strand,phase,feature_id,attributes_json FROM features WHERE feature_type='CDS'"):
        features[row_id] = dict(row_id=row_id,seqid=seqid,start=start,end=end,strand=strand,phase=phase,
                               ident=ident,attrs=json.loads(raw))
    for product,refs in references.items():
        for row_id in sorted(refs):
            source=features[row_id]
            if mode=='ncbi_protein_gff':
                key=('ncbi_feature_id',source['ident']) if source['ident'] is not None else ('missing_unique_feature_id',row_id)
            elif mode=='transcript_gff':
                key=('published_transcript_parent_candidate',product)
            else:
                assert mode=='creolimax_gtf'
                transcripts=source['attrs'].get('transcript_id',[])
                key=('gtf_transcript_candidate',transcripts[0]) if len(transcripts)==1 else ('unresolved_gtf_transcript',row_id)
            grouped[product].setdefault(key,[]).append(dict(source,link_bases=refs[row_id]))
    return grouped


def record_match(target, values, count):
    positions=[n for n,value in enumerate(values) if value['sequence']==target['sequence']]
    status='one_exact_unmodified_genomic_cds_candidate'
    if target['product_id'] is None:status='missing_or_ambiguous_target_product_id'
    elif not values:status='no_genomic_cds_candidate'
    elif all(value['sequence'] is None for value in values):status='all_genomic_candidates_require_review'
    elif not positions:status='no_exact_unmodified_genomic_cds_match'
    elif len(positions)!=1:status='multiple_exact_genomic_candidates_require_review'
    elif count!=1:status='exact_genomic_candidate_multiple_target_records_require_review'
    elif 'gene_level_candidate_requires_review' in values[positions[0]]['link_bases']:status='exact_sequence_gene_level_candidate_requires_review'
    elif 'annotation_exception_retained_not_corrected' in values[positions[0]]['flags']:status='exact_genomic_candidate_with_annotation_exception_requires_review'
    return dict(**{k:v for k,v in target.items() if k!='sequence'},status=status,candidate_count=len(values),
        matching_candidate_indices=positions,exact_sequence_match=bool(positions),
        candidates=[{k:v for k,v in value.items() if k!='sequence'} for value in values])


def check_taxon(entry, report):
    assert report['taxon_id']==entry['taxon_id'] and report['mapping_mode']==entry['mapping_mode']
    verify(entry['pins']);verify(report['source_hashes'])
    root=Path(report['artifact_paths'][0]).parent
    assert report['artifact_paths']==[str(root/name) for name in ['feature_coordinates.tsv.gz','genomic_candidates.jsonl.gz','cds_record_comparisons.jsonl.gz','source_product_dispositions.jsonl.gz']]
    connection=sqlite3.connect(Path(entry['registry_report']['database']).resolve().as_uri()+'?mode=ro',uri=True)
    genome,sizes=genomic_sequences(entry)
    circular=set()
    for seqid,raw in connection.execute("SELECT seqid,attributes_json FROM features WHERE feature_type='region'"):
        if json.loads(raw).get('Is_circular')==['true']:circular.add(seqid)
    coordinate_counts,feature_counts=Counter(),Counter()
    with gzip.open(root/'feature_coordinates.tsv.gz','rt') as handle:
        rows=iter(csv.DictReader(handle,delimiter='\t'))
        for row_id,seqid,kind,start,end,strand in connection.execute('SELECT row_id,seqid,feature_type,start,end,strand FROM features ORDER BY row_id'):
            status=bounds(seqid,start,end,sizes,circular) if genome else 'original_genome_unavailable'
            expected=dict(source_kind='original_annotation_feature',row_id=str(row_id),seqid=seqid,feature_type=kind,
                start=str(start),end=str(end),strand=strand,original_sequence_length=str(sizes.get(seqid,'')),
                circular_annotation=str(seqid in circular),coordinate_status=status)
            assert next(rows)==expected
            feature_counts[kind]+=1;coordinate_counts[kind+':'+status]+=1
        for ordinal,protein,seqid,start,end,strand,status,raw in connection.execute('SELECT * FROM provisional_orfs ORDER BY source_ordinal'):
            observed=bounds(seqid,start,end,sizes,set()) if genome else 'original_genome_unavailable'
            assert next(rows)==dict(source_kind='original_provisional_orf',row_id=str(ordinal),seqid=seqid,feature_type='provisional_orf',
                start=str(start),end=str(end),strand=strand,original_sequence_length=str(sizes.get(seqid,'')),
                circular_annotation='False',coordinate_status=observed)
            feature_counts['provisional_orf_row']+=1;coordinate_counts['provisional_orf:'+observed]+=1
        assert next(rows,None) is None
    assert dict(feature_counts)==entry['registry_report']['feature_counts']==report['feature_counts']
    assert dict(coordinate_counts)==report['coordinate_counts']
    groups=source_groups(connection,entry['mapping_mode']);values={};candidate_counts=Counter()
    with gzip.open(root/'genomic_candidates.jsonl.gz','rt') as handle:
        rows=iter(map(json.loads,handle))
        for product in sorted(groups):
            candidates=[]
            for key,parts in sorted(groups[product].items(),key=lambda x:str(x[0])):
                if genome:expected,sequence=expected_candidate(parts,genome,sizes,circular)
                else:
                    expected=dict(status='original_genome_unavailable',flags=[],sequence_sha256=None,sequence_length=None,
                                  feature_row_ids=[p['row_id'] for p in parts]);sequence=None
                if key[0] in ('missing_unique_feature_id','unresolved_gtf_transcript'):
                    expected.update(status='unresolved_cds_feature_group_requires_review',sequence_sha256=None,sequence_length=None);sequence=None
                if any(p.get('original_orf_status','exact_translation')!='exact_translation' for p in parts):
                    expected['flags'].append('original_provisional_orf_review_status_retained')
                expected.update(group_key=list(key),link_bases=sorted({basis for p in parts for basis in p['link_bases']}))
                assert next(rows)==dict(product_id=product,**expected)
                candidate_counts[expected['status']]+=1;candidates.append(dict(expected,sequence=sequence))
            values[product]=candidates
        assert next(rows,None) is None
    assert dict(candidate_counts)==report['candidate_counts']
    targets=[]
    if entry['target_cds'] is not None:
        src=entry['target_cds']
        for ordinal,(title,sequence) in enumerate(fasta_bytes(src['path']),1):
            title=title.rstrip()  # Producer's Biopython description convention, not genomic header rewriting.
            if src['identifier_mode']=='ncbi_protein_id_header':
                identifiers=re.findall(r'\[protein_id=([^\]]+)\]',title)
                product=identifiers[0] if len(identifiers)==1 else None
            else:
                assert src['identifier_mode']=='exact_fasta_id';product=title.split()[0]
            upper=sequence.upper()
            targets.append(dict(ordinal=ordinal,cds_id=title.split()[0],description=title,product_id=product,
                sequence=upper,sequence_length=len(upper),sequence_sha256=hashlib.sha256(upper).hexdigest()))
        if 'expected_records' in src:assert len(targets)==src['expected_records']
    multiplicity=Counter(t['product_id'] for t in targets);target_counts=Counter();linked=defaultdict(list)
    with gzip.open(root/'cds_record_comparisons.jsonl.gz','rt') as handle:
        rows=iter(map(json.loads,handle))
        for target in targets:
            expected=record_match(target,values.get(target['product_id'],[]),multiplicity[target['product_id']])
            assert next(rows)==expected
            target_counts[expected['status']]+=1
            linked[target['product_id']].append(dict(cds_id=target['cds_id'],ordinal=target['ordinal'],status=expected['status']))
        assert next(rows,None) is None
    assert len(targets)==report['target_records'] and dict(target_counts)==report['target_status_counts']
    assert report['target_cds_source']==entry['target_cds'] and report['missing_target_reason']==entry['missing_target_reason']
    products=selected=0
    with gzip.open(root/'source_product_dispositions.jsonl.gz','rt') as handle:
        rows=iter(map(json.loads,handle))
        for ident,ordinal,description,length,digest,mapping,decision,included in connection.execute('SELECT * FROM products ORDER BY source_ordinal'):
            assert next(rows)==dict(protein_id=ident,source_ordinal=ordinal,sequence_sha256=digest,protein_length=length,
                original_mapping=json.loads(mapping),original_decision=json.loads(decision),selected_representative=bool(included),
                genomic_candidate_count=len(values.get(ident,[])),original_cds_targets=linked.get(ident,[]),
                missing_target_reason=entry['missing_target_reason'] if entry['target_cds'] is None else None)
            products+=1;selected+=included
        assert next(rows,None) is None
    connection.close()
    assert products==report['source_products']==entry['registry_report']['source_products']
    assert selected==report['selected_representatives']==entry['registry_report']['selected_representatives']
    verify(entry['pins'])
    return dict(taxon_id=entry['taxon_id'],source_products=products,selected_representatives=selected,
                target_records=len(targets),feature_counts=dict(feature_counts),coordinate_counts=dict(coordinate_counts),
                candidate_counts=dict(candidate_counts),target_status_counts=dict(target_counts))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True);parser.add_argument('--receipt',type=Path,required=True)
    args=parser.parse_args();plan=json.loads(args.plan.read_text());verify(plan['pins'])
    original=json.loads(Path(plan['producer_plan']).read_text())
    producer=json.loads(Path(plan['producer_receipt']).read_text())
    transport=json.loads(Path(plan['producer_transport']).read_text())
    assert producer['status']=='complete_full_genome_annotation_cds_comparison_pending_independent_readback'
    assert transport['status']=='verified_original_software_wait_and_whole_wrapper_payloads' and transport['original_tool_terminal_exit_code']==0
    verify(transport['source_hashes'])
    reports={r['taxon_id']:r for r in producer['taxa_reports']}
    assert len(reports)==len(original['entries'])==526
    root=Path(plan['output']);root.mkdir(parents=True,exist_ok=False);assert not args.receipt.exists()
    start=time.monotonic();proofs=[]
    with ProcessPoolExecutor(max_workers=plan['cpu']) as pool:
        futures=[pool.submit(check_taxon,e,reports[e['taxon_id']]) for e in original['entries']]
        for future in as_completed(futures):
            proofs.append(future.result())
            state=dict(stage='independent_full_genomic_cds_and_source_replay',completed_taxa=len(proofs),expected_taxa=526,
                target_records=sum(r['target_records'] for r in proofs),elapsed_seconds=time.monotonic()-start)
            temporary=root/'state.partial';temporary.write_text(json.dumps(state,indent=2)+'\n');temporary.replace(root/'state.json')
            print(json.dumps(state),flush=True)
    assert sum(r['source_products'] for r in proofs)==5927745
    assert sum(r['selected_representatives'] for r in proofs)==5815847
    pins=dict(plan['pins'])
    for p,d in producer['source_hashes'].items():bind(pins,p,d)
    for p in (args.plan,plan['producer_receipt'],plan['producer_transport']):bind(pins,p)
    verify(pins)
    result=dict(status='passed_full_independent_genome_annotation_cds_and_source_product_replay',
        checked_utc=datetime.now(timezone.utc).isoformat(),taxa=526,source_products=5927745,selected_representatives=5815847,
        target_records=sum(r['target_records'] for r in proofs),taxa_checks=sorted(proofs,key=lambda r:r['taxon_id']),
        elapsed_seconds=time.monotonic()-start,source_hashes=pins,scientific_eligibility=False,gpu=False,new_predictions=0,
        scope='All original coordinates/candidates/CDS targets/source products and exclusions replayed with own FASTA decoder, '
              'reference grouping, explicit ordering and byte-table IUPAC complement, without producer imports. '
              'Python/gzip/SQLite/hash libraries and qualified original sources are shared. Exact DNA agreement is not gene-copy, '
              'contamination, translation, selection, predictor accuracy or evolutionary acceptance.')
    with args.receipt.open('x') as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_hashes','taxa_checks')},indent=2))


if __name__=='__main__':main()
