#!/usr/bin/env python3
"""Literal genomic joins and complete small-source product/CDS accounting; no corpus pilot."""
import argparse
from collections import Counter
import copy
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import tempfile

from audit_full_genome_annotation_cds_v1 import audit_taxon, compare_target
from build_full_annotation_coordinate_registry_v2 import build_taxon
from check_full_annotation_coordinate_registry_cases_v2 import entry
from genomic_cds_join_v1 import coordinate_status, reconstruct
from reference_measurement_union_sources import bind


def part(number, start, end, strand='+', attrs=None, seqid='chr', phase='0'):
    return dict(row_id=number, seqid=seqid, start=start, end=end, strand=strand,
                phase=phase, attrs={} if attrs is None else attrs)


def literal_joins():
    genomes = {'chr': 'AAAACCCCGGGGTTTT', 'circ': 'AAGC'}
    rows = [part(1,1,4), part(2,9,12,phase='2')]
    report, sequence = reconstruct(rows, genomes, set())
    assert sequence == 'AAAAGGGG' and report['original_phases'] == ['0','2']
    assert 'original_nonzero_or_unknown_phase_not_trimmed' in report['flags']
    checks = ['plus_join_and_nonzero_phase_preserved']
    rows = [part(1,1,4,'-'),part(2,9,12,'-')]
    report, sequence = reconstruct(rows, genomes, set())
    assert sequence == 'CCCCTTTT' and report['ordered_feature_row_ids'] == [2,1]
    checks.append('minus_join_has_biological_coordinate_order')
    report, sequence = reconstruct([part(1,1,4),part(2,3,6)],genomes,set())
    assert sequence == 'AAAAAACC' and 'overlapping_parts_concatenated_without_repair' in report['flags']
    checks.append('overlapping_bases_retained_without_merging')
    report, sequence = reconstruct([part(1,3,6,seqid='circ')],genomes,{'circ'})
    assert sequence == 'GCAA' and report['coordinate_statuses'] == ['documented_circular_virtual_coordinates']
    report, sequence = reconstruct([part(1,3,6,'-',seqid='circ')],genomes,{'circ'})
    assert sequence == 'TTGC'
    checks.append('documented_single_circular_virtual_interval_both_strands')
    report, sequence = reconstruct([part(1,3,4,seqid='circ'),part(2,1,2,seqid='circ')],genomes,{'circ'})
    assert sequence is None and report['status'] == 'circular_multipart_order_requires_review'
    checks.append('circular_multipart_without_explicit_order_not_guessed')
    rows = [part(1,3,4,seqid='circ',attrs={'part':['1/2']}),part(2,1,2,seqid='circ',attrs={'part':['2/2']})]
    report, sequence = reconstruct(rows,genomes,{'circ'})
    assert sequence == 'GCAA' and report['order_source'] == 'complete_explicit_part_X_over_Y'
    checks.append('complete_explicit_part_order_preserved')
    rows[1]['attrs']['part'] = ['1/2']
    report, sequence = reconstruct(rows,genomes,{'circ'})
    assert sequence is None and report['status'] == 'conflicting_explicit_part_order'
    checks.append('duplicate_explicit_part_rank_rejected')
    report, sequence = reconstruct([part(1,1,4,attrs={'part':['1/2']}),part(2,9,12)],genomes,set())
    assert sequence is None and report['status'] == 'incomplete_or_invalid_explicit_part_order'
    checks.append('incomplete_explicit_part_order_retained_for_review')
    report, sequence = reconstruct([part(1,1,4),part(2,9,12,'-')],genomes,set())
    assert sequence is None and report['status'] == 'multiple_sequences_or_strands_without_explicit_order'
    checks.append('mixed_strand_unordered_join_retained_for_review')
    assert coordinate_status('missing',1,2,{k:len(v) for k,v in genomes.items()},set()) == 'missing_genomic_sequence_id'
    assert coordinate_status('circ',3,6,{'circ':4},set()) == 'out_of_bounds_without_circular_annotation'
    assert coordinate_status('circ',1,7,{'circ':4},{'circ'}) == 'circular_span_requires_review'
    checks.append('missing_sequence_non_circular_bounds_and_excess_circular_span')
    report, sequence = reconstruct([part(1,1,4,attrs={'exception':['RNA editing']})],genomes,set())
    assert sequence == 'AAAA' and report['annotation_exceptions'][0]['attributes']['exception'] == ['RNA editing']
    candidate = dict(report, sequence=sequence, link_bases=['explicit_protein_id'])
    target = dict(ordinal=1,cds_id='cds1',description='cds1',product_id='p',sequence='AAAA',sequence_length=4,sequence_sha256=hashlib.sha256(b'AAAA').hexdigest())
    result = compare_target(target,[candidate],1)
    assert result['exact_sequence_match'] and result['status'] == 'exact_genomic_candidate_with_annotation_exception_requires_review'
    checks.append('exact_sequence_with_annotation_exception_keeps_review_status')
    result = compare_target(target,[candidate,candidate],1)
    assert result['matching_candidate_indices'] == [0,1] and result['status'] == 'multiple_exact_genomic_candidates_require_review'
    checks.append('multiple_matching_loci_not_silently_selected')
    result = compare_target(target,[candidate],2)
    assert result['status'] == 'exact_genomic_candidate_multiple_target_records_require_review'
    checks.append('multiple_original_cds_records_not_silently_selected')
    target['sequence'] = 'CCCC'
    assert compare_target(target,[candidate],1)['status'] == 'no_exact_unmodified_genomic_cds_match'
    checks.append('target_mismatch_not_repaired')
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args(); assert not args.receipt.exists()
    checks = literal_joins()
    with tempfile.TemporaryDirectory(prefix='literal-genomic-cds-',dir='results') as folder:
        root = Path(folder)
        for mode in ('ncbi_protein_gff','transcript_gff','creolimax_gtf','verified_orf_coordinates'):
            case = root / mode;case.mkdir()
            source = entry(case,mode)
            index = case/'index';index.mkdir()
            registry = build_taxon(source,index,dict(emergency_free_disk_gib=1,output_allowance_gib=1))
            genome = case/'genome.fna';genome.write_text('>chr original  \n'+'ACGT'*25+'\n')
            metadata = case/'contigs.jsonl';metadata.write_text(json.dumps(dict(sequence_id='chr',description='chr original  ',length=100,
                sequence_sha256=hashlib.sha256(('ACGT'*25).encode()).hexdigest(),uppercase_sequence_sha256=hashlib.sha256(('ACGT'*25).encode()).hexdigest()))+'\n')
            cds = case/'original_cds.fna'
            expected = {'ncbi_protein_gff':'GTACGTACGTACGTACCGTACGT','transcript_gff':'CGTACGTACGTA',
                        'creolimax_gtf':'TACGTACGTACG','verified_orf_coordinates':'ACGTACGTACGT'}[mode]
            cds.write_text('>cds1 [protein_id=p,1]\n'+expected+'\n' if mode == 'ncbi_protein_gff' else '>p,1\n'+expected+'\n')
            inputs = copy.deepcopy(source)
            inputs['registry_report'] = registry
            inputs['genome_report'] = dict(status='verified_publisher_bound_genomic_dna',path=str(genome),contig_metadata=str(metadata))
            inputs['target_cds'] = dict(path=str(cds),identifier_mode='ncbi_protein_id_header' if mode == 'ncbi_protein_gff' else 'exact_fasta_id',expected_records=1)
            inputs['missing_target_reason'] = None
            for path in (genome,metadata,cds,Path(registry['database'])):bind(inputs['pins'],path)
            audit = case/'audit';audit.mkdir()
            report = audit_taxon(inputs,audit,dict(emergency_free_disk_gib=1))
            assert report['source_products'] == 3 and report['selected_representatives'] == 2 and report['target_records'] == 1
            with gzip.open(audit/'literal/cds_record_comparisons.jsonl.gz','rt') as handle:row=json.loads(handle.readline())
            assert row['exact_sequence_match'] and len(row['matching_candidate_indices']) == 1
            with gzip.open(audit/'literal/source_product_dispositions.jsonl.gz','rt') as handle:products=list(map(json.loads,handle))
            assert len(products) == 3 and products[1]['protein_id'] == 'alt' and not products[1]['selected_representative']
            if mode == 'verified_orf_coordinates':
                assert report['coordinate_counts']['provisional_orf:out_of_bounds_without_circular_annotation'] == 1
            checks.append(mode+':complete_literal_genome_source_cds_and_product_accounting')
            # Transcript data are explicitly not used as a coding-DNA target.
            if mode == 'creolimax_gtf':
                inputs['target_cds']=None;inputs['missing_target_reason']='Published transcript contains unqualified coding boundaries; no original CDS target'
                other=case/'no_cds_audit';other.mkdir();r=audit_taxon(inputs,other,dict(emergency_free_disk_gib=1))
                assert r['target_records']==0 and r['source_products']==3 and r['missing_target_reason']
                checks.append('creolimax_transcript_not_substituted_for_cds_target')
    pins={}
    for name in ('genomic_cds_join_v1.py','audit_full_genome_annotation_cds_v1.py',Path(__file__).name):
        bind(pins,Path('scripts')/name)
    result=dict(status='passed_offline_literal_genome_annotation_cds_controls',checked_utc=datetime.now(timezone.utc).isoformat(),
        checks=checks,count=len(checks),source_hashes=pins,scientific_eligibility=False,corpus_pilot=False,gpu=False,new_predictions=0,
        scope='Hand-calculated plus/minus/overlap/circular/order/exception/multiplicity joins and every literal source product '
              'in four source modes. No source-corpus comparison, gene-copy, CDS selection or biological acceptance.')
    with args.receipt.open('x') as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
