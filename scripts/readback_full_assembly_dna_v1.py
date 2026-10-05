#!/usr/bin/env python3
"""Independently reconstruct every acquired genomic FASTA and all taxon dispositions."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from fractions import Fraction
import gzip
import hashlib
import json
import math
import multiprocessing
from pathlib import Path

from Bio import SeqIO

from reference_measurement_union_sources import bind, verify


def digest(path, algorithm='sha256'):
    h=hashlib.new(algorithm)
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def measure(path, contigs):
    """Biopython sequence parsing and independent integer/rational summaries."""
    sizes=[];counts=Counter();seen=set();lower=0;whitespace=0
    opener=gzip.open if Path(path).suffix=='.gz' else open
    with opener(path,'rt') as sequences,opener(path,'rb') as header_file,Path(contigs).open() as saved:
        headers=(line[1:].rstrip(b'\r\n').decode('utf-8') for line in header_file if line.startswith(b'>'))
        rows=(json.loads(line) for line in saved)
        for record in SeqIO.parse(sequences,'fasta'):
            original_header=next(headers);row=next(rows);original=str(record.seq);upper=original.upper()
            assert record.id not in seen and original and not (set(upper)-set('ACGTRYSWKMBDHVN'))
            assert record.id==original_header.split()[0] and row==dict(sequence_id=record.id,description=original_header,
                length=len(original),sequence_sha256=hashlib.sha256(original.encode('ascii')).hexdigest(),
                uppercase_sequence_sha256=hashlib.sha256(upper.encode('ascii')).hexdigest())
            seen.add(record.id);sizes.append(len(original));counts.update(upper)
            lower+=sum(n for c,n in Counter(original).items() if c.islower())
        assert next(rows,None) is None and next(headers,None) is None
    # Original whitespace accounting independently reads all raw sequence lines.
    with opener(path,'rb') as handle:
        for raw in handle:
            line=raw.rstrip(b'\r\n')
            if not line.strip() or line.startswith(b'>'):continue
            whitespace+=sum(1 for byte in line if chr(byte).isspace())
    assert sizes
    canonical=sum(counts[c] for c in 'ACGT');assert canonical>0
    total=sum(sizes);ordered=sorted(sizes,reverse=True);threshold=(total+1)//2
    accumulated=0
    for rank,size in enumerate(ordered,1):
        accumulated+=size
        if accumulated>=threshold:n50,l50=size,rank;break
    return dict(fasta_records=len(sizes),total_length=total,longest_record=ordered[0],record_N50=n50,record_L50=l50,
                N_bases=counts['N'],N_fraction=float(Fraction(counts['N'],total)),
                other_ambiguous_bases=total-canonical-counts['N'],
                gc_percent_canonical=float(Fraction(100*(counts['G']+counts['C']),canonical)),
                lowercase_bases=lower,removed_sequence_whitespace=whitespace,
                base_counts=dict(sorted(counts.items())),
                scope='All deposited FASTA records; no gap splitting, organelle exclusion, taxonomic admission or contamination inference')


def compare_statistics(original,reconstructed):
    assert original.keys()==reconstructed.keys()
    for key,value in reconstructed.items():
        if key in ('N_fraction','gc_percent_canonical'):
            assert math.isfinite(original[key]) and abs(original[key]-value)<=2*math.ulp(value)
        else:assert original[key]==value,key


def check_taxon(entry,report):
    taxon=entry['taxon'];assert report['taxon']==taxon and report['source']==entry['source']
    if report['status']!='verified_publisher_bound_genomic_dna':
        assert report['status']=='genome_acquisition_or_format_error' and report['error_type'] and report['error']
        return dict(taxon_id=taxon['taxon_id'],status='original_error_retained_not_admitted',bases=0)
    path=Path(report['path']);assert digest(path)==report['sha256'] and digest(path,'md5')==report['publisher_md5']
    assert path.stat().st_size==report['bytes'] and report['url']==taxon['genome_url']
    assert digest(report['receipt_path'])==report['receipt_sha256']
    assert json.loads(Path(report['receipt_path']).read_text())=={k:v for k,v in report.items() if k not in ('receipt_path','receipt_sha256')}
    assert digest(report['contig_metadata'])==report['contig_metadata_sha256']
    if entry['source']=='ncbi_exact_assembly_version':
        assert report['checksum_url']==taxon['genome_url'].rsplit('/',1)[0]+'/md5checksums.txt'
        assert digest(report['checksum_listing'])==report['checksum_listing_sha256']
        matches=[];name=taxon['genome_url'].rsplit('/',1)[1]
        for line in Path(report['checksum_listing']).read_text().splitlines():
            fields=line.split()
            if len(fields)==2 and fields[1].removeprefix('./')==name:matches.append(fields[0].lower())
        assert matches==[report['publisher_md5']]
        assert entry['statistics_receipt']['assembly_accession']==taxon['assembly_accession']
        assert name.startswith(taxon['assembly_accession']+'_')
    else:
        original=entry['external_receipt'];assert report['external_receipt']==original
        assert path.resolve()==Path(original['path']).resolve() and report['sha256']==original['sha256']
        assert report['publisher_md5']==original['publisher_md5'] and report['new_download_bytes']==0
    reconstructed=measure(path,report['contig_metadata']);compare_statistics(report['fasta'],reconstructed)
    if entry['source']=='ncbi_exact_assembly_version':
        deposited=int(entry['statistics_receipt']['all_assembly_metrics']['total-length'])
        difference=reconstructed['total_length']-deposited
        assert report['deposited_whole_assembly_length']==deposited
        assert report['fasta_minus_deposited_whole_assembly_length']==difference
        assert report['assembly_scope_length_status']==('matches_deposited_all_scope' if difference==0 else 'requires_assembly_scope_review')
    return dict(taxon_id=taxon['taxon_id'],status='independently_reconstructed_original_genome',
                bases=reconstructed['total_length'],fasta_records=reconstructed['fasta_records'])


def self_test():
    import tempfile
    checks=[]
    with tempfile.TemporaryDirectory(prefix='genome-reader-literal-',dir='results') as folder:
        root=Path(folder);dna=root/'literal.fasta';dna.write_bytes(b'>x original\nacGTNNRY\n>y\nGGCC\n')
        contigs=root/'contigs.jsonl'
        data=[]
        for identifier,description,seq in [('x','x original','acGTNNRY'),('y','y','GGCC')]:
            data.append(dict(sequence_id=identifier,description=description,length=len(seq),
                             sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),
                             uppercase_sequence_sha256=hashlib.sha256(seq.upper().encode()).hexdigest()))
        contigs.write_text(''.join(json.dumps(row)+'\n' for row in data))
        value=measure(dna,contigs)
        assert (value['total_length'],value['record_N50'],value['record_L50'],value['lowercase_bases'],value['N_bases'],value['other_ambiguous_bases'],value['gc_percent_canonical'])==(12,8,1,2,2,2,75)
        compare_statistics(value,value);checks.append('independent_literal_records_hashes_n50_case_ambiguity_and_rational_metrics')
        for key in ['total_length','record_N50','lowercase_bases','N_bases','gc_percent_canonical']:
            bad=dict(value);bad[key]+=1
            try:compare_statistics(bad,value)
            except AssertionError:checks.append('reject_changed_'+key)
            else:raise AssertionError('Accepted changed '+key)
        data[0]['sequence_sha256']='changed';contigs.write_text(''.join(json.dumps(row)+'\n' for row in data))
        try:measure(dna,contigs)
        except AssertionError:checks.append('reject_changed_per_contig_sequence_digest')
        else:raise AssertionError('Accepted changed contig digest')
    return checks


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path);parser.add_argument('--receipt',type=Path,required=True)
    parser.add_argument('--self-test',action='store_true');args=parser.parse_args();assert not args.receipt.exists()
    checks=self_test()
    if args.self_test:
        result=dict(status='passed_literal_full_genome_independent_reader_controls',checks=checks,
                    source_hashes={str(Path(__file__)):digest(Path(__file__))},scientific_eligibility=False)
    else:
        assert args.plan is not None;plan=json.loads(args.plan.read_text());pins=dict(plan['pins']);verify(pins)
        source_plan=json.loads(Path(plan['producer_plan']).read_text())
        producer=json.loads(Path(plan['producer_receipt']).read_text());transport=json.loads(Path(plan['producer_transport']).read_text())
        assert producer['status']=='completed_all_selected_assembly_dna_attempts_pending_independent_readback'
        assert transport['status']=='verified_original_software_wait_and_whole_wrapper_payloads'
        assert transport['validation_sha256']==digest(plan['producer_receipt']) and transport['original_tool_terminal_exit_code']==0
        assert transport['whole_wrapper_initial_and_terminal_payloads_matched']
        assert transport['manager_start_records']==transport['manager_completion_records']==1
        root=Path(plan['output']);assert not root.exists();root.mkdir(parents=True)
        records=[json.loads(line) for line in (Path(source_plan['output'])/'taxon_dispositions.jsonl').read_text().splitlines()]
        by_taxon={r['taxon']['taxon_id']:r for r in records};entries=source_plan['entries']
        assert len(by_taxon)==len(records)==len(entries)==526
        assert by_taxon.keys()=={e['taxon']['taxon_id'] for e in entries}
        proofs=[];counts=Counter();bases=0
        with ProcessPoolExecutor(max_workers=plan['cpu'],mp_context=multiprocessing.get_context('spawn')) as pool:
            futures=[pool.submit(check_taxon,e,by_taxon[e['taxon']['taxon_id']]) for e in entries]
            for future in as_completed(futures):
                proof=future.result();proofs.append(proof);counts[proof['status']]+=1;bases+=proof['bases']
                tmp=root/'state.tmp';tmp.write_text(json.dumps(dict(completed_taxa=len(proofs),expected_taxa=526,counts=dict(counts)),indent=2)+'\n');tmp.replace(root/'state.json')
        assert counts['independently_reconstructed_original_genome']==producer['counts'].get('verified_publisher_bound_genomic_dna',0)
        assert counts['original_error_retained_not_admitted']==producer['counts'].get('genome_acquisition_or_format_error',0)
        assert bases==producer['verified_genome_bases']
        for mapping in (source_plan['pins'],producer['source_hashes'],transport['source_hashes']):
            for path,h in mapping.items():bind(pins,path,h)
        for path in (args.plan,Path(plan['producer_plan']),Path(plan['producer_receipt']),Path(plan['producer_transport']),Path(__file__)):bind(pins,path,digest(path))
        verify(pins)
        result=dict(status='passed_full_selected_genome_disposition_and_original_fasta_reconstruction',
                    checked_utc=datetime.now(timezone.utc).isoformat(),taxa=526,counts=dict(counts),verified_genome_bases=bases,
                    checks=checks,proofs=sorted(proofs,key=lambda p:p['taxon_id']),source_hashes=pins,
                    scientific_eligibility=False,gpu=False,new_predictions=0,
                    scope='Every original taxon/error disposition and available publisher-bound genome, every original '
                          'FASTA record, sequence/case digest, base/length/N50/ambiguity and rational metrics independently '
                          'reconstructed with separate Biopython parser; publisher listings/source version rechecked. '
                          'GFF/CDS/annotation, flanks, contamination/haplotig/taxon identity and biology remain separate. '
                          'No missing genome is excluded or considered available; no source files are rewritten.')
    with args.receipt.open('x') as handle:json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_hashes','proofs')},indent=2))


if __name__=='__main__':main()
