#!/usr/bin/env python3
"""Offline literal tests of assembly DNA checksum, FASTA and source-path contracts."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
from unittest.mock import patch

from retrieve_full_assembly_dna_v1 import Budget, acquire, digest, publisher_md5, read_fasta


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--receipt',type=Path,required=True)
    args=parser.parse_args();assert not args.receipt.exists();checks=[]
    with tempfile.TemporaryDirectory(prefix='genome-literal-',dir='results') as directory:
        root=Path(directory).resolve();dna=b'>contig_one original description\nacGTNNRY\n>contig_two\nGGCC\n'
        fasta=root/'literal.fna.gz';fasta.write_bytes(gzip.compress(dna,mtime=0))
        stats=read_fasta(fasta,root/'contigs.jsonl')
        assert stats['total_length']==12 and stats['fasta_records']==2
        assert (stats['record_N50'],stats['record_L50'])==(8,1)
        assert stats['lowercase_bases']==2 and stats['N_bases']==2 and stats['other_ambiguous_bases']==2
        assert stats['gc_percent_canonical']==75 and stats['removed_sequence_whitespace']==0
        rows=[json.loads(l) for l in (root/'contigs.jsonl').read_text().splitlines()]
        assert rows[0]['sequence_sha256']==hashlib.sha256(b'acGTNNRY').hexdigest()
        assert rows[0]['uppercase_sequence_sha256']==hashlib.sha256(b'ACGTNNRY').hexdigest()
        assert rows[0]['description']=='contig_one original description'
        checks.append('exact_literal_dna_lengths_n50_case_ambiguity_and_hashes')
        for label,content in [('duplicate_id',b'>x\nAC\n>x\nGT\n'),('empty_record',b'>x\n>y\nAC\n'),
                              ('non_iupac',b'>x\nACU\n'),('before_header',b'AC\n>x\nGT\n'),
                              ('all_unknown',b'>x\nNN\n')]:
            p=root/(label+'.fna');p.write_bytes(content)
            try:read_fasta(p,root/(label+'.jsonl'))
            except ValueError:checks.append('reject_'+label)
            else:raise AssertionError('Accepted '+label)
        broken=root/'broken.fna.gz';broken.write_bytes(fasta.read_bytes()[:-4])
        try:read_fasta(broken,root/'broken.jsonl')
        except (EOFError,OSError):checks.append('reject_missing_gzip_trailer')
        else:raise AssertionError('Accepted truncated gzip')
        md5=digest(fasta,'md5');listing=(md5+'  ./literal.fna.gz\n').encode()
        assert publisher_md5(listing,'literal.fna.gz')==md5;checks.append('exact_publisher_checksum_match')
        for value in [b'',listing+listing]:
            try:publisher_md5(value,'literal.fna.gz')
            except ValueError:checks.append('reject_absent_or_duplicate_publisher_checksum')
            else:raise AssertionError('Accepted ambiguous checksum')
        url='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/000/001/GCA_000000001.1_OFFLINE/GCA_000000001.1_OFFLINE_genomic.fna.gz'
        filename=url.rsplit('/',1)[1];payload=fasta.read_bytes()
        checksum=(hashlib.md5(payload).hexdigest()+'  ./'+filename+'\n').encode()
        entry=dict(source='ncbi_exact_assembly_version',taxon=dict(taxon_id='OFFLINE',assembly_accession='GCA_000000001.1',genome_url=url,proteome_url=url.removesuffix('_genomic.fna.gz')+'_protein.faa.gz'),statistics_receipt=dict(assembly_accession='GCA_000000001.1',all_assembly_metrics={'total-length':12}))
        class Response:
            def __init__(self,data):self.data=data;self.position=0;self.headers={'Content-Length':str(len(data))}
            def __enter__(self):return self
            def __exit__(self,*args):return False
            def read(self,size):
                result=self.data[self.position:self.position+size];self.position+=len(result);return result
        def serve(request,timeout):
            return Response(checksum if request.full_url.endswith('/md5checksums.txt') else payload)
        with patch('retrieve_full_assembly_dna_v1.urlopen',side_effect=serve) as requests:
            report=acquire(entry,root,Budget(root,1024*1024,0),1024*1024)
            assert requests.call_count==2 and report['attempted_http_requests']==2
            assert report['status']=='verified_publisher_bound_genomic_dna'
            assert report['sha256']==hashlib.sha256(payload).hexdigest() and report['publisher_md5']==hashlib.md5(payload).hexdigest()
            assert report['assembly_scope_length_status']=='matches_deposited_all_scope'
            assert report['fasta']==stats
            assert json.loads(Path(report['receipt_path']).read_text())=={k:v for k,v in report.items() if k not in ('receipt_path','receipt_sha256')}
            checks.append('literal_two_request_download_checksum_body_receipt_and_absolute_paths')
        budget=Budget(root,4,0);budget.account(4)
        try:budget.account(1)
        except RuntimeError:checks.append('reject_global_download_budget_overrun')
        else:raise AssertionError('Accepted budget overrun')
    pins={str(p):digest(p) for p in (Path(__file__),Path('scripts/retrieve_full_assembly_dna_v1.py'))}
    result=dict(status='passed_offline_full_assembly_dna_contract_controls',checks=checks,live_http_requests=0,
                source_hashes=pins,scientific_eligibility=False,gpu=False,new_predictions=0,
                scope='Synthetic DNA and HTTP responses only; exact streaming record/hash/statistic controls, malformed/truncated data and checksum ambiguity rejection, full receipt/path contract and resource budget. No corpus pilot or live downloads.')
    with args.receipt.open('x') as handle:json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
