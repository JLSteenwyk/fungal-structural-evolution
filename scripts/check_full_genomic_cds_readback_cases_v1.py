#!/usr/bin/env python3
"""Literal full-source replay and deliberately rebound genomic/CDS corruption checks."""
import argparse
import copy
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import tempfile

from audit_full_genome_annotation_cds_v1 import audit_taxon
from build_full_annotation_coordinate_registry_v2 import build_taxon
from check_full_annotation_coordinate_registry_cases_v2 import entry
from readback_full_genome_annotation_cds_v1 import COMPLEMENT, check_taxon
from reference_measurement_union_sources import bind


def literal_case(root, mode, duplicate_gtf=False):
    source=entry(root,mode)
    if duplicate_gtf:
        annotation=Path(source['annotation_path'])
        annotation.write_text(annotation.read_text().replace('gene_id "g"','gene_id "p,1"'))
        source['pins'][str(annotation)]=hashlib.sha256(annotation.read_bytes()).hexdigest()
    index=root/'index';index.mkdir()
    registry=build_taxon(source,index,dict(emergency_free_disk_gib=1,output_allowance_gib=1))
    genome=root/'genome.fna';genome.write_text('>chr original  \n'+'ACGT'*25+'\n')
    metadata=root/'contigs.jsonl'
    digest=hashlib.sha256(('ACGT'*25).encode()).hexdigest()
    metadata.write_text(json.dumps(dict(sequence_id='chr',description='chr original  ',length=100,
        sequence_sha256=digest,uppercase_sequence_sha256=digest))+'\n')
    cds=root/'original_cds.fna'
    expected={'ncbi_protein_gff':'GTACGTACGTACGTACCGTACGT','transcript_gff':'CGTACGTACGTA',
              'creolimax_gtf':'TACGTACGTACG','verified_orf_coordinates':'ACGTACGTACGT'}[mode]
    cds.write_text('>cds1 [protein_id=p,1]\n'+expected+'\n' if mode=='ncbi_protein_gff' else '>p,1\n'+expected+'\n')
    result=dict(source,registry_report=registry,genome_report=dict(status='verified_publisher_bound_genomic_dna',
        path=str(genome),contig_metadata=str(metadata)),target_cds=dict(path=str(cds),
        identifier_mode='ncbi_protein_id_header' if mode=='ncbi_protein_gff' else 'exact_fasta_id',expected_records=1),
        missing_target_reason=None)
    for path in (genome,metadata,cds,Path(registry['database'])):bind(result['pins'],path)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt',type=Path,required=True);args=parser.parse_args();assert not args.receipt.exists()
    checks=[]
    assert b'ACGTRYSWKMBDHVNacgtryswkmbdhvn'.translate(COMPLEMENT)[::-1]==b'nbdhvkmwsryacgtNBDHVKMWSRYACGT'
    checks.append('hand_calculated_all_IUPAC_complements_and_case')
    with tempfile.TemporaryDirectory(prefix='literal-genomic-replay-',dir='results') as folder:
        root=Path(folder);originals={}
        for mode in ('ncbi_protein_gff','transcript_gff','creolimax_gtf','verified_orf_coordinates'):
            case=root/mode;case.mkdir();source=literal_case(case,mode)
            output=case/'audit';output.mkdir()
            report=audit_taxon(source,output,dict(emergency_free_disk_gib=1))
            proof=check_taxon(source,report)
            assert proof['source_products']==3 and proof['selected_representatives']==2 and proof['target_records']==1
            originals[mode]=(source,report)
            checks.append(mode+':every_literal_coordinate_candidate_CDS_product_replayed')
        case=root/'same_gtf_gene_transcript_identifier';case.mkdir();source=literal_case(case,'creolimax_gtf',True)
        output=case/'audit';output.mkdir();report=audit_taxon(source,output,dict(emergency_free_disk_gib=1))
        check_taxon(source,report)
        with gzip.open(output/'literal/genomic_candidates.jsonl.gz','rt') as handle:
            rows=list(map(json.loads,handle))
        candidate=next(r for r in rows if r['product_id']=='p,1')
        assert candidate['feature_row_ids']==[2] and len(candidate['link_bases'])==2
        checks.append('duplicate_gene_and_transcript_reference_does_not_duplicate_exon')
        source,report=originals['ncbi_protein_gff']
        unavailable=copy.deepcopy(source)
        unavailable['genome_report']=dict(status='genome_acquisition_or_format_error',error_type='LiteralUnavailable',error='Literal source exclusion')
        output=root/'missing_genome';output.mkdir();r=audit_taxon(unavailable,output,dict(emergency_free_disk_gib=1))
        p=check_taxon(unavailable,r);assert p['source_products']==3 and p['target_records']==1
        assert r['target_status_counts']=={'all_genomic_candidates_require_review':1}
        checks.append('unavailable_genome_keeps_all_source_and_target_records')
        source_gtf,_=originals['creolimax_gtf'];no_cds=copy.deepcopy(source_gtf)
        no_cds['target_cds']=None;no_cds['missing_target_reason']='No qualified original coding-DNA target; mRNA is not substituted'
        output=root/'unqualified_cds';output.mkdir();r=audit_taxon(no_cds,output,dict(emergency_free_disk_gib=1))
        p=check_taxon(no_cds,r);assert p['target_records']==0 and p['source_products']==3
        checks.append('unqualified_CDS_target_keeps_products_without_transcript_substitution')
        def relevant(rows):
            return next(r for r in rows if r['product_id']=='p,1')
        mutations={
            'coordinate_end':('feature_coordinates.tsv.gz',None),
            'candidate_sequence_hash':('genomic_candidates.jsonl.gz',lambda rows:relevant(rows).update(sequence_sha256='0'*64)),
            'original_phase':('genomic_candidates.jsonl.gz',lambda rows:relevant(rows)['original_phases'].__setitem__(0,'0')),
            'explicit_order':('genomic_candidates.jsonl.gz',lambda rows:relevant(rows).update(ordered_feature_row_ids=list(reversed(relevant(rows)['ordered_feature_row_ids'])))),
            'lost_exception':('genomic_candidates.jsonl.gz',lambda rows:relevant(rows).update(annotation_exceptions=[])),
            'target_sequence_hash':('cds_record_comparisons.jsonl.gz',lambda rows:rows[0].update(sequence_sha256='0'*64)),
            'matching_locus':('cds_record_comparisons.jsonl.gz',lambda rows:rows[0].update(matching_candidate_indices=[])),
            'target_header':('cds_record_comparisons.jsonl.gz',lambda rows:rows[0].update(description='changed source header')),
            'alternative_selection':('source_product_dispositions.jsonl.gz',lambda rows:rows[1].update(selected_representative=True)),
            'missing_unresolved_product':('source_product_dispositions.jsonl.gz',lambda rows:rows.pop())}
        source,report=originals['ncbi_protein_gff'];original_root=Path(report['artifact_paths'][0]).parent
        for name,(filename,mutate) in mutations.items():
            copied=root/('corrupt_'+name);shutil.copytree(original_root,copied)
            changed=copy.deepcopy(report)
            changed['artifact_paths']=[str(copied/Path(p).name) for p in report['artifact_paths']]
            for p in report['artifact_paths']:changed['source_hashes'].pop(p)
            path=copied/filename
            with gzip.open(path,'rt') as handle:lines=handle.readlines()
            if name=='coordinate_end':
                fields=lines[1].rstrip('\n').split('\t');fields[5]=str(int(fields[5])+1);lines[1]='\t'.join(fields)+'\n'
            else:
                rows=[json.loads(line) for line in lines];mutate(rows)
                lines=[json.dumps(row)+'\n' for row in rows]
            with gzip.open(path,'wt') as handle:handle.writelines(lines)
            for p in changed['artifact_paths']:bind(changed['source_hashes'],p)
            try:check_taxon(source,changed)
            except (AssertionError,StopIteration):checks.append('independent_reader_rejects_'+name+'_despite_rebound_hash')
            else:raise AssertionError('Corruption accepted: '+name)
    pins={}
    for name in ('readback_full_genome_annotation_cds_v1.py','audit_full_genome_annotation_cds_v1.py','genomic_cds_join_v1.py',
                 'build_full_annotation_coordinate_registry_v2.py','check_full_annotation_coordinate_registry_cases_v2.py',Path(__file__).name):
        bind(pins,Path('scripts')/name)
    result=dict(status='passed_offline_independent_genomic_CDS_full_source_replay_controls',checked_utc=datetime.now(timezone.utc).isoformat(),
        checks=checks,count=len(checks),source_hashes=pins,scientific_eligibility=False,corpus_pilot=False,gpu=False,new_predictions=0,
        scope='Hand IUPAC complement, all four source modes, duplicate reference/missing source/target retention and10rebound-hash '
              'semantic corruptions. Reader imports no producer decoder, grouping or complement implementation. '
              'Literal software evidence only; no full corpus or biological acceptance.')
    with args.receipt.open('x') as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
