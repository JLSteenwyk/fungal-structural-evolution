#!/usr/bin/env python3
"""Exercise complete synthetic source joins and source-rebound semantic corruptions."""
import argparse
import csv
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
import sqlite3

from coding_structure_source_coupling_v1 import couple_taxon
from readback_full_coding_structure_source_coupling_v1 import check_taxon
from reference_measurement_union_sources import bind


def write_rows(path, rows):
    with gzip.open(path, 'wt') as handle:
        for row in rows:
            handle.write(json.dumps(row) + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    args = p.parse_args()
    assert not args.receipt.exists()
    args.output.mkdir(exist_ok=False)
    target = dict(ordinal=1, cds_id='dna1', product_id='p1', description='dna1 original header',
        sequence_sha256='d' * 64, sequence_length=9, status='one_exact_unmodified_genomic_cds_candidate',
        candidate_count=1, matching_candidate_indices=[0])
    extra = dict(target, ordinal=2, cds_id='dna2', product_id=None, description='dna2 unlinked header',
                 status='missing_or_ambiguous_target_product_id', candidate_count=0, matching_candidate_indices=[])
    products = [dict(protein_id='p1', protein_length=3, sequence_sha256='a' * 64, selected_representative=True,
        original_mapping={'status': 'unresolved_gene'}, original_decision={'decision': 'baseline_only'}, genomic_candidate_count=1,
        original_cds_targets=[dict(cds_id='dna1', ordinal=1, status=target['status'])], missing_target_reason=None),
        dict(protein_id='p2', protein_length=3, sequence_sha256='b' * 64, selected_representative=False,
        original_mapping={'status': 'unresolved_gene'}, original_decision={'decision': 'alternative'}, genomic_candidate_count=0,
        original_cds_targets=[], missing_target_reason=None)]
    original_targets = args.output / 'original_targets.jsonl.gz'
    original_products = args.output / 'original_products.jsonl.gz'
    write_rows(original_targets, [target, extra])
    write_rows(original_products, products)
    translation = args.output / 'translation.tsv'
    columns = ['taxon_id','cds_id','protein_id','status','translation_table','code_source','terminal_stop','gff_flags','cds_header']
    with translation.open('w') as handle:
        w = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerow(dict(taxon_id='literal', cds_id='dna1', protein_id='p1', status='exact_translation',
            translation_table='12', code_source='explicit_gff_transl_table', terminal_stop='False', gff_flags='', cds_header=target['description']))
        w.writerow(dict(taxon_id='literal', cds_id='dna2', protein_id='', status='missing_or_ambiguous_protein_id',
                        cds_header=extra['description']))
    database = args.output / 'source_models.sqlite'
    con = sqlite3.connect(database)
    con.execute('CREATE TABLE proteins(taxon_id TEXT,protein_id TEXT,sequence_sha256 TEXT,length INTEGER,has_afdb INTEGER,has_esmfold INTEGER)')
    con.execute('INSERT INTO proteins VALUES(?,?,?,?,?,?)', ('literal','p1','a'*64,3,1,1))
    con.commit()
    con.close()
    models = args.output / 'source_models.tsv'
    header = 'taxon_id\tprotein_id\tsequence_sha256\tlength\tavailability\tafdb_model_id\tesmfold_model_id\n'
    data = 'literal\tp1\t' + 'a'*64 + '\t3\tboth\tAF_LITERAL\tESM_LITERAL\n'
    models.write_text(header + data)
    span = dict(start=len(header.encode()), end=models.stat().st_size, records=1)
    entry = dict(taxon_id='literal', mapping_mode='ncbi_protein_gff', translation_table=str(translation),
        genomic_targets=str(original_targets), genomic_products=str(original_products), study_role='ingroup',
        species_name='literal software source', lineage='literal', target_records=2, source_products=2, selected_representatives=1)
    produced = args.output / 'coupled'
    produced.mkdir()
    report = couple_taxon(entry, database, produced)
    proof = check_taxon(entry, report, models, header, span)
    assert proof['target_records'] == proof['source_products'] == 2 and proof['selected_representatives'] == 1
    changes = []
    for filename, field, value in [('target_evidence.jsonl.gz','joined_status','incorrect_admission'),
                                    ('product_evidence.jsonl.gz','original_mapping',{'status':'invented_gene'})]:
        path = produced / 'literal' / filename
        original = path.read_bytes()
        with gzip.open(path, 'rt') as handle:
            rows = [json.loads(line) for line in handle]
        rows[0][field] = value
        write_rows(path, rows)
        try:
            check_taxon(entry, report, models, header, span)
        except AssertionError:
            changes.append('reject_self_consistent_hash_' + field)
        else:
            raise AssertionError(field)
        path.write_bytes(original)
    altered_models = args.output / 'changed_models.tsv'
    altered_models.write_text(header + data.replace('both', 'neither'))
    try:
        check_taxon(entry, report, altered_models, header, span)
    except AssertionError:
        changes.append('reject_model_state_inconsistent_with_source_ids')
    else:
        raise AssertionError('model state')
    pins = {}
    for path in [Path(__file__), Path('scripts/readback_full_coding_structure_source_coupling_v1.py'),
                 Path('scripts/coding_structure_source_coupling_v1.py'), *[q for q in args.output.rglob('*') if q.is_file()]]:
        bind(pins, path)
    result = dict(status='passed_separate_coding_structure_full_reader_literal_source_integration',
        checked_utc=datetime.now(timezone.utc).isoformat(), literal_target_records=2, literal_products=2,
        selected_representatives=1, preserved_unlinked_target_records=1, preserved_alternatives=1,
        separate_model_TSV_vs_SQLite=True, rejected_semantic_changes=changes, source_hashes=pins,
        scientific_eligibility=False, gpu=False, new_predictions=0,
        scope='Literal synthetic complete source join plus three semantic source/output corruptions; every real526entry '
              'corpus record still requires full replay. No real taxon pilot, sequence translation, model inference '
              'or biological admission. CSV/JSON/gzip libraries and original evidence remain shared dependencies.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
