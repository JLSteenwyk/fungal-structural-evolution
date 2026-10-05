#!/usr/bin/env python3
"""Offline literal annotation/alternative-product cases and independent corruption checks."""
import argparse
import copy
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import sqlite3
import tempfile

from build_full_annotation_coordinate_registry_v1 import build_taxon
from readback_full_annotation_coordinate_registry_v1 import check_taxon
from reference_measurement_union_sources import bind


def entry(root, mode):
    proteins = [('p,1', 'MKTA'), ('alt', 'MKT'), ('unmapped', 'MA')]
    columns = ['taxon_id', 'protein_id', 'gene_ids_json', 'status', 'partial_cds', 'protein_length']
    mapping, decisions = root / 'mapping.tsv', root / 'decisions.tsv'
    with mapping.open('w') as mh, decisions.open('w') as dh:
        mw = csv.DictWriter(mh, columns, delimiter='\t'); mw.writeheader()
        dw = csv.DictWriter(dh, columns + ['decision'], delimiter='\t'); dw.writeheader()
        for i, (ident, sequence) in enumerate(proteins):
            row = dict(taxon_id='literal', protein_id=ident, gene_ids_json=json.dumps(['g'] if i < 2 else []),
                       status='unique_gene' if i < 2 else 'unmapped', partial_cds='unknown', protein_length=len(sequence))
            mw.writerow(row)
            dw.writerow(dict(row, decision='alternative_product_retained_in_source' if i == 1 else
                             ('longest_per_gene_lexical_tiebreak' if i == 0 else 'retained_unresolved_gene')))
    proteome, selected = root / 'source.faa', root / 'selected.faa'
    proteome.write_text(''.join('>'+ident+' original\n'+sequence+'\n' for ident, sequence in proteins))
    selected.write_text(''.join('>'+ident+' original\n'+sequence+'\n' for i, (ident, sequence) in enumerate(proteins) if i != 1))
    annotation = root / 'annotation.tsv'
    if mode == 'ncbi_protein_gff':
        annotation.write_text('##gff-version 3\n'
          'chr\tR\tregion\t1\t100\t.\t+\t.\tID=r;Is_circular=true\n'
          'chr\tR\tgene\t1\t110\t.\t-\t.\tID=g\n'
          'chr\tR\tmRNA\t1\t110\t.\t-\t.\tID=t;Parent=g\n'
          'chr\tR\tCDS\t95\t110\t.\t-\t2\tID=c;Parent=t;protein_id=p%2C1;part=1/2;exception=ribosomal slippage;transl_except=(pos:98..100%2Caa:Sec)\n'
          'chr\tR\tCDS\t1\t7\t.\t-\t1\tID=c;Parent=t;protein_id=p%2C1;part=2/2;partial=true;start_range=.,1\n'
          'chr\tR\tCDS\t20\t22\t.\t+\t0\tParent=g;protein_id=outside\n')
    elif mode == 'transcript_gff':
        annotation.write_text('##gff-version 3\n'
          'chr\tEVM\tgene\t1\t20\t.\t+\t.\tID=g\n'
          'chr\tEVM\tmRNA\t1\t20\t.\t+\t.\tID=p%2C1;Parent=g\n'
          'chr\tEVM\tCDS\t2\t13\t.\t+\t0\tParent=p%2C1;gene=g\n')
    elif mode == 'creolimax_gtf':
        annotation.write_text('chr\tA\tgene\t1\t20\t.\t-\t.\tbareGene\n'
          'chr\tA\tCDS\t2\t13\t.\t-\t0\tgene_id "g"; transcript_id "p,1";\n')
    elif mode == 'verified_orf_coordinates':
        with annotation.open('w') as handle:
            w = csv.DictWriter(handle, ['protein_id', 'status', 'contig', 'start', 'end', 'strand', 'contig_name_match', 'orf_id'], delimiter='\t')
            w.writeheader()
            for i, (ident, sequence) in enumerate(proteins):
                w.writerow(dict(protein_id=ident, status='coordinates_out_of_bounds' if i == 2 else 'exact_translation',
                    contig='chr', start=1, end=12 if i != 2 else 999, strand='-', contig_name_match='exact', orf_id='ORF:'+ident))
    else:
        raise ValueError(mode)
    pins = {}
    for path in (mapping, decisions, proteome, selected, annotation):
        bind(pins, path)
    return dict(taxon_id='literal', mapping_mode=mode, source_products=3, selected_representatives=2,
                mapping_path=str(mapping), decisions_path=str(decisions), proteome_path=str(proteome),
                representative_path=str(selected), annotation_path=str(annotation), pins=pins)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    assert not args.receipt.exists()
    checks = []
    budget = dict(emergency_free_disk_gib=1, output_allowance_gib=1)
    with tempfile.TemporaryDirectory(prefix='literal-annotation-', dir='results') as folder:
        root = Path(folder)
        originals = {}
        for mode in ('ncbi_protein_gff', 'transcript_gff', 'creolimax_gtf', 'verified_orf_coordinates'):
            case = root / mode; case.mkdir()
            source = entry(case, mode)
            output = case / 'index'; output.mkdir()
            report = build_taxon(source, output, budget)
            replay = check_taxon(source, report)
            assert replay['source_products'] == 3 and replay['selected_representatives'] == 2
            checks.append(mode + ':full_literal_independent_replay')
            originals[mode] = (source, report)
        source, report = originals['ncbi_protein_gff']
        db = sqlite3.connect(report['database'])
        assert db.execute("SELECT start,end,phase,feature_id FROM features WHERE feature_type='CDS' ORDER BY row_id").fetchall() == [(95,110,'2','c'),(1,7,'1','c'),(20,22,'0',None)]
        attrs = json.loads(db.execute('SELECT attributes_json FROM features WHERE row_id=4').fetchone()[0])
        assert attrs['protein_id'] == ['p,1'] and attrs['transl_except'] == ['(pos:98..100,aa:Sec)']
        assert report['candidate_link_counts']['explicit_protein_id:outside_source_proteome'] == 1
        assert report['products_without_direct_cds_candidate'] == 2
        assert db.execute('SELECT selected_representative FROM products WHERE protein_id="alt"').fetchone() == (0,)
        db.close()
        checks.extend(['literal_multipart_phase_circular_virtual_coordinates_retained',
                       'literal_escaped_commas_and_translation_exception_retained',
                       'literal_alternative_unresolved_and_outside_source_products_retained'])
        mutations = {
            'coordinate': 'UPDATE features SET end=100 WHERE row_id=4',
            'phase': "UPDATE features SET phase='0' WHERE row_id=4",
            'multipart_identifier': "UPDATE features SET feature_id='replacement' WHERE row_id=5",
            'candidate_reference': "UPDATE cds_product_refs SET product_id='alt' WHERE feature_row_id=4",
            'parent': "UPDATE feature_parents SET parent_id='wrong' WHERE feature_row_id=4",
            'representative_selection': "UPDATE products SET selected_representative=1 WHERE protein_id='alt'",
            'gene_mapping': "UPDATE products SET mapping_json='{}' WHERE protein_id='p,1'",
            'missing_feature': 'DELETE FROM features WHERE row_id=5',
            'lost_exception': "UPDATE features SET attributes_json='{}' WHERE row_id=4"}
        for name, command in mutations.items():
            corrupt = root / ('corrupt_' + name + '.sqlite')
            shutil.copyfile(report['database'], corrupt)
            con = sqlite3.connect(corrupt); con.execute(command); con.commit(); con.close()
            changed = copy.deepcopy(report)
            changed['database'] = str(corrupt)
            changed['database_bytes'] = corrupt.stat().st_size
            changed['source_hashes'].pop(report['database'])
            bind(changed['source_hashes'], corrupt)
            try:
                check_taxon(source, changed)
            except (AssertionError, StopIteration):
                checks.append('independent_reader_rejects_' + name + '_despite_rebound_hash')
            else:
                raise AssertionError('Corruption admitted: ' + name)
    pins = {}
    for name in ('build_full_annotation_coordinate_registry_v1.py',
                 'readback_full_annotation_coordinate_registry_v1.py', Path(__file__).name):
        bind(pins, Path('scripts') / name)
    result = dict(status='passed_offline_literal_full_annotation_registry_controls',
        checked_utc=datetime.now(timezone.utc).isoformat(), checks=checks, count=len(checks),
        source_hashes=pins, scientific_eligibility=False, corpus_pilot=False, gpu=False,
        scope='All four source modes, complete literal source/product inventory and deliberately corrupted '
              'coordinate/phase/parent/multipart/product/exception/baseline records. No corpus or genome reconstruction qualification.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2); handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
