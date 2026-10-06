#!/usr/bin/env python3
"""Compare the complete focal orthogroup with other proteins from one annotation.

This is an annotation and repeat-risk audit. It does not decide whether a gene
is a transposable element, remove genes, or infer duplication history.
"""
import argparse
import csv
from collections import Counter, defaultdict
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import sqlite3
import statistics

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


FOCAL_SPECIES = '214'
FOCAL_TAXON = 'F181123'
REPEAT_WORDS = ('transpos', 'retrotrans', 'reverse transcriptase', 'integrase',
                'rve', 'gag', 'helitron', 'zinc knuckle', 'dede', 'rvt')


def quantiles(values):
    values = sorted(values)
    if not values:
        return dict(count=0, minimum=None, q25=None, median=None, q75=None, maximum=None)
    def pick(frac):
        return values[round((len(values) - 1) * frac)]
    return dict(count=len(values), minimum=values[0], q25=pick(.25), median=pick(.5),
                q75=pick(.75), maximum=values[-1])


def load_family(path):
    focal = {}
    with Path(path).open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['native_species_id'] == FOCAL_SPECIES:
                if row['gene_id'] in focal:
                    raise ValueError('Repeated focal family identifier')
                focal[row['gene_id']] = row
    if len(focal) != 46787:
        raise ValueError('Unexpected focal family count')
    return focal


def load_ids(path):
    identifiers = {}
    with Path(path).open() as handle:
        for line in handle:
            identifier, protein = line.rstrip('\n').split(': ', 1)
            if identifier.startswith(FOCAL_SPECIES + '_'):
                if identifier in identifiers or not protein:
                    raise ValueError('Repeated or invalid focal native identifier')
                identifiers[identifier] = protein
    if len(identifiers) != 133973:
        raise ValueError('Unexpected focal native-protein grid')
    return identifiers


def database_context(path):
    connection = sqlite3.connect('file:' + str(Path(path).resolve()) + '?mode=ro', uri=True)
    connection.row_factory = sqlite3.Row
    products = {}
    for row in connection.execute('SELECT protein_id, protein_length, sequence_sha256, mapping_json FROM products'):
        mapping = json.loads(row['mapping_json'])
        gene_ids = json.loads(mapping['gene_ids_json'])
        if len(gene_ids) != 1:
            raise ValueError('Unexpected nonunique gene mapping')
        products[row['protein_id']] = dict(protein_length=row['protein_length'],
                                           sequence_sha256=row['sequence_sha256'], gene_id=gene_ids[0])
    if len(products) != 133973:
        raise ValueError('Unexpected product count')
    genes = {}
    contig_lengths = {}
    for row in connection.execute("SELECT seqid, feature_type, start, end, feature_id, attributes_json FROM features WHERE feature_type IN ('gene','region')"):
        if row['feature_type'] == 'region':
            if row['seqid'] in contig_lengths:
                raise ValueError('Repeated region')
            contig_lengths[row['seqid']] = row['end'] - row['start'] + 1
        else:
            if row['feature_id'] in genes:
                raise ValueError('Repeated gene feature')
            genes[row['feature_id']] = dict(contig=row['seqid'], start=row['start'], end=row['end'],
                                             span=row['end'] - row['start'] + 1,
                                             attributes=json.loads(row['attributes_json']))
    if len(genes) != len(products):
        raise ValueError('Product/gene grid differs')
    cds = defaultdict(lambda: dict(segments=0, coding_bases=0, hypothetical=True))
    query = '''SELECT refs.product_id, features.start, features.end, features.attributes_json
               FROM cds_product_refs AS refs JOIN features ON features.row_id=refs.feature_row_id'''
    for row in connection.execute(query):
        item = cds[row['product_id']]
        item['segments'] += 1
        item['coding_bases'] += row['end'] - row['start'] + 1
        attributes = json.loads(row['attributes_json'])
        products_field = attributes.get('product', [])
        if products_field != ['hypothetical protein']:
            item['hypothetical'] = False
    if set(cds) != set(products):
        raise ValueError('CDS/product grid differs')
    for product, item in products.items():
        gene = genes.get(item['gene_id'])
        if gene is None or gene['contig'] not in contig_lengths:
            raise ValueError('Missing gene or contig context')
        item.update(gene, **cds[product], contig_length=contig_lengths[gene['contig']])
    ordered = defaultdict(list)
    for product, item in products.items():
        ordered[item['contig']].append((item['start'], item['end'], product))
    for contig, rows in ordered.items():
        rows.sort()
        for index, (start, end, product) in enumerate(rows):
            previous = rows[index - 1][1] if index else None
            following = rows[index + 1][0] if index + 1 < len(rows) else None
            item = products[product]
            item['scaffold_gene_count'] = len(rows)
            item['left_gap'] = None if previous is None else start - previous - 1
            item['right_gap'] = None if following is None else following - end - 1
            item['overlaps_neighbor'] = any(value is not None and value < 0
                                            for value in (item['left_gap'], item['right_gap']))
    return products, len(contig_lengths)


def load_uniprot(path):
    found = set()
    with Path(path).open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['taxon_id'] != FOCAL_TAXON:
                raise ValueError('Unexpected UniProt taxon')
            found.add(row['protein_id'])
    return found


def load_pfam(paths, hashes):
    hits = defaultdict(list)
    expected = set()
    for path in paths:
        label = Path(path).name.split('.')[0]
        if label in expected:
            raise ValueError('Repeated Pfam chunk label')
        expected.add(label)
        with gzip.open(path, 'rt') as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                identifier = row['sequence_id']
                if identifier.startswith('S') and identifier[1:] in hashes:
                    hits[identifier[1:]].append((row['pfam_accession'], row['pfam_name'], row['description']))
    if len(expected) != 64:
        raise ValueError('Expected all 64 full Pfam chunks')
    return hits


def cohort_summary(records):
    result = dict(proteins=len(records))
    for field in ('protein_length', 'span', 'coding_bases', 'segments', 'scaffold_gene_count', 'contig_length'):
        result[field] = quantiles([record[field] for record in records])
    result['hypothetical_products'] = sum(record['hypothetical'] for record in records)
    result['uniprot_exact_sequence_matches'] = sum(record['uniprot_match'] for record in records)
    result['pfam_detected_proteins'] = sum(bool(record['pfam_hits']) for record in records)
    result['repeat_keyword_pfam_proteins'] = sum(record['repeat_keyword_pfam'] for record in records)
    result['singleton_gene_scaffolds'] = sum(record['scaffold_gene_count'] == 1 for record in records)
    result['neighbor_overlaps'] = sum(record['overlaps_neighbor'] for record in records)
    gaps = [gap for record in records for gap in (record['left_gap'], record['right_gap']) if gap is not None]
    result['neighbor_gap'] = quantiles(gaps)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quality', required=True, type=Path)
    parser.add_argument('--sequence-ids', required=True, type=Path)
    parser.add_argument('--annotation-db', required=True, type=Path)
    parser.add_argument('--uniprot-matches', required=True, type=Path)
    parser.add_argument('--assembly-receipt', required=True, type=Path)
    parser.add_argument('--pfam', required=True, nargs='+', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--receipt', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists() or args.receipt.exists():
        raise FileExistsError('Fresh output and receipt paths are required')
    pins = {}
    for path in [args.quality, args.sequence_ids, args.annotation_db, args.uniprot_matches,
                 args.assembly_receipt, Path(__file__)]:
        bind(pins, path)
    for path in args.pfam:
        bind(pins, path)
    verify(pins)
    family = load_family(args.quality)
    native = load_ids(args.sequence_ids)
    if not set(family) <= set(native):
        raise ValueError('Family member missing native protein mapping')
    native_by_protein = {protein: identifier for identifier, protein in native.items()}
    if len(native_by_protein) != len(native):
        raise ValueError('Repeated focal source protein mapping')
    products, scaffolds = database_context(args.annotation_db)
    if set(native.values()) != set(products):
        raise ValueError('Native protein mapping differs from annotation products')
    family_proteins = {native[identifier] for identifier in family}
    hashes = {products[protein]['sequence_sha256'] for protein in products}
    hash_multiplicity = Counter(record['sequence_sha256'] for record in products.values())
    uniprot = load_uniprot(args.uniprot_matches)
    if not uniprot <= set(products):
        raise ValueError('UniProt protein outside source product grid')
    pfam = load_pfam(args.pfam, hashes)
    output = args.output
    output.mkdir(parents=True)
    rows = []
    pfam_counts = {name: Counter() for name in ('family', 'other')}
    for protein, record in products.items():
        identifier = native_by_protein[protein]
        group = 'family' if protein in family_proteins else 'other'
        source_hits = pfam.get(record['sequence_sha256'], [])
        repeat = any(any(word in ' '.join(hit).lower() for word in REPEAT_WORDS) for hit in source_hits)
        record.update(protein_id=protein, native_id=identifier, group=group,
                      uniprot_match=protein in uniprot, pfam_hits=source_hits,
                      repeat_keyword_pfam=repeat)
        for hit in source_hits:
            pfam_counts[group][hit] += 1
        rows.append(record)
    if sum(row['group'] == 'family' for row in rows) != 46787:
        raise ValueError('Family product count differs')
    if len(rows) != 133973:
        raise ValueError('Source product count differs')
    table = output / 'protein_context.tsv'
    fields = ['protein_id', 'native_id', 'group', 'sequence_sha256', 'protein_length', 'gene_id',
              'contig', 'contig_length', 'start', 'end', 'span', 'coding_bases', 'segments',
              'hypothetical', 'scaffold_gene_count', 'left_gap', 'right_gap', 'overlaps_neighbor',
              'uniprot_match', 'pfam_hit_count', 'repeat_keyword_pfam']
    with table.open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for row in sorted(rows, key=lambda item: (item['group'], item['protein_id'])):
            serializable = {field: row.get(field) for field in fields}
            serializable['pfam_hit_count'] = len(row['pfam_hits'])
            writer.writerow(serializable)
    top = {}
    for name, counts in pfam_counts.items():
        top[name] = [dict(pfam_accession=accession, pfam_name=pfam_name,
                          description=description, hit_rows=count)
                     for (accession, pfam_name, description), count in counts.most_common(100)]
    assembly = json.loads(args.assembly_receipt.read_text())
    result = dict(status='complete_og0000017_source_annotation_repeat_risk_audit',
        checked_utc=datetime.now(timezone.utc).isoformat(), family='OG0000017',
        focal_taxon=FOCAL_TAXON, focal_assembly='GCA_000469055.2',
        focal_assembly_metrics=assembly['all_assembly_metrics'], source_products=len(products),
        focal_family_products=len(family_proteins), other_focal_products=len(products)-len(family_proteins),
        unique_source_sequence_hashes=len(hashes),
        source_products_in_nonunique_sequence_groups=sum(count for count in hash_multiplicity.values() if count > 1),
        nonunique_source_sequence_groups=sum(count > 1 for count in hash_multiplicity.values()),
        source_scaffolds=scaffolds, cohorts=dict(family=cohort_summary([r for r in rows if r['group'] == 'family']),
                                                  other=cohort_summary([r for r in rows if r['group'] == 'other'])),
        top_pfam_hits=top, protein_context_table=str(table), protein_context_table_sha256=sha(table),
        source_hashes=pins, all_source_products_retained=True,
        family_members_excluded=0, repeat_assignment_inferred=False,
        duplication_inferred=False, homology_accepted=False, scientific_eligibility=False,
        scope='All133973 original selected source products from the focal assembly joined to every '
              'OG0000017 member through native IDs, original annotation coordinates, exact UniProt '
              'matches and all64 completed Pfam annotation chunks. Repeat-keyword annotations are '
              'descriptive flags, not TE classification. Fragmentation, gene density, short proteins, '
              'hypothetical labels and repeat-associated Pfam terms are annotation-risk evidence, not '
              'grounds for exclusion or biological duplication claims. No sequences, trees, memberships '
              'or reconciliation inputs are changed.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ('source_hashes', 'top_pfam_hits')}, indent=2))


if __name__ == '__main__':
    main()
