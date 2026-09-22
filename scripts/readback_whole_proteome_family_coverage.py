#!/usr/bin/env python3
"""Independently reconstruct atlas family coverage using Python sets."""
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import sqlite3
import time


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def check(condition, message):
    if not condition:
        raise ValueError(message)


def connect(path):
    return sqlite3.connect(Path(path).resolve().as_uri() + '?mode=ro', uri=True)


def reconstruct(catalog, bridge, observed, coverage, guides):
    """Check all model identities and every family, without producer aggregation."""
    links = {}
    with Path(catalog).open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = (row['taxon_id'], row['protein_id'])
            check(key not in links, 'Duplicate catalog protein')
            links[key] = (*key, row['sequence_sha256'], row['model_id'],
                          int(row['version']), row['model_path'])
    count = len(links)
    expected = {}
    source_proteins = 0
    with connect(bridge) as source:
        for gene, taxon, protein, sequence in source.execute(
                'SELECT native_gene_id,taxon_id,protein_id,sequence_id FROM proteins'):
            source_proteins += 1
            model = links.pop((taxon, protein), None)
            if model is not None:
                check(sequence == 'S' + model[2], 'Source sequence mismatch')
                check(gene not in expected, 'Duplicate source gene')
                expected[gene] = model
        check(not links and len(expected) == count, 'Unmatched catalog proteins')
        del links
        with connect(observed) as result:
            seen = set()
            for row in result.execute('SELECT native_gene_id,taxon_id,protein_id,'
                                      'sequence_sha256,model_id,version,model_path FROM structures'):
                gene = row[0]
                check(gene not in seen, 'Duplicate output gene')
                check(expected.get(gene) == row[1:], 'Output model identity mismatch')
                seen.add(gene)
            check(len(seen) == count, 'Missing output model identities')
        del seen
        summaries = []
        with Path(coverage).open() as handle:
            table = csv.reader(handle, delimiter='\t')
            check(next(table) == ['guide', 'family', 'proteins', 'taxa', 'structure_proteins',
                                  'structure_taxa', 'structure_sequences', 'structure_models'],
                  'Coverage header mismatch')
            for guide in guides:
                summary = dict(guide=guide, families=0, proteins=0, structure_proteins=0,
                               families_with_models_in_multiple_taxa=0,
                               families_with_all_proteins_modeled=0)
                thresholds = {str(n): 0 for n in (2, 4, 10, 25, 50, 100)}
                # Source membership index bounds memory to one family; no output DB join.
                rows = source.execute('''SELECT a.family,a.native_gene_id,p.taxon_id
                    FROM assignments a INDEXED BY assignments_family
                    JOIN proteins p ON p.native_gene_id=a.native_gene_id
                    WHERE a.guide=? ORDER BY a.family''', (guide,))
                for family, members in itertools.groupby(rows, key=lambda row: row[0]):
                    genes, taxa, modeled_taxa, sequences, models = set(), set(), set(), set(), set()
                    covered = 0
                    for _, gene, taxon in members:
                        check(gene not in genes, 'Repeated family member')
                        genes.add(gene); taxa.add(taxon)
                        model = expected.get(gene)
                        if model is not None:
                            covered += 1; modeled_taxa.add(taxon)
                            sequences.add(model[2]); models.add(model[5])
                    values = [len(genes), len(taxa), covered, len(modeled_taxa), len(sequences), len(models)]
                    check(next(table, None) == [guide, family, *map(str, values)],
                          'Family coverage mismatch: ' + guide + '/' + family)
                    summary['families'] += 1
                    summary['proteins'] += len(genes)
                    summary['structure_proteins'] += covered
                    summary['families_with_models_in_multiple_taxa'] += len(modeled_taxa) >= 2
                    summary['families_with_all_proteins_modeled'] += covered == len(genes)
                    for threshold in thresholds:
                        thresholds[threshold] += len(modeled_taxa) >= int(threshold)
                check(summary['proteins'] == source_proteins and summary['structure_proteins'] == count,
                      'Incomplete family partition')
                summary['families_at_minimum_modeled_taxa'] = thresholds
                summaries.append(summary)
                print(json.dumps(summary), flush=True)
            check(next(table, None) is None, 'Extra coverage rows')
    return dict(protein_links=count, source_proteins=source_proteins, guides=summaries)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    started = time.time()
    plan = json.loads(args.plan.read_text()); plan_sha = sha(args.plan)
    output = Path(plan['output'])
    if output.exists():
        raise FileExistsError(output)
    def verify():
        check(sha(args.plan) == plan_sha, 'Changed readback plan')
        for path, digest in plan['pins'].items():
            check(sha(path) == digest, 'Changed pinned file: ' + path)
    verify()
    producer_plan = json.loads(Path(plan['producer_plan']).read_text())
    base = Path(producer_plan['output'])
    receipt = json.loads((base / 'receipt.json').read_text())
    check(receipt['status'] == 'complete_structure_family_coverage_bridge_pending_independent_readback',
          'Producer incomplete')
    check(receipt['plan_sha256'] == sha(plan['producer_plan']), 'Producer plan binding mismatch')
    for name, digest in receipt['artifacts'].items():
        check(sha(base / name) == digest, 'Changed producer artifact')
    catalog = Path(producer_plan['catalog'])
    check(receipt['catalog_receipt_sha256'] == sha(catalog / 'receipt.json'), 'Catalog binding mismatch')
    check(receipt['catalog_readback_sha256'] == sha(producer_plan['readback']), 'Catalog audit binding mismatch')
    result = reconstruct(catalog / 'protein_model_links.tsv', producer_plan['bridge'],
                         base / 'structure_family_bridge.sqlite', base / 'family_structure_coverage.tsv',
                         [g['guide'] for g in receipt['guides']])
    check(result['protein_links'] == receipt['protein_links'], 'Producer link total differs')
    for computed, reported in zip(result['guides'], receipt['guides']):
        check(all(computed[k] == v for k, v in reported.items()), 'Producer summary differs')
    verify()
    result.update(status='passed_independent_full_structure_family_coverage_readback',
                  plan_sha256=plan_sha, producer_receipt_sha256=sha(base / 'receipt.json'),
                  elapsed_seconds=time.time() - started,
                  scope='Every model identity and every family coverage row reconstructed from catalog and source membership using Python sets. Availability only; no confidence qualification, orthology, duplication events or statistical power validation.')
    with output.open('x') as handle:
        json.dump(result, handle, indent=2); handle.write('\n')


if __name__ == '__main__':
    main()
