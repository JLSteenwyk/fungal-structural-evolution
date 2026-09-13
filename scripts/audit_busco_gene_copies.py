#!/usr/bin/env python3
"""Link raw BUSCO copy calls to explicit gene annotations without recalling BUSCO."""
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read_table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def classify(hits, decisions):
    """Preserve annotation uncertainty and distinguish loci from protein products."""
    ids = [r[2] for r in hits if r[1] != 'Missing']
    if len(set(ids)) != len(ids):
        raise ValueError('Repeated BUSCO protein identifier within marker')
    rows = [decisions[pid] for pid in ids]
    genes = set()
    unresolved = 0
    retained = 0
    for row in rows:
        annotated = json.loads(row['gene_ids_json'])
        if row['status'] == 'unique_gene':
            if len(annotated) != 1:
                raise ValueError('Inconsistent unique gene annotation')
            genes.add(annotated[0])
        else:
            unresolved += 1
        decision = row['decision']
        if decision not in {'longest_per_gene_lexical_tiebreak', 'retained_unresolved_gene',
                            'alternative_product_retained_in_source'}:
            raise ValueError('Unknown representative decision')
        retained += decision != 'alternative_product_retained_in_source'
    status = {r[1] for r in hits}
    if len(status) != 1 or not status <= {'Complete', 'Duplicated', 'Fragmented', 'Missing'}:
        raise ValueError('Inconsistent BUSCO status')
    status = next(iter(status))
    if (status == 'Duplicated' and len(ids) < 2) or (status in {'Complete', 'Fragmented'} and len(ids) != 1) or (status == 'Missing' and ids):
        raise ValueError('BUSCO status/count mismatch')
    category = 'not_duplicated'
    if status == 'Duplicated':
        category = ('unresolved_gene_mapping' if unresolved else
                    'one_annotated_gene_multiple_products' if len(genes) == 1 else
                    'multiple_annotated_genes')
    return {'busco_status': status, 'raw_hit_proteins': len(ids),
            'distinct_resolved_gene_ids': len(genes), 'unresolved_hit_proteins': unresolved,
            'hit_proteins_retained_as_representatives': retained,
            'duplicate_annotation_category': category,
            'protein_ids_json': json.dumps(ids), 'resolved_gene_ids_json': json.dumps(sorted(genes))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new immutable output directory')
    manifest_path = ROOT / 'metadata/analysis_manifest.tsv'
    manifest = read_table(manifest_path)
    names = {r['taxon_id'] for r in manifest}
    if len(names) != len(manifest):
        raise ValueError('Duplicate sampling taxa')
    marker_path = ROOT / 'results/phylogeny/markers-full-v1/receipt.json'
    marker = json.loads(marker_path.read_text())
    rep_path = ROOT / 'results/gene_representatives/full-v2/receipt.json'
    rep = json.loads(rep_path.read_text())
    if marker['status'] != 'complete_extraction' or rep['status'] != 'complete':
        raise ValueError('Complete upstream snapshots required')
    if marker['manifest_sha256'] != sha(manifest_path):
        raise ValueError('Sampling manifest changed')
    tables = {r['taxon_id']: r for r in marker['source_tables']}
    representatives = {r['taxon_id']: r for r in rep['taxa']}
    if set(tables) != names or set(representatives) != names:
        raise ValueError('Upstream snapshots differ from full sampling')
    details, summaries, sources = [], [], []
    marker_set = None
    for taxon in manifest:
        name = taxon['taxon_id']
        table_path = ROOT / tables[name]['path']
        decision_path = (ROOT / representatives[name]['path']).with_suffix('.decisions.tsv')
        if sha(table_path) != tables[name]['sha256'] or sha(decision_path) != representatives[name]['decisions_sha256']:
            raise ValueError(f'Changed upstream artifact: {name}')
        decision_rows = read_table(decision_path)
        decisions = {r['protein_id']: r for r in decision_rows}
        if len(decisions) != len(decision_rows) or len(decisions) != representatives[name]['source_proteins']:
            raise ValueError('Representative source count/identifier mismatch')
        if any(r['taxon_id'] != name for r in decision_rows):
            raise ValueError('Wrong taxon in representative decisions')
        hits = defaultdict(list)
        for line in table_path.read_text().splitlines():
            if line and not line.startswith('#'):
                row = line.split('\t')
                hits[row[0]].append(row)
        if marker_set is None:
            marker_set = set(hits)
        if len(hits) != marker['markers'] or set(hits) != marker_set:
            raise ValueError('Marker universe mismatch')
        rows = [dict(taxon_id=name, marker=m, **classify(hits[m], decisions)) for m in sorted(hits)]
        details.extend(rows)
        calls = Counter(r['busco_status'] for r in rows)
        categories = Counter(r['duplicate_annotation_category'] for r in rows)
        summaries.append({'taxon_id': name, 'species_name': taxon['species_name'],
            'study_role': taxon['study_role'], 'lineage_group': taxon['lineage'].split(';')[0],
            'markers': len(rows), 'raw_complete_single': calls['Complete'],
            'raw_duplicated': calls['Duplicated'], 'raw_fragmented': calls['Fragmented'], 'raw_missing': calls['Missing'],
            'duplicated_one_annotated_gene': categories['one_annotated_gene_multiple_products'],
            'duplicated_multiple_annotated_genes': categories['multiple_annotated_genes'],
            'duplicated_unresolved_gene_mapping': categories['unresolved_gene_mapping'],
            'detected_markers_with_no_original_hit_retained': sum(r['raw_hit_proteins'] > 0 and r['hit_proteins_retained_as_representatives'] == 0 for r in rows),
            'duplicated_markers_with_one_original_hit_retained': sum(r['busco_status'] == 'Duplicated' and r['hit_proteins_retained_as_representatives'] == 1 for r in rows)})
        sources.append({'taxon_id': name, 'busco_table_sha256': sha(table_path), 'representative_decisions_sha256': sha(decision_path)})
    args.output.mkdir(parents=True)
    for filename, rows in [('marker_gene_copies.tsv', details), ('taxon_gene_copy_summary.tsv', summaries)]:
        with (args.output / filename).open('w') as handle:
            writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader(); writer.writerows(rows)
    totals = {key: sum(r[key] for r in summaries) for key in summaries[0] if isinstance(summaries[0][key], int)}
    result = {'status': 'complete_raw_busco_gene_annotation_audit', 'taxa': len(summaries),
        'marker_taxon_cells': len(details), 'totals': totals,
        'manifest_sha256': sha(manifest_path), 'marker_extraction_receipt_sha256': sha(marker_path),
        'representative_receipt_sha256': sha(rep_path), 'script_sha256': sha(Path(__file__)),
        'sources': sources,
        'interpretation': 'Annotations partition raw protein-mode BUSCO hits into explicit genes and products. This is not a BUSCO rerun on representatives: newly detectable or threshold-changing hits are not inferred. Distinct annotated genes do not prove biological duplication, ploidy, haplotig redundancy or contamination; unresolved mappings remain separate. No marker or taxon is reclassified in production inputs.',
        'artifacts': {p.name: sha(p) for p in args.output.iterdir()}}
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'sources'}, indent=2))


if __name__ == '__main__':
    main()
