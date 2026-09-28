"""Compare both full family-coverage tables after independent readback.

Resources: single CPU, streaming rows, under 1 GiB memory, up to 1 GiB output;
planning 1–15 minutes. Catalog availability is not confidence qualification.
"""
import csv
import itertools
import json
from collections import Counter, defaultdict
from pathlib import Path
from ancestral_chain_attempt import write_json
from readback_whole_proteome_catalog import sha


def records(path):
    with path.open() as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def main():
    roots = [Path(f'results/structures/whole-proteome-family-coverage-{date}-v1')
             for date in ['20260922', '20260928']]
    sources, receipts = {}, []
    for root in roots:
        receipt_path = root / 'receipt.json'
        receipt = json.loads(receipt_path.read_text())
        audit_path = Path(str(root) + '-readback.json')
        audit = json.loads(audit_path.read_text())
        assert audit['status'] == 'passed_independent_full_structure_family_coverage_readback'
        assert audit['producer_receipt_sha256'] == sha(receipt_path)
        table = root / 'family_structure_coverage.tsv'
        assert sha(table) == receipt['artifacts'][table.name]
        for path in [receipt_path, audit_path, table]:
            sources[str(path)] = sha(path)
        receipts.append(receipt)
    out = Path('results/structures/whole-proteome-family-coverage-change-20260928-v1')
    out.mkdir(parents=True, exist_ok=False)
    fields = ['structure_proteins', 'structure_taxa', 'structure_sequences', 'structure_models']
    columns = ['guide', 'family', 'proteins', 'taxa'] + [f + '_' + s for f in fields for s in ['old', 'new', 'change']]
    columns += ['newly_multi_taxon', 'lost_multi_taxon']
    counts = defaultdict(Counter)
    table = out / 'family_coverage_change.tsv'
    previous = {}
    with table.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        for a, b in itertools.zip_longest(*(records(r / 'family_structure_coverage.tsv') for r in roots)):
            assert a is not None and b is not None
            assert all(a[k] == b[k] for k in ['guide', 'family', 'proteins', 'taxa'])
            guide = a['guide']; family = a['family']
            assert guide not in previous or family > previous[guide]
            previous[guide] = family
            row = {k: a[k] for k in ['guide', 'family', 'proteins', 'taxa']}
            c = counts[guide]; c['families'] += 1
            c['proteins'] += int(a['proteins'])
            for f in fields:
                old, new = int(a[f]), int(b[f])
                row.update({f + '_old': old, f + '_new': new, f + '_change': new - old})
                c[f + '_old'] += old; c[f + '_new'] += new
            old, new = int(a['structure_taxa']), int(b['structure_taxa'])
            row['newly_multi_taxon'] = old < 2 <= new
            row['lost_multi_taxon'] = new < 2 <= old
            c['newly_multi_taxon'] += row['newly_multi_taxon']
            c['lost_multi_taxon'] += row['lost_multi_taxon']
            c['newly_any_model'] += int(a['structure_proteins']) == 0 < int(b['structure_proteins'])
            c['lost_all_models'] += int(b['structure_proteins']) == 0 < int(a['structure_proteins'])
            for minimum in [2, 4, 10, 25, 50, 100]:
                c[f'old_at_least_{minimum}_taxa'] += old >= minimum
                c[f'new_at_least_{minimum}_taxa'] += new >= minimum
            writer.writerow(row)
    for receipt, side in zip(receipts, ['old', 'new']):
        for summary in receipt['guides']:
            c = counts[summary['guide']]
            assert c['families'] == summary['families'] and c['proteins'] == summary['proteins']
            assert c['structure_proteins_' + side] == summary['structure_proteins']
            assert c[side + '_at_least_2_taxa'] == summary['families_with_models_in_multiple_taxa']
    # Full serialized-table arithmetic readback, independent of writer records.
    replay = defaultdict(Counter)
    for row in records(table):
        c = replay[row['guide']]; c['families'] += 1
        for f in fields:
            assert int(row[f + '_new']) - int(row[f + '_old']) == int(row[f + '_change'])
        old, new = int(row['structure_taxa_old']), int(row['structure_taxa_new'])
        for field, truth in [('newly_multi_taxon', old < 2 <= new), ('lost_multi_taxon', new < 2 <= old)]:
            assert row[field] == str(truth); c[field] += truth
    for guide, c in replay.items():
        assert all(value == counts[guide][k] for k, value in c.items())
    result = dict(status='complete_full_family_coverage_change_with_arithmetic_readback',
                  sources=sources, script_sha256=sha(__file__), guides=dict(counts),
                  artifacts={table.name: sha(table)},
                  scope='Every family in each unchanged partition compared across catalogs. Thresholds describe model-bearing taxa, not independent contrasts or statistical power. Sequence/model counts summed across families can repeat entities and are not unique atlas counts. No structural divergence or orthology inference.')
    write_json(out / 'receipt.json', result)
    write_json(Path('metadata/whole_proteome_family_coverage_change_20260928.json'), result)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
