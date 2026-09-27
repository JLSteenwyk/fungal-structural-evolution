#!/usr/bin/env python3
"""Join the full terminal sister inventories by gene pair, preserving guide sensitivity."""
import argparse
import csv
import json
import shutil
import sqlite3
import time
from collections import Counter
from pathlib import Path
import psutil
from run_ortholog_pair_guide_comparison import sha

COPIED = ['family', 'gene_node', 'candidate_status', 'model_coverage',
          'parent_reported_duplication', 'sequence_tip_a', 'sequence_tip_b',
          'sequence_pair_distance', 'model_a', 'version_a', 'model_b', 'version_b']
CANDIDATE = 'cross_taxon_unreported_candidate'


def compare_pair(a, b):
    assert a is not None or b is not None
    exemplar = a if a is not None else b
    genes = (exemplar['gene_a'], exemplar['gene_b'])
    assert genes[0] < genes[1]
    for row in [a, b]:
        if row is not None:
            assert (row['gene_a'], row['gene_b']) == genes
    result = dict(gene_a=genes[0], gene_b=genes[1],
                  taxon_a=exemplar['taxon_a'], taxon_b=exemplar['taxon_b'],
                  presence='both' if a is not None and b is not None else 'profile_only' if a is not None else 'mafft_only')
    for guide, row in [('profile', a), ('mafft', b)]:
        for field in COPIED:
            result[guide + '_' + field] = row[field] if row is not None else ''
    if a is not None and b is not None:
        for field in ['taxon_a', 'taxon_b', 'model_coverage', 'model_a', 'version_a', 'model_b', 'version_b']:
            assert a[field] == b[field], 'Frozen model/taxon join differs: ' + field
        result['candidate_class_agrees'] = int(a['candidate_status'] == b['candidate_status'])
        result['family_label_agrees'] = int(a['family'] == b['family'])
    else:
        result['candidate_class_agrees'] = ''
        result['family_label_agrees'] = ''
    def modeled_candidate(row):
        return row is not None and row['candidate_status'] == CANDIDATE and row['model_coverage'] in ['two_distinct_models', 'identical_model']
    result['modeled_candidate_in_either_guide'] = int(modeled_candidate(a) or modeled_candidate(b))
    result['modeled_candidate_in_both_guides'] = int(modeled_candidate(a) and modeled_candidate(b))
    result['both_parents_unreported'] = int(a is not None and b is not None and a['parent_reported_duplication'] == b['parent_reported_duplication'] == '0')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args(); plan = json.loads(args.plan.read_text()); ph = sha(args.plan)
    def verify():
        assert sha(args.plan) == ph
        for path, digest in plan['pins'].items():
            assert sha(path) == digest, path
    verify()
    dep = plan['dependency']
    while True:
        try:
            proc = psutil.Process(dep['pid'])
            if proc.create_time() != dep['created'] or proc.status() == psutil.STATUS_ZOMBIE:
                break
            assert proc.cmdline() == dep['cmdline']
        except psutil.NoSuchProcess:
            break
        time.sleep(30)
    verify()
    audit = json.loads(Path(plan['readback']).read_text()); ah = sha(plan['readback'])
    root = Path(plan['source']); receipt = json.loads((root / 'receipt.json').read_text())
    assert audit['status'] == 'passed_full_terminal_sister_inventory_readback'
    assert audit['producer_receipt_sha256'] == sha(root / 'receipt.json')
    assert audit['guides'] == receipt['guides']
    for name, digest in receipt['artifacts'].items():
        assert sha(root / name) == digest
    out = Path(plan['output']); assert shutil.disk_usage(out.parent).free >= plan['resources']['minimum_free_disk_bytes']
    out.mkdir(exist_ok=False)
    database = out / 'pair_sources.sqlite'; con = sqlite3.connect(database)
    con.execute('PRAGMA cache_size=-65536')
    con.execute('PRAGMA temp_store=FILE')
    for guide in ['profile', 'mafft']:
        con.execute(f'CREATE TABLE {guide} (a TEXT, b TEXT, record TEXT NOT NULL, PRIMARY KEY(a,b)) WITHOUT ROWID')
        count = 0; batch = []
        with (root / (guide + '_terminal_sisters.tsv')).open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                assert row['guide'] == guide and row['gene_a'] < row['gene_b']
                batch.append((row['gene_a'], row['gene_b'], json.dumps(row, separators=(',', ':'))))
                count += 1
                if len(batch) == 10000:
                    con.executemany(f'INSERT INTO {guide} VALUES (?,?,?)', batch); con.commit(); batch.clear()
            con.executemany(f'INSERT INTO {guide} VALUES (?,?,?)', batch); con.commit()
        assert count == next(r['pairs'] for r in receipt['guides'] if r['guide'] == guide)
        print(guide, count, 'rows loaded', flush=True)
    query = '''SELECT p.a,p.b,p.record,m.record FROM profile p LEFT JOIN mafft m ON p.a=m.a AND p.b=m.b
               UNION ALL SELECT m.a,m.b,NULL,m.record FROM mafft m LEFT JOIN profile p ON p.a=m.a AND p.b=m.b WHERE p.a IS NULL
               ORDER BY 1,2'''
    counts = Counter(); candidate_counts = Counter(); classes = Counter(); total = 0; selected = 0
    with (out / 'guide_comparison.tsv').open('w') as full, (out / 'modeled_candidate_union.tsv').open('w') as subset:
        writer = candidate_writer = None
        for ga, gb, left, right in con.execute(query):
            a, b = json.loads(left) if left else None, json.loads(right) if right else None
            row = compare_pair(a, b)
            if writer is None:
                writer = csv.DictWriter(full, list(row), delimiter='\t', lineterminator='\n'); writer.writeheader()
                candidate_writer = csv.DictWriter(subset, list(row), delimiter='\t', lineterminator='\n'); candidate_writer.writeheader()
            writer.writerow(row); total += 1; counts[row['presence']] += 1
            classes[row['profile_candidate_status'] + '|' + row['mafft_candidate_status']] += 1
            if row['modeled_candidate_in_either_guide']:
                candidate_writer.writerow(row); selected += 1
                candidate_counts[row['presence']] += 1
                candidate_counts['both_candidate'] += row['modeled_candidate_in_both_guides']
                candidate_counts['both_candidate_parents_unreported'] += row['modeled_candidate_in_both_guides'] * row['both_parents_unreported']
            if total % 100000 == 0:
                print(total, 'union pairs exported', flush=True)
    for guide in ['profile', 'mafft']:
        assert counts['both'] + counts[guide + '_only'] == next(r['pairs'] for r in receipt['guides'] if r['guide'] == guide)
    con.close(); verify(); assert sha(plan['readback']) == ah
    for name, digest in receipt['artifacts'].items():
        assert sha(root / name) == digest
    result = dict(status='complete_terminal_sister_guide_comparison_pending_readback', plan_sha256=ph,
                  source_receipt_sha256=sha(root / 'receipt.json'), source_readback_sha256=ah,
                  union_pairs=total, presence=dict(counts), classification_cross_tab=dict(classes),
                  modeled_candidate_union=selected, modeled_candidate_counts=dict(candidate_counts),
                  artifacts={p.name: sha(p) for p in out.iterdir()},
                  scope='Full union keyed by lexical gene IDs, not guide-dependent node or family labels. Every source field retained in SQLite; all pairs exported and modeled cross-taxon unreported candidates separately flagged, including identical models. Shared guides are dependent alternatives. No orthology membership, matched controls or duplication-effect test established.')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n'); print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
