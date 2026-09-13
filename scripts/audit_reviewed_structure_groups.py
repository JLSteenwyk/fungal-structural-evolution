#!/usr/bin/env python3
"""Read back derived membership against raw exact-conversion scores and provenance."""
import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha, read_table


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['groups', 'clusters', 'annotations', 'output']:
        p.add_argument('--' + name, required=True, type=Path)
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument('--exact', type=Path)
    source.add_argument('--controls', type=Path)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use immutable output')
    receipt = checked_receipt(a.groups)
    for name in ['clusters', 'annotations']:
        checked_receipt(getattr(a, name))
        assert receipt['source_receipts'][name] == sha(getattr(a, name) / 'receipt.json')
    if a.controls:
        control = checked_receipt(a.controls)
        config = json.loads((a.controls / 'config.json').read_text())
        assert control['config_sha256'] == sha(a.controls / 'config.json')
        command = config['commands'][1]
        score_path = a.controls / 'inclusive-span.tsv'
        score_receipt_path = a.controls / 'receipt.json'
        assert Path(command[5]).resolve() == score_path.resolve()
        assert Path(command[0]).name == 'foldseek-inclusive-span'
    else:
        config = json.loads((a.exact / 'config.json').read_text())
        completion = json.loads((a.exact / 'completion.json').read_text())
        assert completion['alignment_table_sha256'] == sha(a.exact / 'alignments.tsv')
        assert completion['config_sha256'] == sha(a.exact / 'config.json')
        command = config['command']
        score_path = a.exact / 'alignments.tsv'
        score_receipt_path = a.exact / 'completion.json'
    assert command[command.index('--exact-tmscore') + 1] == '1'
    fields = command[command.index('--format-output') + 1].split(',')
    scores = {}
    for values in csv.reader(score_path.open(), delimiter='\t'):
        assert len(values) == len(fields)
        row = dict(zip(fields, values)); key = row['query'], row['target']
        assert key not in scores
        numeric = {k: float(row[k]) for k in fields[2:]}
        bounded = ['qcov', 'tcov', 'alntmscore', 'qtmscore', 'ttmscore', 'lddt']
        valid = all(math.isfinite(v) and v >= 0 for v in numeric.values()) and all(numeric[k] <= 1 for k in bounded)
        scores[key] = valid and numeric['evalue'] <= 1e-5 and min(numeric['qcov'], numeric['tcov']) >= .8 and numeric['alntmscore'] >= .5
    original = read_table(a.clusters / 'model_cluster_membership.tsv')
    expected_pairs = {(r['representative_input_id'], r['cluster_input_id']) for r in original if r['representative_input_id'] != r['cluster_input_id']}
    assert set(scores) == expected_pairs | {(b, a) for a, b in expected_pairs}
    retained = read_table(a.groups / 'retained_membership.tsv')
    deferred = read_table(a.groups / 'deferred_members.tsv')
    by_id = {r['cluster_input_id']: r for r in retained + deferred}
    retained_ids = {r['cluster_input_id'] for r in retained}
    assert len(by_id) == len(retained) + len(deferred) == len(original)
    assert set(by_id) == {r['cluster_input_id'] for r in original}
    for row in original:
        rep, model = row['representative_input_id'], row['cluster_input_id']
        derived = by_id[model]
        assert all(derived[k] == v for k, v in row.items())
        passed = rep == model or (scores[rep, model] and scores[model, rep])
        assert (model in retained_ids) == passed
        assert derived['reviewed_group_id'] == ('RG_' + rep if passed else '')
    links = read_table(a.groups / 'model_taxon_marker_links.tsv')
    oldlinks = read_table(a.annotations / 'cluster_taxon_marker_links.tsv')
    keys = list(oldlinks[0])
    assert Counter(tuple(r[k] for k in keys) for r in links) == Counter(tuple(r[k] for k in keys) for r in oldlinks)
    for row in links:
        assert row['reviewed_group_id'] == by_id[row['cluster_input_id']]['reviewed_group_id']
        assert row['disposition'] == by_id[row['cluster_input_id']]['disposition']
    summaries = read_table(a.groups / 'group_summary.tsv')
    counts = Counter(r['reviewed_group_id'] for r in retained)
    assert {r['reviewed_group_id']: int(r['models']) for r in summaries} == counts
    a.output.mkdir(parents=True)
    result = {'status': 'passed_complete_membership_score_and_provenance_readback',
              'group_receipt_sha256': sha(a.groups / 'receipt.json'),
              'score_source_receipt_sha256': sha(score_receipt_path),
              'script_sha256': sha(Path(__file__)), 'models_checked': len(original),
              'directed_scores_checked': len(scores), 'taxon_marker_links_checked': len(links),
              'groups_checked': len(summaries),
              'interpretation': 'Membership eligibility recalculated from raw numeric scores, independently of edge-review booleans; complete source fields, partitions, group model counts and annotation multiset checked. Not an independent alignment or TM-score implementation.'}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
