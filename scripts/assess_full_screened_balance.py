#!/usr/bin/env python3
"""Describe complete fixed-control balance, retention shifts and reuse after QC."""
import argparse
import csv
import gzip
import itertools
import json
import shutil
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
from full_screened_balance_statistics import FEATURES, feature_vector, balance
from full_screened_balance_sources import load, target_flags, KEY, MASKS, REUSE_FIELDS
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha

SHIFT_FIELDS = {'baseline_targets': 'baseline_targets', 'baseline_target_mean': 'target_mean', 'selection_mean_shift': 'mean_shift', 'selection_shift_in_baseline_sd': 'shift_in_baseline_sd', 'selection_shift_status': 'shift_status'}


def run(plan_path):
    plan_path = Path(plan_path); plan = json.loads(plan_path.read_text()); source, bindings = load(plan, plan_path)
    assert plan['features'] == FEATURES; nodes = source['nodes']; flags = target_flags(source, plan['policies'])
    indices = {kind: {ident: i for i, ident in enumerate(rows)} for kind, rows in nodes.items()}
    records = {kind: list(rows.values()) for kind, rows in nodes.items()}
    vectors = {kind: np.array([feature_vector(n) for n in rows]).reshape(-1, len(FEATURES)) for kind, rows in records.items()}
    baseline = {g: vectors['target'][[i for i, n in enumerate(records['target']) if n['guide'] == g]] for g in plan['guides']}
    eligible_baseline = {(g, m, screen['id']): vectors['target'][[i for i, n in enumerate(records['target']) if n['guide'] == g and flags[n['node_id']][j] & (1 << bit)]] for g in plan['guides'] for j, m in enumerate(MASKS) for bit, screen in enumerate(plan['screens'])}
    groups = defaultdict(list); background_flags = {}; selected = 0
    with gzip.open(source['root'] / 'selected_pair_coverage.tsv.gz', 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            tid, bid = row['target_id'], row['background_id']; t, b = nodes['target'][tid], nodes['background'][bid]
            key = row['guide'], row['policy'], row['scenario_id']; assert t['guide'] == b['guide'] == key[0]
            tv, bv, jv = [[int(row[k + '_' + m + '_pass_bits']) for m in MASKS] for k in ['target', 'control', 'joint']]
            assert tv == flags[tid] and all(0 <= v < 64 for v in bv) and bv[2] == bv[0] & bv[1] and jv == [a & b for a, b in zip(tv, bv)]
            assert bid not in background_flags or background_flags[bid] == bv; background_flags[bid] = bv
            if b['same_model']: assert bv == [0, 0, 0]
            groups[key].append((indices['target'][tid], indices['background'][bid], *jv)); selected += 1
    assert selected == plan['expected']['selected_records']
    out = Path(plan['output']); assert shutil.disk_usage(out.parent).free >= plan['resources']['minimum_free_disk_gib'] * 2**30; out.mkdir(exist_ok=False)
    coverage_rows = balance_rows = reuse_rows = retained_cells = strata = 0
    with (out / 'coverage.tsv').open('w') as cf, (out / 'balance.tsv').open('w') as bf, gzip.open(out / 'retained_control_reuse.tsv.gz', 'wt', compresslevel=1) as rf:
        cw = bw = None; rw = csv.DictWriter(rf, REUSE_FIELDS, delimiter='\t', lineterminator='\n'); rw.writeheader()
        for group in itertools.product(plan['guides'], plan['policies'], [s['scenario_id'] for s in source['scenarios']]):
            rows = np.array(groups.get(group, []), dtype=np.int64).reshape(-1, 5)
            assert len(set(rows[:, 0])) == len(rows)
            x, y = vectors['target'][rows[:, 0]], vectors['background'][rows[:, 1]]; strata += 1
            for mi, mask in enumerate(MASKS):
                for bit, screen in enumerate(plan['screens']):
                    key = (*group, mask, screen['id']); attr = source['attrition'][key]
                    keep = np.flatnonzero(rows[:, mi + 2] & (1 << bit)); rr = rows[keep]; xx, yy = x[keep], y[keep]; n = len(rr); retained_cells += n
                    assert n == int(attr['joint_pass_matched_records']) and len(rows) == int(attr['matched_records']) and len(baseline[group[0]]) == int(attr['all_target_records'])
                    eb = eligible_baseline[group[0], mask, screen['id']]; assert len(eb) == int(attr['target_pass_all_records'])
                    tn = [records['target'][i] for i in rr[:, 0]]; bn = [records['background'][i] for i in rr[:, 1]]
                    assert not any(v['same_model'] for v in tn + bn)
                    reuse = Counter(v['node_id'] for v in bn); physical = Counter(v['pair_key'] for v in bn)
                    node_weights = [1 / reuse[v['node_id']] for v in bn]; pair_weights = [1 / physical[v['pair_key']] for v in bn]
                    row = dict(zip(KEY, key)); row.update(all_target_records=len(baseline[group[0]]), metadata_matched_records=len(rows), metadata_unmatched_records=int(attr['unmatched_records']), target_eligible_all_records=len(eb), retained_matches=n, screen_excluded_matches=len(rows)-n, retained_target_taxa=len({v['taxon_id'] for v in tn}), retained_target_families=len({v['family'] for v in tn}), distinct_target_genes=len({v['gene_'+s] for v in tn for s in ['a','b']}), distinct_target_physical_pairs=len({v['pair_key'] for v in tn}), distinct_control_nodes=len(reuse), distinct_control_physical_pairs=len(physical), distinct_control_genes=len({v['gene_'+s] for v in bn for s in ['a','b']}), distinct_control_taxa=len({v['taxon_'+s] for v in bn for s in ['a','b']}), maximum_control_node_reuse=max(reuse.values(),default=0), maximum_control_physical_pair_reuse=max(physical.values(),default=0), top_five_control_node_fraction=sum(sorted(reuse.values(),reverse=True)[:5])/n if n else '', top_five_control_physical_pair_fraction=sum(sorted(physical.values(),reverse=True)[:5])/n if n else '', reciprocal_node_weight_sum=float(sum(node_weights)), reciprocal_physical_pair_weight_sum=float(sum(pair_weights)), node_weight_kish_concentration=(sum(node_weights)**2/sum(v*v for v in node_weights)) if n else '', physical_pair_weight_kish_concentration=(sum(pair_weights)**2/sum(v*v for v in pair_weights)) if n else '', zero_sequence_distance_targets=int(sum(xx[:,0]==0)), zero_sequence_distance_controls=int(sum(yy[:,0]==0)))
                    if cw is None: cw = csv.DictWriter(cf, list(row), delimiter='\t', lineterminator='\n'); cw.writeheader()
                    cw.writerow(row); coverage_rows += 1
                    for bid, count in sorted(reuse.items()):
                        pair = nodes['background'][bid]['pair_key']; rw.writerow(dict(zip(KEY,key),background_id=bid,background_pair_key=pair,retained_target_records=count,reciprocal_node_reuse_weight=1/count,retained_records_using_physical_pair=physical[pair],reciprocal_physical_pair_reuse_weight=1/physical[pair])); reuse_rows += 1
                    for col, feature in enumerate(FEATURES):
                        stats = balance(xx[:, col], yy[:, col], baseline[group[0]][:, col])
                        extra = {}
                        for label, base in [('metadata_matched', x), ('target_eligible_all', eb)]:
                            shift = balance(xx[:, col], yy[:, col], base[:, col])
                            extra.update({label + '_' + new: shift[old] for old, new in SHIFT_FIELDS.items()})
                        value = dict(zip(KEY,key), feature=feature, retained_feature_excluded_pairs=n-stats['pairs'], **stats, **extra)
                        if bw is None: bw = csv.DictWriter(bf,list(value),delimiter='\t',lineterminator='\n'); bw.writeheader()
                        bw.writerow(value); balance_rows += 1
            if strata % 12 == 0: print('Complete screened balance strata', strata, '/', len(plan['guides'])*len(plan['policies'])*plan['expected']['scenarios'], flush=True)
    assert coverage_rows == len(source['attrition']) == strata*3*len(plan['screens']) and balance_rows == coverage_rows*len(FEATURES)
    verify(bindings)
    result = dict(status='complete_full_screened_balance_pending_independent_readback',plan_sha256=sha(plan_path),matching_receipt_sha256=sha(source['prior_receipt']),target_nodes=len(nodes['target']),background_nodes=len(nodes['background']),selected_records=selected,scenarios=plan['expected']['scenarios'],strata=strata,coverage_rows=coverage_rows,balance_rows=balance_rows,reuse_rows=reuse_rows,retained_selection_screen_cells=retained_cells,features=FEATURES,screens=plan['screens'],masks=MASKS,guides=plan['guides'],policies=plan['policies'],source_hashes=bindings,artifacts={name:sha(out/name) for name in ['coverage.tsv','balance.tsv','retained_control_reuse.tsv.gz']},scientific_eligibility=False,scope=plan['scope'])
    with (out/'receipt.json').open('x') as handle: handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']},indent=2)); return result


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);run(parser.parse_args().plan)
