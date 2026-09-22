#!/usr/bin/env python3
"""Summarize audited annotation variation with exact-membership guide controls."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import sqlite3
import time


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8388608), b''):
            h.update(b)
    return h.hexdigest()


def metrics(row):
    counts = Counter(families=1)
    for key in ['proteins', 'no_ga_hit_proteins', 'policy_disagreement_proteins',
                'rank_tie_proteins', 'unresolved_overlap_proteins', 'candidate_nested_proteins',
                'alignment_overlap_proteins', 'partial_hmm_proteins']:
        counts[key] = int(row[key])
    for prefix in ['observed_', 'conservative_']:
        counts[prefix+'proteins'] = int(row[prefix+'proteins'])
        multi = int(row[prefix+'taxa']) >= 2
        counts[prefix+'families_with_at_least_two_annotated_taxa'] = int(multi)
        for field in ['ordered_signatures', 'multisets', 'model_type_sets']:
            counts[prefix+'families_with_multiple_'+field] = int(int(row[prefix+field]) > 1)
        for field in ['multisets_with_order_variation', 'model_type_sets_with_multiplicity_variation']:
            varied = int(row[prefix+field]) > 0
            counts[prefix+'families_with_'+field] = int(varied)
            counts[prefix+'multi_taxon_families_with_'+field] = int(multi and varied)
    return counts


def compare(folder, db, summaries):
    guides = [r['guide'] for r in summaries]
    if len(guides) != 2 or len(set(guides)) != 2:
        raise ValueError('Expected two distinct guide partitions')
    membership = {}
    for guide in guides:
        membership[guide] = dict(db.execute('SELECT family,membership_sha256 FROM families WHERE guide=?', (guide,)))
        if len(membership[guide]) != next(r['families'] for r in summaries if r['guide'] == guide):
            raise ValueError('Wrong family universe')
        if len(set(membership[guide].values())) != len(membership[guide]):
            raise ValueError('Repeated membership within partition')
    shared = set(membership[guides[0]].values()) & set(membership[guides[1]].values())
    previous = {}
    totals = defaultdict(Counter)
    policies = {'alignment_evalue','alignment_bitscore','envelope_evalue','envelope_bitscore'}
    checked = 0
    for guide in guides:
        seen = set()
        with gzip.open(folder/(guide+'_family_architectures.tsv.gz'), 'rt') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                key = (row['family'], row['policy'])
                if (row['guide'] != guide or key in seen or row['family'] not in membership[guide]
                        or row['policy'] not in policies):
                    raise ValueError('Unexpected or duplicate row')
                seen.add(key)
                identity = membership[guide][row['family']]
                group = 'exact_shared_membership' if identity in shared else 'guide_specific_membership'
                summary = metrics(row)
                for scope in ['all_families', group]:
                    totals[guide, row['policy'], scope].update(summary)
                if identity in shared:
                    values = {k:v for k,v in row.items() if k not in ['guide','family']}
                    digest = hashlib.sha256(json.dumps(values, sort_keys=True).encode()).hexdigest()
                    shared_key = (identity, row['policy'])
                    if guide == guides[0]:
                        previous[shared_key] = digest
                    elif previous.pop(shared_key, None) != digest:
                        raise ValueError('Different annotations for identical family membership')
                    else:
                        checked += 1
        if len(seen) != len(membership[guide])*4:
            raise ValueError('Missing family/policy rows')
    if previous or checked != len(shared)*4:
        raise ValueError('Incomplete matched comparison')
    rows = [dict(guide=g, policy=p, membership_scope=s, **dict(c)) for (g,p,s),c in sorted(totals.items())]
    return rows, len(shared), checked


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', required=True, type=Path)
    a = ap.parse_args(); plan = json.loads(a.plan.read_text())
    pins = {str(a.plan): sha(a.plan), **plan['pins']}
    def verify():
        for path, digest in pins.items():
            if sha(path) != digest:
                raise ValueError('Changed input: '+path)
    verify()
    out = Path(plan['output']); out.mkdir(exist_ok=False)
    start = time.time()
    def state(stage):
        p = out/'state.tmp';p.write_text(json.dumps(dict(stage=stage, elapsed_seconds=time.time()-start))+'\n');p.replace(out/'state.json')
    state('waiting_for_exact_summary_auditor')
    proc = Path('/proc')/str(plan['auditor_pid'])/'stat'
    while proc.exists():
        try: fields = proc.read_text().rsplit(')',1)[1].split()
        except FileNotFoundError: break
        if fields[19] != str(plan['auditor_start_ticks']) or fields[0] == 'Z': break
        time.sleep(20)
    audit_plan = json.loads(Path(plan['auditor_plan']).read_text())
    source = json.loads(Path(audit_plan['producer_plan']).read_text())
    folder = Path(source['output']); receipt_path = folder/'receipt.json'
    audit_path = Path(audit_plan['output'])/'receipt.json'
    audit = json.loads(audit_path.read_text()); receipt = json.loads(receipt_path.read_text())
    if (audit['status'] != 'passed_full_family_architecture_variation_readback'
            or audit['producer_receipt_sha256'] != sha(receipt_path)
            or audit['plan_sha256'] != sha(plan['auditor_plan'])):
        raise ValueError('Unverified family summaries')
    pins.update(receipt['input_hashes'])
    for p in [receipt_path,audit_path]: pins[str(p)] = sha(p)
    pins.update({str(folder/n): d for n,d in receipt['artifacts'].items()})
    verify()
    if shutil.disk_usage(out).free < plan['resources']['minimum_free_disk_gib']*2**30:
        raise ValueError('Insufficient free disk')
    state('comparing_full_partitions')
    db = sqlite3.connect('file:'+str(Path(source['bridge_database']).resolve())+'?mode=ro', uri=True)
    rows, shared, checked = compare(folder, db, receipt['guides']);db.close()
    path = out/'guide_policy_architecture_summary.tsv'
    with path.open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n');w.writeheader();w.writerows(rows)
    verify()
    result = dict(status='complete_audited_family_architecture_guide_summary',
                  exact_shared_memberships=shared, matched_family_policy_rows_checked=checked,
                  summary_rows=len(rows), guides=receipt['guides'], plan_sha256=sha(a.plan),
                  input_hashes=pins, artifacts={path.name:sha(path)}, elapsed_seconds=time.time()-start,
                  scope='Descriptive annotation variation under both full guide partitions, with exact-membership controls. Shared-family metrics must agree exactly. Guide-specific membership is not a domain evolutionary event. Protein counts repeat across policy and scope rows; multi-taxon variation does not establish independent transitions, duplication, branch direction or statistical significance.')
    (out/'receipt.json').write_text(json.dumps(result, indent=2)+'\n');state(result['status'])


if __name__ == '__main__':
    main()
