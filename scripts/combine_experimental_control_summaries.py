#!/usr/bin/env python3
"""Combine disjoint, audited protein control cohorts while retaining source tiers."""
import argparse
import json
from pathlib import Path
from statistics import median
from audit_joint_path_uncertainty import checked, rows, sha, quantile
from prepare_paired_phylogenetic_inputs import write_table
from summarize_experimental_predictor_controls import METRICS, DELTAS


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--summaries', nargs='+', type=Path, required=True)
    ap.add_argument('--readbacks', nargs='+', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    if len(a.summaries) != len(a.readbacks):
        raise ValueError('One readback required per summary')
    if a.output.exists():
        raise FileExistsError(a.output)
    proteins, dispositions, sources, seen = [], [], [], set()
    for folder, audit in zip(a.summaries, a.readbacks):
        receipt = checked(folder)
        readback = json.loads(audit.read_text())
        if readback['status'] != 'passed_hierarchical_control_summary_readback' or readback['summary_receipt_sha256'] != sha(folder/'receipt.json'):
            raise ValueError('Matching successful summary audit required')
        p = list(rows(folder/'protein_summary.tsv'))
        d = list(rows(folder/'protein_disposition.tsv'))
        universe = {x['sequence_sha256'] for x in d}
        if seen & universe:
            raise ValueError('Overlapping sequences cannot be counted twice')
        expected = {(s,c) for s in universe for c in ['0','70','90']}
        identity = lambda x: (x['sequence_sha256'], x['joint_predicted_plddt_cutoff'])
        if len(d) != len(expected) or {identity(x) for x in d} != expected or len(universe) != receipt['control_sequences']:
            raise ValueError('Incomplete disposition universe')
        included = {identity(x) for x in d if x['status'] == 'included'}
        if len(p) != len(included) or {identity(x) for x in p} != included:
            raise ValueError('Protein/disposition mismatch')
        seen.update(universe)
        for values, destination in [(p,proteins),(d,dispositions)]:
            for row in values:
                destination.append(dict(row, source_tier=folder.name))
        sources.append({'summary':str(folder), 'summary_receipt_sha256':sha(folder/'receipt.json'), 'readback':str(audit), 'readback_sha256':sha(audit), 'control_sequences':len(universe)})
    summaries = []
    for tier in ['combined'] + [f.name for f in a.summaries]:
        tier_proteins = [x for x in proteins if tier == 'combined' or x['source_tier'] == tier]
        common = set.intersection(*[{x['sequence_sha256'] for x in tier_proteins if x['joint_predicted_plddt_cutoff'] == c} for c in ['0','70','90']])
        for cohort in ['all_available','same_proteins_all_thresholds']:
            for cutoff in ['0','70','90']:
                subset = [x for x in tier_proteins if x['joint_predicted_plddt_cutoff'] == cutoff and (cohort == 'all_available' or x['sequence_sha256'] in common)]
                for metric in METRICS + DELTAS:
                    vals = [float(x[metric]) for x in subset if x[metric] != '']
                    summaries.append({'source_tier':tier,'cohort':cohort,'joint_predicted_plddt_cutoff':cutoff,'metric':metric,'proteins':len(subset),'proteins_with_metric':len(vals),'median':median(vals) if vals else '', 'q25':quantile(vals,.25) if vals else '', 'q75':quantile(vals,.75) if vals else '', 'positive':sum(x>1e-8 for x in vals) if metric in DELTAS else '', 'negative':sum(x< -1e-8 for x in vals) if metric in DELTAS else '', 'approximately_zero':sum(abs(x)<=1e-8 for x in vals) if metric in DELTAS else ''})
    a.output.mkdir(parents=True)
    for name, values in [('protein_summary.tsv',proteins),('protein_disposition.tsv',dispositions),('cohort_summary.tsv',summaries)]:
        write_table(a.output/name, values)
    result = {'status':'complete_disjoint_control_summary_combination','control_sequences':len(seen),'sources':sources,'script_sha256':sha(Path(__file__)), 'interpretation':'Equal canonical-protein descriptive summaries of previously audited within-protein hierarchical medians. Source tiers remain explicit. Sequence disjointness is not statistical or evolutionary independence. Confidence filters change residues and eligible proteins; no causal length effect, significance or general accuracy claim. Additional pending controls and alternative chunk settings are not included.', 'artifacts':{p.name:sha(p) for p in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__ == '__main__':
    main()
