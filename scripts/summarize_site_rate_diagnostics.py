#!/usr/bin/env python3
"""Summarize audited rate fits without treating warnings as model adequacy tests."""
import argparse
import json
from pathlib import Path
from statistics import median
from collections import Counter, defaultdict
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import read_table, sha
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--audits', nargs='+', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    summary, flags, warnings, pins = [], [], [], {}
    for folder in a.audits:
        r = checked_receipt(folder)
        if r['status'] != 'passed_full_site_rate_output_audit':
            raise ValueError('Completed rate audit required')
        pins[str(folder)] = sha(folder / 'receipt.json')
        fits = read_table(folder / 'fit_summary.tsv')
        if len(fits) != r['fits']:
            raise ValueError('Fit count mismatch')
        groups = defaultdict(list)
        for row in fits:
            groups[row['fit']].append(row)
            delta = float(row['log_likelihood_change'])
            if abs(delta) > .1:
                flags.append(dict(cohort=folder.name, marker=row['marker'], fit=row['fit'],
                    log_likelihood_change=delta, flag='absolute_export_refit_change_above_0.1'))
        for label, rows in sorted(groups.items()):
            deltas = [float(x['log_likelihood_change']) for x in rows]
            summary.append(dict(cohort=folder.name, fit=label, markers=len(rows),
                site_rate_rows=sum(int(x['sites']) for x in rows),
                fits_with_warnings=sum(int(x['warnings']) > 0 for x in rows),
                minimum_log_likelihood_change=min(deltas), median_log_likelihood_change=median(deltas),
                maximum_log_likelihood_change=max(deltas),
                export_refits_worse_by_more_than_0_1=sum(x < -.1 for x in deltas),
                export_refits_better_by_more_than_0_1=sum(x > .1 for x in deltas)))
        warning_rows = read_table(folder / 'warnings.tsv') if (folder / 'warnings.tsv').exists() else []
        for message, count in Counter(x['warning'] for x in warning_rows).most_common():
            warnings.append(dict(cohort=folder.name, warning=message, fit_occurrences=count))
    a.output.mkdir(parents=True)
    for name, rows in [('model_summary.tsv', summary), ('likelihood_change_flags.tsv', flags), ('warning_frequencies.tsv', warnings)]:
        if rows:
            write_table(a.output / name, rows)
    result = dict(status='complete_descriptive_site_rate_diagnostics', source_audits=pins,
        script_sha256=sha(Path(__file__)), fits=sum(x['markers'] for x in summary),
        flagged_likelihood_changes=len(flags),
        interpretation='Likelihood changes compare repeated exports of each fit within a cohort, not predictors or different models. Absolute 0.1 is a descriptive review threshold, not a significance cutoff. Warnings are retained, not automatic exclusions. Audit success does not establish convergence, model adequacy or calibrated rate uncertainty; pending optimization and rate-model sensitivity remain required.',
        artifacts={f.name: sha(f) for f in a.output.iterdir()})
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
