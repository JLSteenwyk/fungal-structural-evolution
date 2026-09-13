#!/usr/bin/env python3
"""Join audited conditional site rates to exposure, composition and coverage."""
import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha, read_table
from prepare_paired_phylogenetic_inputs import write_table

LABELS = ['aa', '3di_af', '3di_af_empirical', '3di_llm']
AA = 'ACDEFGHIKLMNPQRSTVWY'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['rates', 'rate_readback', 'exposure', 'inputs', 'output']:
        p.add_argument('--' + name.replace('_', '-'), type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use a new immutable frame')
    receipts = {k: checked_receipt(getattr(a, k)) for k in ['rates', 'exposure', 'inputs']}
    rb = json.loads(a.rate_readback.read_text())
    if rb['status'] != 'passed_full_rate_comparison_readback' or rb['comparison_receipt_sha256'] != sha(a.rates / 'receipt.json'):
        raise ValueError('Matching full rate readback required')
    for name in ['rates', 'exposure']:
        if receipts[name]['source_receipts']['inputs'] != sha(a.inputs / 'receipt.json'):
            raise ValueError('Input snapshot differs')
    ready = {r['marker']: r for r in read_table(a.inputs / 'marker_summary.tsv') if r['status'] == 'ready_for_inference'}
    rates = defaultdict(dict)
    for row in read_table(a.rates / 'sites.tsv'):
        key = row['marker'], int(row['paired_column_1based']); label = row['fit']
        if label not in LABELS or label in rates[key] or key[0] not in ready:
            raise ValueError('Duplicate or unexpected rate identity')
        for field in ['gamma_rate', 'freerate_rate']:
            value = float(row[field])
            if not math.isfinite(value) or value < 0:
                raise ValueError('Invalid rate')
        rates[key][label] = row
    summaries = read_table(a.rates / 'fit_summary.tsv')
    if len(summaries) != 4 * len(ready) or {(r['marker'], r['fit']) for r in summaries} != {(m, f) for m in ready for f in LABELS}:
        raise ValueError('Fit summary grid differs')
    exposure_rows = read_table(a.exposure / 'site_parsimony_exposure.tsv')
    exposures = {(r['marker'], int(r['paired_column_1based'])): r for r in exposure_rows}
    if len(exposures) != len(exposure_rows) or set(exposures) != set(rates):
        raise ValueError('Exposure/rate site grid differs')
    output = []; seen = set(); observations = 0
    for marker, summary in ready.items():
        folder = a.inputs / marker
        aa = {r.id: str(r.seq) for r in SeqIO.parse(folder / 'aa.faa', 'fasta')}
        columns = read_table(folder / 'columns.tsv')
        if len(aa) != int(summary['eligible_taxa']) or len(columns) != int(summary['retained_columns']):
            raise ValueError('Input dimensions differ')
        for i, col in enumerate(columns):
            key = marker, i + 1
            base = exposures[key]
            if base['matrix_column_1based'] != col['matrix_column_1based'] or set(rates[key]) != set(LABELS):
                raise ValueError('Site correspondence differs')
            counts = Counter(s[i] for s in aa.values() if s[i] != '?'); n = sum(counts.values())
            if not set(counts) <= set(AA) or n != int(base['observed_taxa']) or len(counts) != int(base['aa_distinct_states']):
                raise ValueError('Observed amino-acid composition differs')
            row = dict(base, marker_taxa=len(aa), marker_columns=len(columns),
                       observed_fraction=n / len(aa),
                       aa_entropy_nats=-sum((c/n)*math.log(c/n) for c in counts.values()),
                       **{'aa_count_' + state: counts[state] for state in AA})
            for label in LABELS:
                for field in ['gamma_rate', 'freerate_rate']:
                    row[label + '_' + field] = rates[key][label][field]
            output.append(row); seen.add(key); observations += n
    if seen != set(rates) or len(seen) * 4 != receipts['rates']['site_comparisons'] or observations != receipts['exposure']['taxon_site_observations']:
        raise ValueError('Incomplete full site/rate/observation grid')
    a.output.mkdir(parents=True)
    write_table(a.output / 'site_rate_exposure.tsv', output)
    write_table(a.output / 'fit_diagnostics.tsv', summaries)
    result = {'status': 'complete_site_rate_exposure_analysis_frame', 'markers': len(ready),
              'sites': len(output), 'site_rate_values': len(output)*8, 'taxon_site_observations': observations,
              'source_receipts': {k: sha(getattr(a, k) / 'receipt.json') for k in receipts},
              'rate_readback_sha256': sha(a.rate_readback), 'script_sha256': sha(Path(__file__)),
              'artifacts': {p.name: sha(p) for p in a.output.iterdir()},
              'interpretation': 'All paired sites retained with eight conditional model-relative rate estimates, extant exposure quantiles under two scales, amino-acid composition, coverage and copy-review status. Rates derive from phylogenetic likelihood fits but exposure summaries are not phylogenetically adjusted or ancestral. Model labels are not predictors; rates across models are not a shared physical scale. This is an analysis frame, not a coupling/exposure significance test. Family effects, topology and rate uncertainty, nonlocal feature dependence, prediction confidence, source circularity and missingness require further controls.'}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
