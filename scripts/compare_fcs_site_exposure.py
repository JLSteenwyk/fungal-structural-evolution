#!/usr/bin/env python3
"""Compare extant exposure quantiles before/after FCS omissions on identical sites."""
import argparse
import csv
import gzip
import json
import math
from collections import defaultdict
from pathlib import Path
import numpy as np
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha, read_table
from prepare_paired_phylogenetic_inputs import write_table


def quantile(values, q):
    values = sorted(values)
    if not values:
        return ''
    x = (len(values) - 1) * q
    lo, hi = math.floor(x), math.ceil(x)
    answer = values[lo] + (x - lo) * (values[hi] - values[lo])
    if not math.isclose(answer, float(np.quantile(values, q)), abs_tol=1e-12, rel_tol=1e-12):
        raise ValueError('Independent quantile implementations differ')
    return answer


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['baseline', 'sensitivity', 'inputs', 'baseline-inputs', 'readback', 'output']:
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    receipts = {k: checked_receipt(getattr(a, k)) for k in ['baseline', 'sensitivity', 'inputs', 'baseline_inputs']}
    rb = json.loads(a.readback.read_text())
    if rb['status'] != 'passed_full_accessibility_normalization_readback' or rb['normalized_receipt_sha256'] != sha(a.sensitivity / 'receipt.json'):
        raise ValueError('Sensitivity normalization readback differs')
    if receipts['inputs']['source_input_receipt_sha256'] != sha(a.baseline_inputs / 'receipt.json'):
        raise ValueError('Baseline input differs')
    # Bind each normalized table to the input receipt through its projection.
    for label, input_path in [('baseline', a.baseline_inputs), ('sensitivity', a.inputs)]:
        target = receipts[label]['source_receipt_sha256']
        root = getattr(a, label).parent
        matches = []
        for folder in root.glob('paired-accessibility-*'):
            rp = folder / 'receipt.json'
            if rp.exists() and sha(rp) == target:
                pr = checked_receipt(folder)
                if pr['source_receipts']['inputs'] != sha(input_path / 'receipt.json'):
                    raise ValueError('Normalized/input binding differs')
                matches.append(folder)
        if len(matches) != 1:
            raise ValueError('Expected exactly one matching projection receipt')
    ready = {r['marker'] for r in read_table(a.inputs / 'marker_summary.tsv') if r['status'] == 'ready_for_inference'}
    tables = {}
    rows_scanned = {}
    for label in ['baseline', 'sensitivity']:
        grouped = defaultdict(dict)
        n = 0
        with gzip.open(getattr(a, label) / 'normalized_paired_sites.tsv.gz', 'rt') as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                n += 1
                if row['marker'] not in ready:
                    continue
                key = row['marker'], int(row['paired_column_1based'])
                taxon = row['taxon_id']
                if taxon in grouped[key]:
                    raise ValueError('Duplicate site/taxon')
                grouped[key][taxon] = row
        if n != receipts[label]['rows']:
            raise ValueError('Normalized row count differs')
        tables[label] = grouped
        rows_scanned[label] = n
    if set(tables['baseline']) != set(tables['sensitivity']):
        raise ValueError('Changed site universe')
    output = []
    retained = omitted = 0
    summaries = []
    fields = ['rsa_tien2013_theoretical', 'rsa_miller1987', 'ca_plddt']
    for marker in sorted(ready):
        alignments = {}
        for label, folder in [('baseline', a.baseline_inputs), ('sensitivity', a.inputs)]:
            alignments[label] = {r.id: str(r.seq) for r in SeqIO.parse(folder / marker / 'aa.faa', 'fasta')}
        cols = read_table(a.inputs / marker / 'columns.tsv')
        start = len(output)
        for index, col in enumerate(cols, 1):
            key = marker, index
            before, after = tables['baseline'].pop(key), tables['sensitivity'].pop(key)
            for label, local in [('baseline', before), ('sensitivity', after)]:
                expected = {t for t, s in alignments[label].items() if s[index - 1] != '?'}
                if set(local) != expected:
                    raise ValueError('Observed alignment grid differs')
                if any(r['amino_acid'] != alignments[label][t][index - 1] or r['matrix_column_1based'] != col['matrix_column_1based'] for t, r in local.items()):
                    raise ValueError('Source residue/column differs')
            if not set(after) <= set(before) or any(row != before[t] for t, row in after.items()):
                raise ValueError('Retained source values changed')
            retained += len(after)
            omitted += len(before) - len(after)
            row = {'marker': marker, 'paired_column_1based': index,
                   'matrix_column_1based': col['matrix_column_1based'],
                   'baseline_observed_taxa': len(before), 'sensitivity_observed_taxa': len(after),
                   'omitted_observed_taxa': len(before) - len(after)}
            for field in fields:
                for stat, q in [('q25', .25), ('median', .5), ('q75', .75)]:
                    values = {}
                    for label, local in [('baseline', before), ('sensitivity', after)]:
                        v = [float(r[field]) for r in local.values() if r[field] != '']
                        if any(not math.isfinite(x) for x in v):
                            raise ValueError('Nonfinite observation')
                        values[label] = quantile(v, q)
                        row[label + '_' + field + '_' + stat] = values[label]
                    row['delta_' + field + '_' + stat] = values['sensitivity'] - values['baseline'] if all(v != '' for v in values.values()) else ''
            output.append(row)
        local = output[start:]
        summary = {'marker': marker, 'sites': len(local),
                   'sites_with_omitted_observation': sum(r['omitted_observed_taxa'] > 0 for r in local),
                   'omitted_residues': sum(r['omitted_observed_taxa'] for r in local)}
        for field in fields:
            for stat in ['q25', 'median', 'q75']:
                ds = [abs(r['delta_' + field + '_' + stat]) for r in local if r['delta_' + field + '_' + stat] != '']
                summary['max_abs_delta_' + field + '_' + stat] = max(ds) if ds else ''
        summaries.append(summary)
    if tables['baseline'] or tables['sensitivity'] or retained != receipts['sensitivity']['rows']:
        raise ValueError('Incomplete site comparison')
    a.output.mkdir(parents=True)
    write_table(a.output / 'site_exposure_changes.tsv', output)
    write_table(a.output / 'marker_summary.tsv', summaries)
    result = {'status': 'complete_matched_site_exposure_omission_comparison',
              'markers': len(ready), 'sites': len(output), 'retained_observations': retained,
              'omitted_observations': omitted,
              'sites_with_omitted_observation': sum(r['omitted_observed_taxa'] > 0 for r in output),
              'rows_scanned': rows_scanned,
              'max_abs_median_rsa_change_tien': max(abs(r['delta_rsa_tien2013_theoretical_median']) for r in output),
              'max_abs_median_rsa_change_miller': max(abs(r['delta_rsa_miller1987_median']) for r in output),
              'source_receipts': {k: sha(getattr(a, k) / 'receipt.json') for k in receipts},
              'normalization_readback_sha256': sha(a.readback), 'script_sha256': sha(Path(__file__)),
              'artifacts': {f.name: sha(f) for f in a.output.iterdir()},
              'interpretation': 'All matched extant-site exposure and confidence quantiles compared with exact retained-row and full alignment-grid checks. Quantiles independently checked by NumPy. Descriptive composition sensitivity, not phylogenetically adjusted exposure, ancestral reconstruction, rate inference, or coupling test.'}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
