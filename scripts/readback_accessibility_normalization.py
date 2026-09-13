#!/usr/bin/env python3
"""Check every normalized row against projected ASA and rebuild scale summaries."""
import argparse
import csv
import gzip
import json
import math
from collections import Counter
from itertools import zip_longest
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha, read_table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['projection', 'normalized', 'snapshot', 'config', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    a = parser.parse_args()
    if a.output.exists():
        raise FileExistsError('Use an immutable readback output')
    source = checked_receipt(a.projection)
    result = checked_receipt(a.normalized)
    checked_receipt(a.snapshot)
    config = json.loads(a.config.read_text())
    if (result['status'] != 'complete_reference_normalization'
            or result['source_receipt_sha256'] != sha(a.projection / 'receipt.json')
            or result['snapshot_receipt_sha256'] != sha(a.snapshot / 'receipt.json')
            or source['source_receipts']['snapshot'] != result['snapshot_receipt_sha256']
            or result['normalization_config_sha256'] != sha(a.config)
            or sha(Path(config['source_path'])) != config['source_sha256']):
        raise ValueError('Source linkage differs')
    models = json.loads((a.snapshot / 'model_provenance.json').read_text())
    lengths = {m['model_id']: m['length'] for m in models}
    if len(lengths) != len(models):
        raise ValueError('Duplicate model identity')
    scales = config['scales']; thresholds = config['diagnostic_thresholds']
    if len(scales) != 2:
        raise ValueError('Expected two reference scales')
    counts = Counter(); discordant = Counter(); n = terminal = 0
    with gzip.open(a.projection / 'paired_site_accessibility.tsv.gz', 'rt') as src, gzip.open(a.normalized / 'normalized_paired_sites.tsv.gz', 'rt') as dst:
        left = csv.DictReader(src, delimiter='\t'); right = csv.DictReader(dst, delimiter='\t')
        if set(right.fieldnames) != set(left.fieldnames) | {'normalization_status'} | {'rsa_' + s for s in scales}:
            raise ValueError('Unexpected normalized columns')
        for before, after in zip_longest(left, right):
            if before is None or after is None or any(after[k] != v for k, v in before.items()):
                raise ValueError('Source row count, order or fields changed')
            pos = int(before['protein_residue_1based']); length = lengths[before['model_id']]
            area = float(before['sasa_angstrom_squared'])
            if not 1 <= pos <= length or not math.isfinite(area) or area < 0:
                raise ValueError('Invalid source position or area')
            n += 1
            if pos in (1, length):
                terminal += 1
                if after['normalization_status'] != 'terminal_not_normalized' or any(after['rsa_' + s] != '' for s in scales):
                    raise ValueError('Incorrect terminal treatment')
            else:
                if after['normalization_status'] != 'internal_residue_normalized':
                    raise ValueError('Incorrect internal-residue status')
                classifications = []
                for scale, maxima in scales.items():
                    denominator = maxima[before['amino_acid']]
                    value = float(after['rsa_' + scale])
                    if not math.isfinite(value) or value < 0 or not math.isclose(value * denominator, area, rel_tol=1e-12, abs_tol=1e-12):
                        raise ValueError('Normalized value does not recover source ASA')
                    expected = area / denominator
                    counts[scale, 'normalized_rows'] += 1
                    counts[scale, 'above_one'] += expected > 1
                    labels = [expected < t for t in thresholds]
                    classifications.append(labels)
                    for t, label in zip(thresholds, labels):
                        counts[scale, 'below_' + str(t)] += label
                for t, x, y in zip(thresholds, *classifications):
                    discordant[t] += x != y
            if n % 500000 == 0:
                print('Read back', n, 'rows', flush=True)
    if n != source['observed_sites_linked'] or n != result['rows'] or terminal != result['terminal_rows_not_normalized']:
        raise ValueError('Aggregate row counts differ')
    summaries = read_table(a.normalized / 'scale_summary.tsv')
    if len(summaries) != len(scales) or {r['scale'] for r in summaries} != set(scales):
        raise ValueError('Scale summary grid differs')
    for row in summaries:
        for field in ['normalized_rows', 'above_one'] + ['below_' + str(t) for t in thresholds]:
            if int(row[field]) != counts[row['scale'], field]:
                raise ValueError('Scale summary count differs')
    sensitivity = read_table(a.normalized / 'threshold_sensitivity.tsv')
    if len(sensitivity) != len(thresholds) or {float(r['threshold']) for r in sensitivity} != set(thresholds):
        raise ValueError('Threshold grid differs')
    for row in sensitivity:
        t = float(row['threshold'])
        if int(row['eligible_rows']) != n-terminal or int(row['scale_disagreements']) != discordant[t] or not math.isclose(float(row['fraction_disagreement']), discordant[t] / (n-terminal), rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError('Threshold summary differs')
    receipt = {'status': 'passed_full_accessibility_normalization_readback', 'rows': n,
               'terminal_rows': terminal, 'projection_receipt_sha256': sha(a.projection / 'receipt.json'),
               'normalized_receipt_sha256': sha(a.normalized / 'receipt.json'),
               'script_sha256': sha(Path(__file__)),
               'scope': 'Every source row preserved, normalized value recovers projected ASA, terminal treatment and all scale/threshold totals checked. Does not independently verify projection against raw residue tables, recompute solvent area or validate biological exposure.'}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt), flush=True)


if __name__ == '__main__':
    main()
