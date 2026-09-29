#!/usr/bin/env python3
"""Recheck the historically discrepant setting with explicit compensated sums."""
import csv
import json
import math
from pathlib import Path

import pandas as pd
from diagnose_whole_protein_record_readback import METRICS, KEYS, SETTING
from screen_duplication_alignment_reuse import sha


def average(values):
    finite = [float(x) for x in values if not pd.isna(x)]
    assert all(math.isfinite(x) for x in finite)
    return math.fsum(finite) / len(finite) if finite else None


def main():
    root = Path('results/structural_comparisons/matched-whole-protein-records-20260928-v1')
    receipt = json.loads((root / 'receipt.json').read_text())
    bindings = {str(root / 'receipt.json'): sha(root / 'receipt.json'),
                **receipt['source_hashes']}
    bindings.update({str(root / k): v for k, v in receipt['artifacts'].items()})
    for path, digest in bindings.items():
        assert sha(path) == digest, path
    setting = ('plddt70', 'same_mask', 'n50_c70', '0', '0')
    parts = json.loads((root / 'partition_manifest.json').read_text())
    part, = [p for p in parts if tuple(p[k] for k in SETTING) == setting]
    path = root / part['path']
    bindings[str(path)] = part['sha256']
    assert sha(path) == part['sha256']
    values = pd.read_parquet(path)
    selected = pd.read_csv('results/orthology/background-control-selection-20260927-v1/selections.tsv.gz',
                           sep='\t', usecols=['target_id', 'background_id', 'policy', 'scenario_id'])
    records = selected.merge(values, on=['target_id', 'background_id'], validate='many_to_one')
    groups = records.groupby(KEYS).indices
    rows = [r for r in csv.DictReader((root / 'record_summary.tsv').open(), delimiter='\t')
            if tuple(r[k] for k in SETTING) == setting]
    assert len(rows) == 432
    checked = 0
    max_error = 0.0
    original = None
    for row in rows:
        key = tuple(row[k] for k in KEYS)
        frame = records.iloc[groups.get(key, [])]
        assert len(frame) == int(row['retained_records'])
        for metric in METRICS:
            estimates = {'record_mean': average(frame[metric])}
            for label, suffix in [('family', 'family_equal_mean'), ('focal_taxon', 'taxon_equal_mean')]:
                estimates[suffix] = average([average(sub[metric]) for _, sub in frame.groupby(label)])
            for suffix, estimate in estimates.items():
                stored = row[metric + '_' + suffix]
                if estimate is None:
                    assert stored == '', (key, metric, suffix)
                else:
                    error = abs(float(stored) - estimate)
                    max_error = max(max_error, error)
                    assert math.isclose(float(stored), estimate, rel_tol=1e-9, abs_tol=1e-10), (key, metric, suffix, stored, estimate)
                checked += 1
            if key == ('mafft', 'envelope_bitscore', 'S41') and metric == 'rmsd_difference':
                original = dict(records=len(frame), stored=float(row[metric + '_record_mean']), fsum_mean=estimates['record_mean'])
    assert original is not None and checked == 432 * 36
    for path, digest in bindings.items():
        assert sha(path) == digest, path
    result = dict(status='passed_compensated_sum_check_of_all_failed_setting_summaries',
                  setting=setting, groups=len(rows), cells_checked=checked,
                  maximum_absolute_error=max_error, original_failure_group=original,
                  source_hashes=bindings, script_sha256=sha(__file__),
                  scope='All 432 strata and 12 metrics under three weighting schemes in the historically discrepant setting. Explicit math.fsum within each record/family/taxon group; no pandas or SQL mean. This does not identify the cause of the historical failure or validate other settings.')
    with Path('metadata/whole_protein_failed_setting_fsum_completed_20260928.json').open('x') as out:
        json.dump(result, out, indent=2)
        out.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}), flush=True)


if __name__ == '__main__':
    main()
