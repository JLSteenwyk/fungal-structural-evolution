#!/usr/bin/env python3
"""Describe signed-distance sensitivity across all settings and tied references."""
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from screen_duplication_domain_alignment_coverage import sha

TOLERANCE = 1e-8  # Numerical zero only; not a biological effect-size threshold.


def direction(lo, hi):
    if lo > TOLERANCE:
        return 'a_farther'
    if hi < -TOLERANCE:
        return 'b_farther'
    if lo >= -TOLERANCE and hi <= TOLERANCE:
        return 'numerically_zero'
    return 'setting_or_reference_sensitive'


def main():
    source = Path('results/structural_comparisons/whole-protein-common-fit-event-links-20260927-v1')
    sr = source/'receipt.json'; r = json.loads(sr.read_text())
    assert r['status'] == 'complete_whole_protein_common_fit_event_eligibility'
    for path, digest in r['source_hashes'].items():
        assert sha(path) == digest
    for name, digest in r['artifacts'].items():
        assert sha(source/name) == digest
    fitroot = Path('results/structural_comparisons/whole-protein-common-core-fits-20260927-v1')
    fr = json.loads((fitroot/'receipt.json').read_text()); table = fitroot/'common_residue_fits.tsv'
    assert sha(table) == fr['artifacts'][table.name]
    cols = ['triad_id', 'rmsd_ar_minus_br', 'sequence_identity_ar', 'sequence_identity_br']
    frame = pd.read_csv(table, sep='\t', usecols=cols)
    frame['sequence_divergence_a_minus_b'] = frame.sequence_identity_br-frame.sequence_identity_ar
    groups = frame.groupby('triad_id', sort=True)
    ranges = groups.agg(structural_min=('rmsd_ar_minus_br','min'), structural_max=('rmsd_ar_minus_br','max'),
                        sequence_min=('sequence_divergence_a_minus_b','min'), sequence_max=('sequence_divergence_a_minus_b','max'),
                        numeric_configurations=('rmsd_ar_minus_br','count'))
    assert len(ranges) == 17619 and (groups.size() == 32).all()
    # Independently reconstruct every finite range from the full serialized table.
    scalar = defaultdict(lambda: [[], []])
    with table.open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            values = scalar[row['triad_id']]
            if row['rmsd_ar_minus_br']:
                values[0].append(float(row['rmsd_ar_minus_br']))
                values[1].append(float(row['sequence_identity_br'])-float(row['sequence_identity_ar']))
    for key, (structural, sequence) in scalar.items():
        got = ranges.loc[key]
        assert int(got.numeric_configurations) == len(structural)
        for prefix, values in [('structural',structural),('sequence',sequence)]:
            for label, function in [('min',min),('max',max)]:
                actual = got[prefix+'_'+label]
                assert math.isclose(actual, function(values), rel_tol=1e-12, abs_tol=1e-12) if values else math.isnan(actual)
    events = pd.read_csv(source/'event_reference_eligibility.tsv', sep='\t')
    links = events.merge(ranges.reset_index(), on='triad_id', how='left', validate='many_to_one')
    assert len(links) == 36944
    screens = ['n30_c50','n30_c70','n30_c90','n50_c50','n50_c70','n50_c90']
    keys = ['guide','family','gene_node','gene_a','gene_b']
    output = []
    for key, group in links.groupby(keys, sort=True):
        for screen in screens:
            eligible = group[(group.distinct_models == 3) & (group[screen+'_all_configurations_pass'] == 1)]
            row = dict(zip(keys,key), screen=screen, reference_links=len(group), eligible_reference_links=len(eligible),
                       all_references_eligible=int(len(eligible)==len(group)), structural_min='',structural_max='',
                       sequence_min='',sequence_max='', structural_direction='no_eligible_reference', sequence_direction='no_eligible_reference')
            if len(eligible):
                assert (eligible.numeric_configurations == 32).all()
                for prefix in ['structural','sequence']:
                    lo = float(eligible[prefix+'_min'].min()); hi = float(eligible[prefix+'_max'].max())
                    row[prefix+'_min']=lo; row[prefix+'_max']=hi; row[prefix+'_direction']=direction(lo,hi)
            output.append(row)
    result = pd.DataFrame(output)
    assert len(result) == (17448+17461)*6 and not result.duplicated(keys+['screen']).any()
    summary = result.groupby(['guide','screen','all_references_eligible','structural_direction'],sort=True).size().reset_index(name='events')
    assert summary.events.sum() == len(result)
    out = Path('results/structural_comparisons/whole-protein-reference-sensitivity-20260927-v1'); out.mkdir(exist_ok=False)
    artifacts = {}
    for name, data in [('triad_ranges.tsv',ranges.reset_index()),('event_reference_ranges.tsv',links),('event_sensitivity.tsv',result),('summary.tsv',summary)]:
        path = out/name; data.to_csv(path,sep='\t',index=False)
        back = pd.read_csv(path,sep='\t',keep_default_na=False)
        pd.testing.assert_frame_equal(back.fillna(''),data.fillna('').reset_index(drop=True),check_dtype=False,check_exact=False,rtol=1e-12,atol=1e-12)
        artifacts[name]=sha(path)
    for name,digest in r['artifacts'].items():
        assert sha(source/name)==digest
    assert sha(table)==fr['artifacts'][table.name]
    receipt = dict(status='complete_descriptive_whole_protein_reference_sensitivity',
                   source_hashes={str(sr):sha(sr),str(fitroot/'receipt.json'):sha(fitroot/'receipt.json')},
                   script_sha256=sha(__file__), numeric_zero_tolerance=TOLERANCE,
                   triads=len(ranges),event_reference_links=len(links),event_screen_rows=len(result),artifacts=artifacts,
                   scope='Ranges cover all 32 configurations within eligible references. All reference links retained; all_references_eligible separates partial coverage. Structural contrast is RMSD(A,R)-RMSD(B,R); sequence contrast is identity(B,R)-identity(A,R). Numerical sign robustness is descriptive, not biological significance, prediction uncertainty, an evolutionary rate or evidence of selection. Ranges are not confidence intervals.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__ == '__main__':
    main()
