#!/usr/bin/env python3
"""Verify every threshold count, empirical quantile and descriptive rank correlation."""
import argparse
from collections import Counter
import csv
from datetime import datetime,timezone
import json
import math
from pathlib import Path

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind,verify

METRICS=['ca_superposition_rmsd_angstrom','all_distance_rms_change_angstrom',
         'local_distance_rms_change_angstrom','sequence_local_distance_rms_change_angstrom']


def ranks(values):
    counts=Counter(values);result={};position=1
    for value in sorted(counts):
        count=counts[value];result[value]=position+(count-1)/2;position+=count
    return [result[value] for value in values]


def rho(x,y):
    if len(x)<2 or min(x)==max(x) or min(y)==max(y):return ''
    xx,yy=ranks(x),ranks(y);mx=math.fsum(xx)/len(xx);my=math.fsum(yy)/len(yy)
    xy=math.fsum((a-mx)*(b-my) for a,b in zip(xx,yy))
    sx=math.fsum((a-mx)**2 for a in xx);sy=math.fsum((b-my)**2 for b in yy)
    return xy/math.sqrt(sx*sy)


def quantile(values,p):
    if not values:return ''
    order=sorted(values);rank=(len(order)-1)*p;lo=math.floor(rank);hi=math.ceil(rank)
    return order[lo]+(rank-lo)*(order[hi]-order[lo])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['summary','transport','receipt']:p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();assert not a.receipt.exists()
    s=json.loads(a.summary.read_text());t=json.loads(a.transport.read_text())
    assert s['status']=='complete_descriptive_full_overlap_coordinate_summary_and_figure'
    assert t['validation_sha256']==sha(a.summary) and t['original_tool_terminal_exit_code']==0
    pins=dict(s['source_hashes'])
    for q in [a.summary,a.transport,Path(__file__)]:bind(pins,q)
    verify(pins)
    rows=list(csv.DictReader(Path('results/phylogeny/overlap-coordinate-benchmark-20261004-v1/all_model_pair_coordinates.tsv').open(),delimiter='\t'))
    summaries=list(csv.DictReader(Path('results/phylogeny/overlap-coordinate-summary-20261004-v1/threshold_geometry_summary.tsv').open(),delimiter='\t'))
    corrs=list(csv.DictReader(Path('results/phylogeny/overlap-coordinate-summary-20261004-v1/descriptive_geometry_correlations.tsv').open(),delimiter='\t'))
    assert len(rows)==7716 and len(summaries)==12 and len(corrs)==96
    summary_keys=set();corr_keys=set();maximum_difference=0.
    def compare(saved,expected):
        nonlocal maximum_difference
        if expected=='':assert saved==''
        else:
            value=float(saved);assert math.isfinite(value)
            maximum_difference=max(maximum_difference,abs(value-expected))
            assert math.isclose(value,expected,rel_tol=1e-12,abs_tol=1e-12),(saved,expected)
    for summary in summaries:
        cutoff,pae=summary['plddt_cutoff'],summary['pae_cutoff'];key=(cutoff,pae)
        assert key not in summary_keys;summary_keys.add(key)
        group=[r for r in rows if (r['plddt_cutoff'],r['pae_cutoff'])==key]
        good=[r for r in group if r['status']=='coordinates_compared'];assert len(group)==643
        counts=dict(declared_model_pairs=643,compared_model_pairs=len(good),
            insufficient_common_coverage=sum(r['status']=='insufficient_common_coverage' for r in group),
            coordinate_validation_rejected=sum(r['status']=='coordinate_validation_rejected' for r in group),
            retained_residues=sum(int(r['retained_residues']) for r in good),
            state_mismatches=sum(int(r['state_mismatches']) for r in good))
        assert all(int(summary[k])==v for k,v in counts.items())
        for metric in METRICS:
            v=[float(r[metric]) for r in good]
            for label,q in [('p025',.025),('median',.5),('p975',.975)]:compare(summary[metric+'_'+label],quantile(v,q))
        for corr in [c for c in corrs if (c['plddt_cutoff'],c['pae_cutoff'])==key]:
            tag=key+(corr['x'],corr['y']);assert tag not in corr_keys;corr_keys.add(tag)
            assert int(corr['model_pairs'])==len(good)
            assert corr['x'] in ['state_mismatch_fraction','partner_change_fraction'] and corr['y'] in METRICS
            field='state_mismatches' if corr['x']=='state_mismatch_fraction' else 'partner_changes'
            xx=[int(r[field])/int(r['retained_residues']) for r in good]
            yy=[float(r[corr['y']]) for r in good]
            compare(corr['descriptive_spearman_rho'],rho(xx,yy))
            assert corr['scope']=='Dependent selected predictor pairs; no p-value or evolutionary attribution'
    assert summary_keys=={(c,p) for c in ['0','70','90'] for p in ['unfiltered','5','10','15']}
    assert len(corr_keys)==96;verify(pins)
    result=dict(status='passed_full_independent_coordinate_summary_readback',checked_utc=datetime.now(timezone.utc).isoformat(),
        summary_receipt=str(a.summary),summary_receipt_sha256=sha(a.summary),source_rows=7716,threshold_summary_rows=12,
        descriptive_correlation_rows=96,maximum_summary_absolute_difference=maximum_difference,
        summary_relative_tolerance=1e-12,summary_absolute_tolerance=1e-12,source_hashes=pins,
        scientific_eligibility=False,all_eight_aims_incomplete=True,
        scope='Every coverage count/empirical quantile/rank correlation independently verified with sorted linear ranks, '
              'explicit tied midranks and scalar math.fsum moments. Figure/source checksums verified; visual inspection separate. '
              'Summary roundoff tolerances do not change weighted numerical guards or qualify biological inference.')
    with a.receipt.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
