#!/usr/bin/env python3
"""Verify every plotted value against the audited raw comparison tables."""
import csv
import hashlib
import json
import math
from pathlib import Path


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def table(p):
    with Path(p).open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    result = {r['case_id']: r for r in rows}
    assert len(result) == len(rows)
    return result


def main():
    prefix = Path('docs/figures/codon_alignment_sensitivity_20260927')
    receipt = json.loads(prefix.with_suffix('.receipt.json').read_text())
    for name, digest in receipt['artifacts'].items(): assert sha(prefix.parent/name) == digest
    actual = table(prefix.with_suffix('.tsv'))
    base = table('results/cds/codon-alignment-divergence-comparison-20260927-v1/cases.tsv')
    base = {k:v for k,v in base.items() if v['comparison_status']=='matched'}
    trees = table('results/cds/codon-alignment-tree-comparison-20260927-v1/cases.tsv')
    assert set(actual) == set(base) and len(actual) == 1625
    counts = {}
    for case, a in actual.items():
        b, t = base[case], trees[case]
        for key in ['historical_review_flags','selection_eligibility']: assert a[key] == b[key]
        assert (a['topology_changed']=='True') == (int(t['rf_distance'])>0)
        assert math.isclose(float(a['median_pair_retention']),float(t['alignment_median_pair_retention']),rel_tol=0,abs_tol=1e-15)
        for metric in ['global_omega','tree_ds_equal_alternative','tree_dn_equal_alternative']:
            x,y = float(b['original_'+metric]),float(b['local_'+metric])
            v = a[metric+'_log2_local_over_original']
            if x>0 and y>0:
                assert math.isclose(float(v),math.log2(y/x),rel_tol=1e-10,abs_tol=1e-10)
                counts[metric] = counts.get(metric,0)+1
            else:
                assert metric=='global_omega' and x==y==0 and v=='' and a['global_omega_ratio_disposition']=='zero_both'
        assert a['global_omega_ratio_disposition'] == ('positive_both' if float(b['original_global_omega'])>0 and float(b['local_global_omega'])>0 else 'zero_both')
    assert counts == dict(global_omega=1624,tree_ds_equal_alternative=1625,tree_dn_equal_alternative=1625)
    proof = dict(status='passed_full_codon_alignment_figure_data_readback',cases=1625,ecdf_defined_cases=counts,topology_changed=sum(a['topology_changed']=='True' for a in actual.values()),figure_receipt_sha256=sha(prefix.with_suffix('.receipt.json')),script_sha256=sha(__file__),scope='Every exported ratio recomputed from audited comparison values; retention, flags, exact case grid and topology labels checked against source tables. Retention tolerance 1e-15 accommodates dataframe decimal serialization. Does not test significance or causal inference; rendered figure requires separate visual inspection.')
    Path('metadata/codon_alignment_figure_readback_20260927.json').write_text(json.dumps(proof,indent=2)+'\n')
    print(json.dumps(proof,indent=2))


if __name__=='__main__':main()
