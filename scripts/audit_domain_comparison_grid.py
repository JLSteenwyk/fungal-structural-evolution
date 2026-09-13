#!/usr/bin/env python3
"""Account for every baseline pair/shared eligible domain and validate result fields."""
import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha, read_table


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['domains', 'comparisons', 'output']:
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use immutable output')
    r = checked_receipt(a.domains); checked_receipt(a.comparisons)
    if r['sources']['comparisons']['receipt_sha256'] != sha(a.comparisons / 'receipt.json'):
        raise ValueError('Baseline receipt mismatch')
    baseline = {}
    for row in read_table(a.comparisons / 'pairwise_metrics.tsv'):
        if row['plddt_cutoff'] != '70':
            continue
        key = row['marker'], row['taxon_a'], row['taxon_b']
        if key in baseline:
            raise ValueError('Duplicate baseline pair')
        baseline[key] = row
    if len(baseline) != r['full_marker_pairs_screened']:
        raise ValueError('Baseline count mismatch')
    hits = read_table(a.domains / 'domain_eligibility.tsv')
    eligible = defaultdict(dict)
    for hit in hits:
        if hit['status'] == 'eligible':
            key = hit['marker'], hit['taxon_id']; family = hit['pfam_accession']
            if family in eligible[key]:
                raise ValueError('Multiple eligible copies')
            eligible[key][family] = hit
    if len(hits) != r['mapped_annotation_hits'] or sum(len(v) for v in eligible.values()) != r['eligible_annotation_hits']:
        raise ValueError('Annotation count mismatch')
    expected = set()
    for marker, ta, tb in baseline:
        families = set(eligible[marker, ta]) & set(eligible[marker, tb])
        for family in families or {''}:
            expected.add((marker, ta, tb, family))
    seen = set(); exclusions = Counter(); accepted_count = 0; accepted_pairs = set(); families_seen = set()
    for accepted, filename in [(True, 'domain_comparisons.tsv'), (False, 'excluded_pairs.tsv')]:
        for row in read_table(a.domains / filename):
            key = row['marker'], row['taxon_a'], row['taxon_b'], row['pfam_accession']
            if key not in expected or key in seen:
                raise ValueError('Unexpected or duplicate domain disposition')
            seen.add(key)
            if not accepted:
                expected_reason = 'no_shared_eligible_single_instance_domain' if not key[-1] else 'fewer_than_30_qualified_or_half_domain_shared_positions'
                if row['reason'] != expected_reason:
                    raise ValueError('Incorrect exclusion class')
                exclusions[row['reason']] += 1
                continue
            if not key[-1]:
                raise ValueError('Accepted missing domain')
            accepted_count += 1; accepted_pairs.add(key[:3]); families_seen.add(key[-1])
            for suffix, taxon in [('a', key[1]), ('b', key[2])]:
                hit = eligible[key[0], taxon][key[3]]
                if row['hit_' + suffix] != hit['hit_id'] or row['domain_' + suffix + '_start'] != hit['alignment_start'] or row['domain_' + suffix + '_end'] != hit['alignment_end']:
                    raise ValueError('Domain correspondence mismatch')
            n = int(row['domain_compared_residues']); shared = int(row['domain_shared_positions'])
            if not 30 <= n <= shared or 2*n < shared:
                raise ValueError('Domain coverage mismatch')
            b = baseline[key[:3]]
            if int(row['whole_marker_compared_residues']) != int(b['compared_residues']) or float(row['whole_marker_rmsd_angstrom']) != float(b['ca_superposition_rmsd_angstrom']):
                raise ValueError('Whole-marker baseline mismatch')
            for field, value in row.items():
                if value != '' and ('angstrom' in field or 'fraction' in field or field == 'domain_uncorrected_sequence_difference'):
                    number = float(value)
                    if not math.isfinite(number) or number < -1e-7:
                        raise ValueError('Invalid metric')
            own = float(row['domain_own_fit_rmsd_angstrom']); whole = float(row['domain_sites_under_whole_marker_fit_rmsd_angstrom'])
            if own > whole + 1e-7 or abs((whole-own)-float(row['domain_fit_improvement_angstrom'])) > 1e-7:
                raise ValueError('Domain fit relation mismatch')
            local, pae = int(row['local_distance_pairs']), int(row['pae10_local_pairs'])
            if not 0 <= pae <= local <= n*(n-1)//2:
                raise ValueError('Residue-pair count mismatch')
            if local and abs(float(row['pae10_local_pair_fraction']) - pae/local) > 1e-12:
                raise ValueError('PAE fraction mismatch')
    if seen != expected or accepted_count != r['domain_comparisons'] or sum(exclusions.values()) != r['excluded_pair_or_domain_rows'] or len(accepted_pairs) != r['distinct_marker_taxon_pairs'] or len(families_seen) != r['distinct_pfam_domains']:
        raise ValueError('Incomplete result grid or summary mismatch')
    a.output.mkdir(parents=True)
    result = {'status': 'passed_complete_domain_disposition_grid_and_metric_readback',
              'domain_receipt_sha256': sha(a.domains / 'receipt.json'), 'script_sha256': sha(Path(__file__)),
              'baseline_pairs': len(baseline), 'candidate_dispositions': len(expected), 'accepted_domains': accepted_count,
              'exclusion_counts': dict(exclusions),
              'interpretation': 'Every candidate accounted for against reported eligible annotations; identities, bounds, baseline fields, fit relations and counts checked. Does not independently validate annotation eligibility, reconstruct residue masks or recompute coordinate geometry.'}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
