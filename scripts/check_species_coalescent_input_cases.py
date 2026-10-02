#!/usr/bin/env python3
"""Check contraction-before-pruning and missing/exact-cutoff software cases."""
import argparse
from collections import Counter
import io
import json
from pathlib import Path

from Bio import Phylo
from prepare_species_coalescent_inputs import transform, topology
from readback_species_coalescent_inputs import parse, expected
from run_ortholog_pair_guide_comparison import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    text = '(((A:.1,B:.1)10:.2,C:.1)9.9:.2,(D:.1,E:.1)80:.2,(F:.1,(G:.1,H:.1):.2)79.9:.2);'
    index = {taxon: i for i, taxon in enumerate('ABCDEFGH')}
    _, source_mask, source = parse(text, index)
    rows = [dict(side=frozenset(taxon for taxon, bit in index.items() if mask & (1 << bit)), support=value)
            for mask, value in source.items()]
    cases = []
    cohorts = [set('ABCDEFGH'), set('BCDEFGH'), set('ABCDEF'), set('ACDEFG'), set('ABDEGH')]
    for i, retained in enumerate(cohorts):
        for cutoff in [None, 10, 80]:
            tree, counts = transform(Phylo.read(io.StringIO(text), 'newick'), retained, cutoff)
            produced = topology(tree)
            tips, mask, splits = parse(produced, index, clean=True)
            wanted, reconstructed = expected(rows, mask, cutoff, index)
            assert set(tips) == retained and set(splits) == wanted and counts == reconstructed
            assert counts['original_internal'] == 5
            assert counts['contracted_missing_support'] == (0 if cutoff is None else 1)
            assert counts['contracted_below_support'] == (0 if cutoff is None else 1 if cutoff == 10 else 3)
            cases.append(dict(cohort=i, cutoff=cutoff, taxa=len(tips), internal_splits=len(splits), text=produced, counts=counts))
    rejected = []
    for i, case in enumerate(cases):
        modified = case['text'].replace(next(t for t in 'ABCDEFGH' if t in case['text']), 'Z', 1)
        try:
            parse(modified, index, clean=True)
        except AssertionError: rejected.append('foreign_taxon_' + str(i))
        else: raise AssertionError('Accepted foreign taxon')
    for cutoff, changed in [(10, 10.01), (80, 80.01)]:
        tree, _ = transform(Phylo.read(io.StringIO(text), 'newick'), cohorts[0], changed)
        _, mask, actual = parse(topology(tree), index, clean=True)
        wanted, _ = expected(rows, mask, cutoff, index)
        assert set(actual) != wanted
        rejected.append('incorrect_exact_cutoff_' + str(cutoff))
    result = dict(status='passed_full_five_cohort_three_policy_contraction_pruning_cases',
                  cases=cases, false_inputs_rejected=len(rejected), rejected=rejected,
                  script_sha256=sha(__file__), producer_sha256=sha('scripts/prepare_species_coalescent_inputs.py'),
                  reader_sha256=sha('scripts/readback_species_coalescent_inputs.py'), scientific_eligibility=False,
                  scope='Synthetic eight-tip tests of all five cohort/three support settings, missing support and exact SH-aLRT cutoffs. Independent DendroPy split reconstruction versus Bio contraction/pruning. Not a biological pilot or real input completion proof.')
    with (args.output / 'receipt.json').open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(status=result['status'], cases=len(cases), false_inputs_rejected=len(rejected))), flush=True)


if __name__ == '__main__':
    main()
