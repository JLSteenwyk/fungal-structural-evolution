#!/usr/bin/env python3
"""Check tie-set relations without treating a lexical choice as independent evidence."""
import json
from compare_duplication_sister_references import reference_relation

def row(genes,status='provisional_reference_available'):return dict(status=status,nearest_reference_genes=json.dumps(genes))
assert reference_relation(None,row(['a']))=='mafft_only_candidate'
assert reference_relation(row(['a']),None)=='profile_only_candidate'
assert reference_relation(row([],'no_modeled_nonfocal_sister'),row(['a']))=='not_provisional_in_both'
assert reference_relation(row(['a','b']),row(['b','a']))=='same_nearest_reference_set'
assert reference_relation(row(['a','b']),row(['b','c']))=='overlapping_nearest_reference_sets'
assert reference_relation(row(['a']),row(['b']))=='disjoint_nearest_reference_sets'
try:reference_relation(row([]),row(['b']))
except ValueError:pass
else:raise AssertionError('Empty provisional set accepted')
print('Six reference-set relations passed; invalid empty provisional set rejected.')
