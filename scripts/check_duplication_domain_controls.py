#!/usr/bin/env python3
"""Check architecture categories and conservative repeat-aware domain matching."""
from inventory_duplication_domain_controls import architecture_class,single_copy_matches
A=('PF1.1','Domain');B=('PF2.1','Domain')
assert architecture_class([],[])=='neither_annotated'
assert architecture_class([A],[])=='one_unannotated'
assert architecture_class([A,B],[A,B])=='same_ordered_annotations'
assert architecture_class([A,B],[B,A])=='same_content_different_order'
assert architecture_class([A,A],[A])=='different_annotation_content'
assert architecture_class([A],[('PF1.2','Domain')])=='different_annotation_content'
def hit(name='PF1.1',kind='Domain',eligible=1,ident='a'):
    return dict(pfam_accession=name,pfam_type=kind,domain_interval_candidate=eligible,hit_id=ident)
x=hit();y=hit(ident='b')
assert single_copy_matches([x],[y])==[('PF1.1',x,y)]
assert not single_copy_matches([x,hit(eligible=0,ident='partial-repeat')],[y])
assert not single_copy_matches([hit(kind='Family')],[y])
assert not single_copy_matches([hit(eligible=0)],[y])
assert not single_copy_matches([x],[hit(name='PF1.2')])
assert not single_copy_matches([],[y])
print('Passed five architecture categories, repeat/version preservation, Family exclusion and ineligible-copy repeat blocking.')
