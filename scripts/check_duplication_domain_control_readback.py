#!/usr/bin/env python3
"""Known repeat-aware reconstruction cases plus broad producer/readback agreement fixtures."""
import random
from readback_duplication_domain_controls import pair_summary
from inventory_duplication_domain_controls import architecture_class,single_copy_matches

def hit(a='PF1.1',eligible=1,ident='x',kind='Domain'):
    return dict(pfam_accession=a,pfam_type=kind,hit_id=ident,domain_interval_candidate=eligible,conservative_architecture=eligible)
a=hit();b=hit(ident='b')
r,m=pair_summary([a],[b]);assert r['annotation_class']=='same_ordered_annotations' and len(m)==1
r,m=pair_summary([a,hit(eligible=0,ident='repeat')],[b]);assert r['annotation_class']=='different_annotation_content' and not m
r,m=pair_summary([],[]);assert r['annotation_class']=='neither_annotated' and r['both_conservative']==0
rng=random.Random(260926)
for case in range(1000):
    rows=[[hit(a='PF'+str(rng.randrange(5))+'.1',eligible=rng.randrange(2),ident=str(i),kind=rng.choice(['Domain','Family'])) for i in range(rng.randrange(8))] for _ in range(2)]
    result,matches=pair_summary(*rows)
    assert result['annotation_class']==architecture_class(*[[(h['pfam_accession'],h['pfam_type']) for h in hits] for hits in rows])
    assert matches==single_copy_matches(*rows)
print('Known empty/repeat cases and 1000 independently reconstructed category/domain-match fixtures passed.')
