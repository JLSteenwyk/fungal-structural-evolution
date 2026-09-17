#!/usr/bin/env python3
"""Exercise competition chains, nesting, ambiguity and decimal precision."""
from itertools import permutations
from domain_clan_competition import compete


def hit(name,start,end,e='1e-20',bits='100',clan='CL1',accession=None,env=None):
    return dict(hit_id=name,sequence_id='query',pfam_accession=accession or name,
                pfam_clan=clan,alignment_start=start,alignment_end=end,
                envelope_start=(env or (start,end))[0],envelope_end=(env or (start,end))[1],
                independent_evalue=e,domain_score=bits)


def main():
    # A-B-C chain must preserve nonoverlapping A/C; no component-wide winner.
    chain=[hit('a',1,10,'1e-30'),hit('b',8,18,'1e-20'),hit('c',16,25,'1e-10')]
    expected=compete(chain,set())
    assert expected['retained_hits']==['a','c']
    for order in permutations(chain):
        assert compete(order,set())==expected
    # Real decimal values below binary floating-point range remain distinct.
    precise=[hit('a',1,10,'1e-1000'),hit('b',1,10,'1e-900')]
    assert compete(precise,set())['retained_hits']==['a']
    assert not compete(precise,set())['primary_rank_ties']
    tied=[hit('a',1,10,'0','100'),hit('b',1,10,'0','200')]
    assert compete(tied,set())['primary_rank_ties']==[['a','b']]
    assert compete(tied,set())['retained_hits']==['b']
    # Cross-clan and missing-clan overlaps are unresolved, never discarded.
    for clan in ['CL2','']:
        result=compete([hit('a',1,10),hit('b',5,15,clan=clan)],set())
        assert len(result['retained_hits'])==2 and result['unresolved_overlap_pairs']==[['a','b']]
    pair=[hit('a',1,100),hit('b',20,40)]
    assert compete(pair,{('a','b')})['candidate_nested_pairs']==[['a','b']]
    assert len(compete(pair,{('b','a')})['retained_hits'])==1
    equal=[hit('a',1,100),hit('b',1,100)]
    assert not compete(equal,{('a','b')})['candidate_nested_pairs']
    env=[hit('a',2,10,'1e-30',env=(1,12)),hit('b',11,20,env=(10,22))]
    assert len(compete(env,set(),'alignment')['retained_hits'])==2
    assert len(compete(env,set(),'envelope')['retained_hits'])==1
    opposite=[hit('a',1,10,'1e-30','50'),hit('b',1,10,'1e-10','100')]
    assert compete(opposite,set(),ranking='evalue')['retained_hits']==['a']
    assert compete(opposite,set(),ranking='bitscore')['retained_hits']==['b']
    assert compete([],set())['retained_hits']==[]
    for invalid in [[hit('a',0,10)], [hit('a',1,10,e='NaN')], [hit('a',1,10),hit('a',2,12)]]:
        try:
            compete(invalid,set())
        except ValueError:
            pass
        else:
            raise AssertionError('Invalid input accepted')
    print('Passed: chain/order invariance, decimal precision, ties, clan ambiguity, directed strict nesting, coordinate/ranking sensitivities, empty and invalid inputs')


if __name__=='__main__':
    main()
