#!/usr/bin/env python3
"""Compare production and independent graphs on reproducible synthetic cases."""
import random
from domain_clan_competition import compete
from readback_full_domain_competition import independent_policy
from check_domain_clan_competition import hit


def main():
    rng=random.Random(20260917)
    tests=0
    for case in range(500):
        hits=[]
        for i in range(rng.randrange(0,25)):
            start=rng.randrange(1,80)
            end=start+rng.randrange(1,60)
            hits.append(hit(str(i),start,end,rng.choice(['0','1e-1000','1e-40','1e-5']),
                            str(rng.randrange(10,150)),rng.choice(['','CL1','CL2']),
                            rng.choice(['PF1','PF2','PF3']),env=(max(1,start-2),end+3)))
        for coord in ['alignment','envelope']:
            for rank in ['evalue','bitscore']:
                assert compete(hits,{('PF1','PF2')},coord,rank)==independent_policy(hits,{('PF1','PF2')},coord,rank)
                tests+=1
    print(tests,'randomized complete policy comparisons passed')


if __name__=='__main__':
    main()
