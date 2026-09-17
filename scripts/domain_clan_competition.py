#!/usr/bin/env python3
"""Deterministic candidate clan competition; preserve exclusions and ambiguity.

This produces annotation hypotheses, not validated biological architectures.
Nested metadata plus strict coordinate containment exempts a competing pair,
but does not establish physical nesting or discontinuous match-state segments.
"""
from decimal import Decimal
from itertools import combinations


def overlaps(a, b, coordinates):
    return max(0, min(int(a[coordinates+'_end']), int(b[coordinates+'_end']))
               - max(int(a[coordinates+'_start']), int(b[coordinates+'_start'])) + 1)


def nesting_candidate(a, b, nested):
    def contains(outer, inner):
        return (outer['pfam_accession'], inner['pfam_accession']) in nested and (
            int(outer['alignment_start']) < int(inner['alignment_start'])
            and int(outer['alignment_end']) > int(inner['alignment_end']))
    return contains(a,b) or contains(b,a)


def compete(hits, nested, coordinates='alignment', ranking='evalue'):
    if coordinates not in ['alignment','envelope'] or ranking not in ['evalue','bitscore']:
        raise ValueError('Unknown competition policy')
    if len({h['hit_id'] for h in hits}) != len(hits):
        raise ValueError('Duplicate hit IDs within query')
    if len({h['sequence_id'] for h in hits}) > 1:
        raise ValueError('Competition must stay within one sequence')
    scores = {}
    for h in hits:
        e, b = Decimal(h['independent_evalue']), Decimal(h['domain_score'])
        if not e.is_finite() or e<0 or not b.is_finite():
            raise ValueError('Invalid score')
        if not (1<=int(h['envelope_start'])<=int(h['alignment_start'])
                <=int(h['alignment_end'])<=int(h['envelope_end'])):
            raise ValueError('Invalid coordinates')
        scores[h['hit_id']] = (e,-b) if ranking=='evalue' else (-b,e)

    def related(a,b):
        return bool(a['pfam_clan']) and a['pfam_clan']==b['pfam_clan']

    def conflicts(a,b):
        return related(a,b) and overlaps(a,b,coordinates)>0 and not nesting_candidate(a,b,nested)

    # Secondary score and ID provide a reproducible representative, while a tie
    # on the primary criterion is always reported as an ambiguity.
    ranked=sorted(hits,key=lambda h:(*scores[h['hit_id']],h['hit_id']))
    retained=[]
    decisions=[]
    for h in ranked:
        winners=[w for w in retained if conflicts(h,w)]
        if winners:
            decisions.append(dict(hit_id=h['hit_id'],disposition='suppressed_clan_competitor',
                                  competing_retained_hits=[w['hit_id'] for w in winners]))
        else:
            retained.append(h)
            decisions.append(dict(hit_id=h['hit_id'],disposition='retained_candidate',competing_retained_hits=[]))
    ties=[]
    for a,b in combinations(hits,2):
        if conflicts(a,b) and scores[a['hit_id']][0]==scores[b['hit_id']][0]:
            ties.append(sorted([a['hit_id'],b['hit_id']]))
    unresolved=[]
    nesting=[]
    for a,b in combinations(retained,2):
        if not overlaps(a,b,coordinates):
            continue
        pair=sorted([a['hit_id'],b['hit_id']])
        (nesting if nesting_candidate(a,b,nested) else unresolved).append(pair)
    ordered=sorted(retained,key=lambda h:(int(h['alignment_start']),int(h['alignment_end']),h['hit_id']))
    return dict(retained_hits=[h['hit_id'] for h in ordered],
                decisions=sorted(decisions,key=lambda d:d['hit_id']),
                primary_rank_ties=sorted(ties),unresolved_overlap_pairs=sorted(unresolved),
                candidate_nested_pairs=sorted(nesting),coordinates=coordinates,ranking=ranking)
