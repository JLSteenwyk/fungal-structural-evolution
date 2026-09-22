#!/usr/bin/env python3
"""Represent candidate annotation order without resolving overlaps or model types."""
from collections import Counter
from decimal import Decimal
import hashlib
import json

POLICIES = ('alignment_evalue', 'alignment_bitscore', 'envelope_evalue', 'envelope_bitscore')


def compact(value):
    return json.dumps(value, separators=(',', ':'), sort_keys=True)


def describe(hits, policies):
    by_id = {h['hit_id']: h for h in hits}
    if len(by_id) != len(hits) or set(policies) != set(POLICIES):
        raise ValueError('Duplicate hit or unexpected policy universe')
    alternatives, states, seen = [], {}, {}
    for policy in POLICIES:
        result = policies[policy]
        ids = result['retained_hits']
        if len(set(ids)) != len(ids) or not set(ids) <= set(by_id):
            raise ValueError('Invalid retained-hit identities')
        selected = sorted((by_id[i] for i in ids),
                          key=lambda h: (int(h['alignment_start']), int(h['alignment_end']), h['hit_id']))
        key = tuple(h['hit_id'] for h in selected)
        if key not in seen:
            annotations = []
            for h in selected:
                annotations.append({k: h[k] for k in ['hit_id', 'pfam_accession', 'pfam_type', 'pfam_clan',
                                    'alignment_start', 'alignment_end', 'envelope_start', 'envelope_end', 'hmm_coverage']})
            tokens = [(h['pfam_accession'], h['pfam_type']) for h in selected]
            counts = Counter(tokens)
            overlap_pairs = sum(int(a['alignment_end']) >= int(b['alignment_start'])
                                for i, a in enumerate(selected) for b in selected[i + 1:])
            alternatives.append(dict(annotations=annotations,
                                     ordered_model_tokens=tokens,
                                     model_multiplicities=[[acc, kind, n] for (acc, kind), n in sorted(counts.items())],
                                     token_signature_sha256=hashlib.sha256(compact(tokens).encode()).hexdigest(),
                                     alignment_overlap_pairs=overlap_pairs,
                                     partial_hmm_hits_below_070=sum(Decimal(h['hmm_coverage']) < Decimal('0.70') for h in selected),
                                     interpretation='Coordinate-sorted candidate annotations; repeated models and all Pfam types retained. Signature is not a homology assignment or validated biological architecture.'))
            seen[key] = len(alternatives) - 1
        states[policy] = dict(alternative_index=seen[key],
                              primary_rank_tie_pairs=len(result['primary_rank_ties']),
                              unresolved_overlap_pairs=len(result['unresolved_overlap_pairs']),
                              candidate_nested_pairs=len(result['candidate_nested_pairs']))
    return dict(alternatives=alternatives, policies=states,
                status='candidate_architecture_requires_validation' if hits else 'no_GA_hit_not_proven_absence')
