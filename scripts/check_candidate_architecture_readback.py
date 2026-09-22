#!/usr/bin/env python3
"""Exercise the independent reader against seeded candidates and corrupted exports."""
from copy import deepcopy
import json
import random
from candidate_domain_architecture import POLICIES, describe
from readback_candidate_domain_architectures import check_payload


def main():
    rng = random.Random(20260922)
    for case in range(500):
        hits = []
        for i in range(rng.randrange(12)):
            start = rng.randrange(1, 100)
            end = start + rng.randrange(1, 80)
            hits.append(dict(hit_id=str(i), pfam_accession='PF' + str(rng.randrange(4)),
                             pfam_type=rng.choice(['Domain', 'Family', 'Repeat', 'Motif']),
                             pfam_clan='', alignment_start=start, alignment_end=end,
                             envelope_start=start, envelope_end=end + 1,
                             hmm_coverage=rng.choice(['0.5', '0.70', '1'])))
        policies = {key: dict(retained_hits=[h['hit_id'] for h in hits if rng.randrange(2)],
                              primary_rank_ties=[], unresolved_overlap_pairs=[], candidate_nested_pairs=[])
                    for key in POLICIES}
        payload = json.loads(json.dumps(describe(hits, policies)))
        check_payload(payload, hits, policies)
    hits = [dict(hit_id='a', pfam_accession='PF1', pfam_type='Repeat', pfam_clan='',
                 alignment_start=1, alignment_end=20, envelope_start=1, envelope_end=22, hmm_coverage='0.5')]
    policies = {key: dict(retained_hits=['a'], primary_rank_ties=[], unresolved_overlap_pairs=[], candidate_nested_pairs=[])
                for key in POLICIES}
    base = json.loads(json.dumps(describe(hits, policies)))
    mutations = [lambda p: p['alternatives'][0]['annotations'][0].update(pfam_type='Domain'),
                 lambda p: p['alternatives'][0].update(model_multiplicities=[['PF1', 'Repeat', 2]]),
                 lambda p: p['alternatives'][0].update(partial_hmm_hits_below_070=0),
                 lambda p: p['policies']['alignment_evalue'].update(primary_rank_tie_pairs=1),
                 lambda p: p['policies']['alignment_evalue'].update(alternative_index=5),
                 lambda p: p.update(status='no_GA_hit_not_proven_absence')]
    for mutate in mutations:
        bad = deepcopy(base)
        mutate(bad)
        try:
            check_payload(bad, hits, policies)
        except ValueError:
            continue
        raise AssertionError('Corrupted export accepted')
    print('Passed 500 seeded four-policy cases and six deliberate corruption checks')


if __name__ == '__main__':
    main()
