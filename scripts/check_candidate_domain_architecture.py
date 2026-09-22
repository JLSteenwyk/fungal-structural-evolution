#!/usr/bin/env python3
"""Check repeat, annotation-type, uncertainty and missing-hit preservation."""
from copy import deepcopy
from candidate_domain_architecture import POLICIES, describe


def policy(ids):
    return dict(retained_hits=ids, primary_rank_ties=[], unresolved_overlap_pairs=[], candidate_nested_pairs=[])


def main():
    hits = [dict(hit_id=str(i), pfam_accession=acc, pfam_type=kind, pfam_clan='',
                 alignment_start=start, alignment_end=end, envelope_start=start,
                 envelope_end=end, hmm_coverage=coverage)
            for i, acc, kind, start, end, coverage in [
                (1, 'PF1', 'Repeat', 1, 20, '1'), (2, 'PF1', 'Repeat', 30, 50, '0.69'),
                (3, 'PF2', 'Family', 40, 80, '1')]]
    states = {key: policy(['3', '2', '1']) for key in POLICIES}
    states['envelope_bitscore'] = policy(['1', '3'])
    states['alignment_evalue']['unresolved_overlap_pairs'] = [['2', '3']]
    result = describe(hits, states)
    assert len(result['alternatives']) == 2
    first = result['alternatives'][0]
    assert first['model_multiplicities'] == [['PF1', 'Repeat', 2], ['PF2', 'Family', 1]]
    assert first['ordered_model_tokens'] == [('PF1', 'Repeat'), ('PF1', 'Repeat'), ('PF2', 'Family')]
    assert first['alignment_overlap_pairs'] == 1 and first['partial_hmm_hits_below_070'] == 1
    assert result['policies']['alignment_evalue']['unresolved_overlap_pairs'] == 1
    assert describe([], {key: policy([]) for key in POLICIES})['status'] == 'no_GA_hit_not_proven_absence'
    bad = deepcopy(states)
    bad['alignment_evalue']['retained_hits'] = ['missing']
    try:
        describe(hits, bad)
    except ValueError:
        pass
    else:
        raise AssertionError('Accepted an unknown retained hit')
    print('Passed repeat/type preservation, overlap/partial ambiguity, alternatives, missing-hit and invalid-ID checks')


if __name__ == '__main__':
    main()
