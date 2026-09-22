#!/usr/bin/env python3
"""Check repeat/order separation, no-hit handling and conservative exclusions."""
import json
from summarize_family_architectures import POLICIES, summarize, Variation


def row(taxon,seq,tokens,agree=1,partial=0):
    payload=dict(alternatives=[dict(ordered_model_tokens=tokens,alignment_overlap_pairs=0,
                                   partial_hmm_hits_below_070=partial)],
                 policies={p:dict(alternative_index=0,primary_rank_tie_pairs=0,
                                  unresolved_overlap_pairs=0,candidate_nested_pairs=0) for p in POLICIES})
    return taxon,seq,len(tokens),agree,json.dumps(payload)


if __name__=='__main__':
    a=['PF1','Domain']; b=['PF2','Repeat']
    rows=[row('F1','s1',[a,b]),row('F2','s2',[b,a]),row('F2','s3',[a,a,b]),
          row('F3','s4',[]),row('F4','s5',[a,b],partial=1),row('F5','s6',[a,b],agree=0)]
    for result in summarize(iter(rows)):
        expected=dict(proteins=6,taxa=5,unique_sequences=6,no_ga_hit_proteins=1,
                      taxa_with_multiple_family_members=1,maximum_members_per_taxon=2,
                      observed_proteins=5,observed_taxa=4,observed_ordered_signatures=3,
                      observed_multisets=2,observed_model_type_sets=1,
                      observed_multisets_with_order_variation=1,
                      observed_model_type_sets_with_multiplicity_variation=1,
                      observed_taxa_with_multiple_signatures=1,
                      conservative_proteins=3,conservative_taxa=2,
                      partial_hmm_proteins=1,policy_disagreement_proteins=1)
        for k,v in expected.items():
            assert result[k]==v,(k,result[k],v)
    for result in summarize([row('F1','s1',[]),row('F2','s2',[])]):
        assert result['observed_ordered_signatures']==0
        assert result['conservative_proteins']==0
    v=Variation()
    v.add('F1',[['PF1','Domain']])
    v.add('F2',[['PF1','Family']])
    assert v.result('')['model_type_sets']==2
    print('PASS: order versus multiplicity, taxon copies, no-hit exclusion, type preservation and uncertainty exclusions')
