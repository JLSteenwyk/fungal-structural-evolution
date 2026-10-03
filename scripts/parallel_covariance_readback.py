"""Complete cohort-parallel replay using the frozen independent arithmetic."""
from collections import Counter
import itertools
import numpy as np

from covariance_basis_context import ComponentKernelProducts
from covariance_basis_independent import IndependentKernelProducts
from full_covariance_qualification_sources import cohort_rows, design_matrix, folded_operators, audit_id, covariance_status
from full_entity_operator_sources import MODES
from full_expanded_model_design_sources import AXES, DEGREES, digest
from full_expanded_model_input_sources import ORDERS
from readback_full_covariance_qualification import numeric


def verify_cohort(source,plan,cohort,designs,records):
    """Every original record is checked; only cohort scheduling is parallel."""
    assert len(designs)==30 and len(records)==30*len(MODES)*len(plan['trees'])
    rows=cohort_rows(source,cohort);contexts={};fallback={};seen=set()
    if len(rows):
        factors={tree:source['factors'][tree][rows] for tree in plan['trees']}
        for mode in MODES:
            operators=folded_operators(source,rows,mode)
            independent=IndependentKernelProducts(source['labels'][rows],operators)
            contexts[mode]={tree:independent.tree(factors[tree]) for tree in plan['trees']}
            fallback[mode]=(operators,factors)
    counts=Counter();index=[];reviews=0;maximum_error=0.;position=0
    for design in designs:
        matrix=design_matrix(source,cohort,rows,design)
        assert design['cohort_id']==cohort['cohort_id']
        seen.add((design['order_contrast'],design['sequence_axis'],design['degree']))
        for mode,tree in itertools.product(MODES,plan['trees']):
            record=records[position];position+=1;identifier=audit_id(source['contract'],design['design_id'],mode,tree)
            expected=dict(audit_id=identifier,source_contract=source['contract'],cohort_id=cohort['cohort_id'],
                design_id=design['design_id'],loading_mode=mode,tree=tree,residual_diagonal='uniform_one',
                records=len(rows),source_design_disposition=design['disposition'],active_column_indices=design['active_column_indices'],
                raw_design_sha256=design['raw_design_sha256'],folded_terms={
                'target_node':'combined_with_uniform_residual','family':'zero_signed_endpoint_kernel' if mode=='signed' else
                'combined_variance_family_intercept_plus_four_times_endpoint_family'})
            assert all(record[k]==v for k,v in expected.items())
            if design['disposition']!='full_rank_design':
                assert record['disposition']==design['disposition'] and record['numerical_audit'] is None
            elif record['disposition']=='numerical_covariance_qualification_requires_review':
                operators,factors=fallback[mode]
                try:ComponentKernelProducts(source['labels'][rows],operators,np.ones(len(rows))).tree(factors[tree]).audit(matrix)
                except (ValueError,ArithmeticError) as error:
                    assert record['error_type']==type(error).__name__ and record['error_message']==str(error)
                else:raise AssertionError('False numerical-failure disposition')
                assert record['numerical_audit'] is None
            else:
                reference=contexts[mode][tree].project(matrix);value=record['numerical_audit']
                assert value['fixed_effect_columns']==matrix.shape[1] and value['residual_dimension']==len(rows)-matrix.shape[1]
                assert value['family_components']==len(contexts[mode][tree].bank.parts)
                assert value['exactly_zero_incidence_names']==[k for k,z in fallback[mode][0].items() if not z.nnz]
                # This is the original frozen reciprocal-condition, complete
                # Gram/envelope, gesvd-rank and conservative-review validator.
                conservative,error=numeric(value,reference,contexts[mode][tree].bank.names,len(rows))
                reviews+=int(conservative);maximum_error=max(maximum_error,error)
                assert record['disposition']==covariance_status(value)
            counts[record['disposition']]+=1
            index.append([identifier,design['design_id'],mode,tree,record['disposition']])
    assert seen==set(itertools.product(ORDERS,AXES,DEGREES)) and position==len(records)
    return dict(cohort_id=cohort['cohort_id'],records=len(rows),design_rows=30,audit_rows=len(records),
        source_design_batch_sha256=digest(designs),source_audit_batch_sha256=digest(records),
        audit_status_counts=dict(counts),audit_index=index,conservative_independent_classification_differences=reviews,
        maximum_absolute_projected_gram_error=maximum_error,scientific_eligibility=False)
