"""Pair order may change; the gene-to-model assignment must not."""
from check_expanded_duplication_candidate_export import identity

row = dict(taxon_id='T', gene_a='a', gene_b='b', model_a='M1', model_b='M2', same_model='0')
reversed_row = dict(row, gene_a='b', gene_b='a', model_a='M2', model_b='M1')
assert identity(row) == identity(reversed_row)
assert identity(row) != identity(dict(row, model_a='M2', model_b='M1'))
assert identity(row) != identity(dict(row, gene_b='c'))
assert identity(row) != identity(dict(row, taxon_id='U'))
assert identity(row) != identity(dict(row, same_model='1'))
expanded = {('new_' + k if k.startswith('model_') or k == 'same_model' else k): v for k, v in reversed_row.items()}
assert identity(row) == identity(expanded, 'new_')
print('Pair-order invariance and gene/model/taxon/flag mismatch checks passed.')
