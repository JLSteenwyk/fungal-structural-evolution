"""Independent fixed-code/IUPAC codon lookup from the pinned NCBI codebook."""
from itertools import product


DNA = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': 'AG', 'Y': 'CT', 'S': 'CG',
       'W': 'AT', 'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}


def lookup(codebook, code):
    letters = codebook[str(code)]['amino_acids']
    assert len(letters) == 64
    canonical = {''.join(codon): aa for codon, aa in zip(product('TCAG', repeat=3), letters)}
    table = {}
    ambiguities = {frozenset('DN'): 'B', frozenset('EQ'): 'Z', frozenset('IL'): 'J'}
    for symbols in product(DNA, repeat=3):
        amino_acids = {canonical[''.join(codon)] for codon in product(*(DNA[s] for s in symbols))}
        if len(amino_acids) == 1:
            amino_acid = next(iter(amino_acids))
        else:
            amino_acid = ambiguities.get(frozenset(amino_acids), 'X')
        table[''.join(symbols)] = amino_acid
    return table


def translate(dna, table):
    assert len(dna) % 3 == 0
    return ''.join(table[dna[i:i + 3].upper()] for i in range(0, len(dna), 3))
