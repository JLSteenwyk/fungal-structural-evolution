#!/usr/bin/env python3
"""Check independent FASTA validation, including updated-content corruption."""
from pathlib import Path
from tempfile import TemporaryDirectory
from readback_full_codon_group_realignments import fasta,check_sequences
source={'a':'ACDX','b':'ADX'}
assert check_sequences(source,{'a':'ACDX','b':'A-DX'})==(2,4,7)
for bad in [{'b':'A-DX','a':'ACDX'},{'a':'ACDA','b':'A-DX'},{'a':'ACDX-','b':'A-DX-'},{'a':'ACDX','b':'ADX'},{'a':'ACDX'}]:
 try:check_sequences(source,bad)
 except AssertionError:pass
 else:raise AssertionError('Corruption accepted')
with TemporaryDirectory() as d:
 p=Path(d)/'x.faa'
 for content in ['>a\nACD\n>a\nACD\n','ACD\n>a\nACD\n']:
  p.write_text(content)
  try:fasta(p)
  except ValueError:pass
  else:raise AssertionError('Malformed FASTA accepted')
print('Passed independent alignment preservation and seven corruption cases')
