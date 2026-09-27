#!/usr/bin/env python3
"""Exercise exact sequence/order preservation and rejection of malformed outputs."""
from pathlib import Path
from tempfile import TemporaryDirectory
from run_full_codon_group_realignments import validate
with TemporaryDirectory() as d:
 p=Path(d);source=p/'input.faa';target=p/'output.faa';source.write_text('>a\nACDX\n>b\nADX\n')
 target.write_text('>a\nACDX\n>b\nA-DX\n');assert validate(source,target)==(2,4)
 for bad in ['>a\nACDX\n>a\nA-DX\n','>b\nA-DX\n>a\nACDX\n','>a\nACDA\n>b\nA-DX\n','>a\nACDX-\n>b\nA-DX-\n','>a\nACDX\n>b\nADX\n']:
  target.write_text(bad)
  try:validate(source,target)
  except AssertionError:pass
  else:raise AssertionError('Invalid alignment accepted')
print('Passed exact gapped sequence preservation and five malformed alignment rejections')
