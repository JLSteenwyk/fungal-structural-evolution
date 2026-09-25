#!/usr/bin/env python3
"""Real-coordinate equivalence, resume, and corrupted-checkpoint checks."""
import sys,tempfile,json,csv,shutil,contextlib,io
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,'scripts')
import prepare_paired_site_geometry as old
import prepare_checkpointed_paired_site_geometry as new
from Bio import SeqIO
base=Path('results/phylogeny/paired-inputs-esmfold-all-completed-20260922-v1')
snapshot=Path('results/structural_markers/esmfold-all-completed-20260922-v1');pae=Path('results/structural_pae/esmfold-all-completed-20260923-v1')
row=next(r for r in csv.DictReader((base/'marker_summary.tsv').open(),delimiter='\t') if r['status']=='ready_for_inference')
with tempfile.TemporaryDirectory(prefix='geometry-check-') as tmp:
 root=Path(tmp);inputs=root/'inputs';inputs.mkdir();folder=inputs/row['marker'];folder.mkdir()
 rows=list(SeqIO.parse(base/row['marker']/'aa.faa','fasta'));selected={r.id for r in rows[:2]}
 for name in ['aa.faa','3di.faa']:
  SeqIO.write([r for r in SeqIO.parse(base/row['marker']/name,'fasta') if r.id in selected],folder/name,'fasta')
 shutil.copy2(base/row['marker']/'columns.tsv',folder/'columns.tsv');row['eligible_taxa']='2'
 with (inputs/'marker_summary.tsv').open('w') as f:
  w=csv.DictWriter(f,row.keys(),delimiter='\t');w.writeheader();w.writerow(row)
 receipt=json.loads((base/'receipt.json').read_text());receipt['ready_markers']=1
 (inputs/'receipt.json').write_text(json.dumps(receipt))
 def checked(path):return json.loads((path/'receipt.json').read_text())
 def run(module,out):
  argv=['geometry','--inputs',str(inputs),'--snapshot',str(snapshot),'--pae',str(pae),'--output',str(out)]
  if module is new:argv+=['--checkpoints',str(root/'checkpoints')]
  with patch.object(sys,'argv',argv),patch.object(module,'checked_receipt',checked),contextlib.redirect_stdout(io.StringIO()):module.main()
 run(old,root/'old');run(new,root/'new');run(new,root/'resumed')
 for file in (root/'old').glob('*.tsv'):
  assert file.read_bytes()==(root/'new'/file.name).read_bytes()==(root/'resumed'/file.name).read_bytes()
 ck=root/'checkpoints'/(row['marker']+'.json');saved=json.loads(ck.read_text());saved['payload_sha256']='invalid';ck.write_text(json.dumps(saved))
 try:run(new,root/'corrupt')
 except ValueError as e:assert 'payload changed' in str(e)
 else:raise AssertionError('Corruption accepted')
 print('PASS: actual two-taxon marker output byte-identical to original, resume identical, corrupted checkpoint rejected. Fixture source receipts mocked only to restrict cohort; actual coordinates/PAE validated.')
