#!/usr/bin/env python3
"""Create a lossless full-atlas Foldseek-alias identity map."""
import argparse, csv, gzip, hashlib, json, os, tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()
def alias(model, version):
 return f'{model}-v{version}' if model and version else ''
def lookup_members(p):
 out=set()
 with open(p) as f:
  for n,line in enumerate(f,1):
   x=line.rstrip('\n').split('\t')
   if len(x)!=3 or not x[1] or x[1] in out: raise ValueError(f'invalid lookup row {n}')
   out.add(x[1])
 return out
def build(links, lookup, output):
 expected=lookup_members(lookup); seen=set(); counts=Counter(); identity_rows=0; taxa=set()
 output.parent.mkdir(parents=True,exist_ok=True)
 tmp=output.with_suffix(output.suffix+'.tmp')
 with open(links,newline='') as inp, gzip.open(tmp,'wt',newline='') as out:
  reader=csv.DictReader(inp,delimiter='\t'); fields=['foldseek_alias','taxon_id','protein_id','source','model_id','model_version','sequence_sha256','length','availability','paired_model_alias']
  out.write('\t'.join(fields)+'\n')
  for row in reader:
   af, es=alias(row['afdb_model_id'],row['afdb_version']),alias(row['esmfold_model_id'],row['esmfold_version'])
   for source, current, paired, mid, ver in [('AFDB',af,es,row['afdb_model_id'],row['afdb_version']),('ESMFold',es,af,row['esmfold_model_id'],row['esmfold_version'])]:
    if not current: continue
    if current not in expected: raise ValueError('source alias absent from lookup: '+current)
    if current not in seen:
     seen.add(current); counts[source]+=1
    identity_rows+=1; taxa.add(row['taxon_id'])
    out.write('\t'.join([current,row['taxon_id'],row['protein_id'],source,mid,ver,row['sequence_sha256'],row['length'],row['availability'],paired])+'\n')
 if seen!=expected: raise ValueError(f'identity map coverage differs: {len(seen)} of {len(expected)}')
 os.replace(tmp,output)
 return dict(models=len(seen),identity_rows=identity_rows,counts_by_source=dict(counts),taxa=len(taxa))
def test():
 with tempfile.TemporaryDirectory() as d:
  d=Path(d); l=d/'l'; q=d/'q'; o=d/'o.tsv.gz'
  q.write_text('0\tAF-X-v1\t0\n1\tES-Y-v1\t1\n'); l.write_text('taxon_id\tprotein_id\tsequence_sha256\tlength\tavailability\tafdb_model_id\tafdb_version\tafdb_path\tesmfold_model_id\tesmfold_version\tesmfold_path\tesmfold_prediction_config_sha256\nT\tP\th\t4\tboth\tAF-X\t1\t\tES-Y\t1\t\t\n')
  assert build(l,q,o)['models']==2
  q.write_text('0\tAF-X-v1\t0\n')
  try: build(l,q,d/'bad.tsv.gz')
  except ValueError: pass
  else: raise AssertionError('missing lookup alias accepted')
def main():
 p=argparse.ArgumentParser(description=__doc__); p.add_argument('--links',type=Path);p.add_argument('--lookup',type=Path);p.add_argument('--output',type=Path);p.add_argument('--receipt',type=Path);p.add_argument('--self-test',action='store_true');a=p.parse_args();test()
 if a.self_test: print(json.dumps({'status':'passed_full_atlas_identity_map_controls','checks':['complete_lookup_coverage','missing_alias_rejected']}));return
 if not all((a.links,a.lookup,a.output,a.receipt)): p.error('all paths required')
 if a.output.exists() or a.receipt.exists(): raise FileExistsError('fresh output and receipt required')
 result=build(a.links,a.lookup,a.output); result.update(status='completed_full_atlas_cluster_identity_map',checked_utc=datetime.now(timezone.utc).isoformat(),output=str(a.output),output_sha256=sha(a.output),links=str(a.links),links_sha256=sha(a.links),lookup=str(a.lookup),lookup_sha256=sha(a.lookup),scientific_eligibility=False,scope='Exact full lookup alias to source/taxon/representative-protein map with paired predictor aliases retained; no homology or evolutionary claim.')
 a.receipt.parent.mkdir(parents=True,exist_ok=True);a.receipt.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__': main()
