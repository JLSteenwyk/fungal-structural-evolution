#!/usr/bin/env python3
"""Compare frozen queue source identities; never authorize alignment reuse."""
import argparse,csv,hashlib,json
from collections import Counter
from pathlib import Path

def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(2**20),b''):h.update(block)
 return h.hexdigest()

def load(root):
 receipt=root/'receipt.json';r=json.loads(receipt.read_text());bindings={str(receipt):sha(receipt)}
 for name in ['models.jsonl','model_pairs.tsv']:
  p=root/name;bindings[str(p)]=sha(p)
  if bindings[str(p)]!=r['artifacts'][name]:raise ValueError('Changed queue artifact')
 models={}
 for line in (root/'models.jsonl').open():
  row=json.loads(line);key=(row['model_id'],row['version'])
  if key in models:raise ValueError('Repeated model')
  models[key]={k:row[k] for k in ['sha256','sequence_sha256','length']}
 pairs={}
 for row in csv.DictReader((root/'model_pairs.tsv').open(),delimiter='\t'):
  ends=tuple(sorted([(row['model_a'],int(row['version_a'])),(row['model_b'],int(row['version_b']))]))
  key=hashlib.sha256(json.dumps(ends,separators=(',',':')).encode()).hexdigest()
  if key!=row['pair_key'] or key in pairs or ends[0]==ends[1] or not all(e in models for e in ends):raise ValueError('Invalid pair')
  pairs[key]=ends
 if len(models)!=r['unique_models'] or len(pairs)!=r['unique_distinct_model_pairs']:raise ValueError('Queue count differs')
 return models,pairs,bindings

def classify(key,ends,oldmodels,oldpairs,newmodels):
 if key not in oldpairs:return 'new_pair_requires_alignment'
 if ends!=oldpairs[key]:raise ValueError('Pair endpoints changed under same key')
 if any(oldmodels[e]!=newmodels[e] for e in ends):return 'shared_pair_changed_source_requires_alignment'
 return 'shared_pair_matching_catalog_source_pending_input_and_result_checks'

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 for name in ['old','new','output']:ap.add_argument('--'+name,type=Path,required=True)
 args=ap.parse_args();om,op,ob=load(args.old);nm,np,nb=load(args.new)
 out=args.output;out.mkdir(parents=True,exist_ok=False);counts=Counter()
 with (out/'pair_reuse_screen.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['pair_key','disposition'])
  for key,ends in sorted(np.items()):
   status=classify(key,ends,om,op,nm);counts[status]+=1;w.writerow([key,status])
 bindings={**ob,**nb}
 for p,h in bindings.items():
  if sha(p)!=h:raise ValueError('Source changed during screen')
 result=dict(status='complete_catalog_source_reuse_screen_not_reuse_authorization',script_sha256=sha(__file__),source_hashes=bindings,old_pairs=len(op),new_pairs=len(np),old_pairs_absent=len(set(op)-set(np)),counts=dict(counts),artifacts={'pair_reuse_screen.tsv':sha(out/'pair_reuse_screen.tsv')},scope='All new pairs compared against frozen old model/version, catalog coordinate checksum, sequence checksum and length. Coordinate files are not rehashed here. Matching catalog identities are only reuse candidates: independently audited materialized PDB bytes, residue mapping, confidence mask, executable/options, input order and previous checkpoint results must match before any reuse. Original output and discrepancy flags must remain traceable. No alignment result copied or admitted.')
 (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
