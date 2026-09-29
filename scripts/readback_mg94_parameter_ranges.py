#!/usr/bin/env python3
"""Independently parse all saved codon fits and reconstruct every parameter range."""
import csv,hashlib,json,math,re,subprocess
from pathlib import Path

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def parse(text):
 result={}
 for line in text.splitlines():
  match=re.match(r'^(?:global )?([A-Za-z0-9_.]+)(:=|=)([^;]+);',line)
  if not match:continue
  name,operator,raw=match.groups()
  relevant=(line.startswith('global ') and '.model_MGREV.' in name and operator=='=') or ('.tree_0.' in name and name.endswith('.t'))
  if not relevant:continue
  assert operator=='=',('fixed branch',name)
  value=float(raw);assert math.isfinite(value) and value>=0 and name not in result
  result[name]=value
 assert result
 return result

def main():
 root=Path('results/cds/local-mg94-unconstrained-multistarts-20260927-v1');out=Path('results/cds/local-mg94-multistart-parameter-ranges-20260928-v1');audit=Path('results/cds/local-mg94-unconstrained-multistart-audit-20260927-v1')
 bindings={}
 def load(path):bindings[str(path)]=sha(path);return json.loads(path.read_text())
 r=load(root/'receipt.json');q=load(out/'receipt.json');a=load(audit/'receipt.json');assert q['source_receipt_sha256']==a['source_receipt_sha256']==sha(root/'receipt.json') and q['source_audit_sha256']==sha(audit/'receipt.json')
 assert a['status']=='passed_full_local_mg94_multistart_artifact_and_numeric_audit' and q['status']=='complete_full_codon_multistart_parameter_ranges'
 states={}
 for unit in ['fungal-local-mg94-unconstrained-multistarts-20260927.service','fungal-local-mg94-unconstrained-multistart-audit-20260927.service','fungal-mg94-parameter-ranges-20260928.service']:
  state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',unit,'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines());assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0');states[unit]=state
 path=out/'all_parameter_ranges.tsv';assert sha(path)==q['artifacts'][path.name];bindings[str(path)]=sha(path)
 table={}
 for row in csv.DictReader(path.open(),delimiter='\t'):
  key=row['case_id'],row['parameter'];assert key not in table;table[key]=row
 seen=set();fits=0;cases=set();parameter_values=0
 for entry in r['case_receipts']:
  case=entry['case_id'];assert case not in cases;cases.add(case);cp=root/case/'receipt.json';assert sha(cp)==entry['receipt_sha256'];rec=load(cp)
  proofs={x['start_label']:x for x in rec['proofs']};rows=rec['rows'];assert len(rows)==len(proofs)==8 and len({x['start_label'] for x in rows})==8
  best=sorted(float(x['log_likelihood']) for x in rows)[-1];values={};near={};near_count=0
  for row in rows:
   fit=root/case/row['start_label']/'fit.bf';assert sha(fit)==proofs[row['start_label']]['artifacts']['fit.bf'];params=parse(fit.read_text());bindings[str(fit)]=sha(fit);fits+=1;parameter_values+=len(params)
   if values:assert set(params)==set(values)
   selected=best-float(row['log_likelihood'])<=q['likelihood_tolerance'];near_count+=int(selected)
   for name,v in params.items():
    values.setdefault(name,[]).append(v)
    if selected:near.setdefault(name,[]).append(v)
  assert 1<=near_count<=8 and set(values)==set(near)
  for name,allvalues in values.items():
   key=case,name;seen.add(key);actual=table[key];ordered=sorted(allvalues);selected=sorted(near[name]);expected=dict(all_start_min=ordered[0],all_start_max=ordered[-1],near_best_min=selected[0],near_best_max=selected[-1],near_best_range=selected[-1]-selected[0])
   assert int(actual['near_best_starts'])==near_count
   for field,value in expected.items():assert math.isclose(float(actual[field]),value,rel_tol=1e-12,abs_tol=1e-12),(key,field)
  if len(cases)%200==0:print('Verified parameter ranges',len(cases),'/1632',flush=True)
 assert seen==set(table) and len(cases)==q['cases']==r['cases']==1632 and fits==r['fits']==13056 and len(seen)==q['parameter_rows']==28094 and parameter_values==a['parameters_checked']
 for path,h in bindings.items():assert sha(path)==h,path
 result=dict(status='passed_full_independent_mg94_parameter_range_readback',cases=len(cases),fits=fits,parameter_values=parameter_values,range_rows=len(seen),source_receipt_sha256=sha(root/'receipt.json'),range_receipt_sha256=sha(out/'receipt.json'),audit_receipt_sha256=sha(audit/'receipt.json'),checker_sha256=sha(__file__),terminal_states=states,scope='Every saved fit hash and parameter parsed with separate parser; all-start and near-best extrema, spans, counts, identities and full membership independently reconstructed. Descriptive numerical ranges, not confidence intervals, global optima or selection evidence.')
 with Path('metadata/mg94_parameter_ranges_completed_readback_20260928.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps(result),flush=True)
if __name__=='__main__':main()
