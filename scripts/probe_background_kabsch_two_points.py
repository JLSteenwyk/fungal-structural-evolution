#!/usr/bin/env python3
"""Probe the local Kabsch routine on all flagged background two-point mappings."""
import csv,json,math,subprocess
from pathlib import Path
from screen_duplication_alignment_reuse import sha
source=Path('/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/US-align/USalign.cpp')
review=Path('metadata/background_two_point_discrepancy_review_20260928.json');proof=json.loads(review.read_text())
root=Path('results/structural_comparisons/background-alignments-20260927-v1')
out=Path('results/structural_comparisons/background-kabsch-two-point-probe-20260928-v1');out.mkdir(exist_ok=False)
text=source.read_text();start=text.index('bool Kabsch(');brace=text.index('{',start);depth=1;end=brace+1
while depth:
 depth+=(text[end]=='{')-(text[end]=='}');end+=1
routine=text[start:end]
harness='#include <cmath>\n#include <iostream>\n#include <iomanip>\nusing namespace std;\n'+routine+'''
int main(){double a[2][3],b[2][3];double* x[2]={a[0],a[1]};double* y[2]={b[0],b[1]};
while(cin>>a[0][0]){for(int i=1;i<6;i++)cin>>a[i/3][i%3];for(int i=0;i<6;i++)cin>>b[i/3][i%3];
double rms,t[3],u[3][3];Kabsch(x,y,2,0,&rms,t,u);cout<<setprecision(17)<<sqrt(rms/2)<<"\\n";} }
'''
(out/'probe.cpp').write_text(harness)
command=['g++','-O2','-o',str(out/'probe'),str(out/'probe.cpp')];subprocess.run(command,check=True)
rows=[];inputs=[]
for item in proof['checked']:
 pair,mask,order=item['pair_key'],item['mask'],item['order'];path=root/f'pairs/{pair[:2]}/{pair}-{mask}-{order}.json';assert sha(path)==item['checkpoint_sha256']
 record=json.loads(path.read_text());texts=[record['metrics']['alignment_left'],record['metrics']['alignment_right']];positions=[i for i,(a,b) in enumerate(zip(*texts)) if a!='-' and b!='-'];assert len(positions)==2
 coords=[]
 for inp,alignment in zip(record['inputs'],texts):
  assert sha(inp['path'])==inp['sha256']
  atoms=[l for l in Path(inp['path']).read_text().splitlines() if l.startswith('ATOM  ')]
  xyz=[[float(l[a:b]) for a,b in [(30,38),(38,46),(46,54)]] for l in atoms]
  coords.append([xyz[sum(c!='-' for c in alignment[:i+1])-1] for i in positions])
 analytic=abs(math.dist(*coords[0])-math.dist(*coords[1]))/2;assert abs(analytic-item['analytic_rmsd'])<1e-10
 for variant in ['original','centered','translated_1000']:
  selected=[]
  for side in coords:
   center=[(side[0][k]+side[1][k])/2 for k in range(3)]
   selected.append([[v-(center[k] if variant=='centered' else 0)+(1000 if variant=='translated_1000' else 0) for k,v in enumerate(atom)] for atom in side])
  inputs.append(' '.join(format(v,'.17g') for side in selected for atom in side for v in atom))
  rows.append(dict(pair_key=pair,mask=mask,order=order,variant=variant,analytic_rmsd=analytic,original_native_rmsd=item['native_rmsd']))
input_text='\n'.join(inputs)+'\n';(out/'coordinates.txt').write_text(input_text)
run=subprocess.run([str(out/'probe')],input=input_text,text=True,capture_output=True,check=True);(out/'stdout.txt').write_text(run.stdout)
values=[float(x) for x in run.stdout.splitlines()];assert len(values)==len(rows)
for row,value in zip(rows,values):
 assert math.isfinite(value);row.update(probe_rmsd=value,error_from_analytic=abs(value-row['analytic_rmsd']))
with (out/'probe_results.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
result=dict(status='complete_isolated_local_source_kabsch_probe',source_sha256=sha(source),review_sha256=sha(review),script_sha256=sha(__file__),compiler=subprocess.check_output(['g++','--version'],text=True).splitlines()[0],compile_command=command,cases=len(rows),artifacts={name:sha(out/name) for name in ['probe.cpp','probe','coordinates.txt','stdout.txt','probe_results.tsv']},maximum_error_by_variant={v:max(r['error_from_analytic'] for r in rows if r['variant']==v) for v in ['original','centered','translated_1000']},scope='Exact local Kabsch function extracted and compiled separately at O2; original, centered and translated two-point coordinates tested. Translation preserves analytic RMSD mathematically. This is not a reproduction of the full native alignment pipeline, which calls Kabsch on transformed coordinates; source-to-installed-binary build identity is unproven. No production executable or results modified. RMSD discrepancy flags remain.')
(out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
