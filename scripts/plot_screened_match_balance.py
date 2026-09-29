#!/usr/bin/env python3
"""Plot audited strict-scenario balance across every whole-protein coverage screen."""
import argparse,csv,json,math,subprocess,time
from pathlib import Path
import numpy as np
import psutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_ortholog_pair_guide_comparison import sha
from readback_background_control_balance import FEATURES
LABELS=['Sequence distance','Log positive distance','Mean log length','Log length asymmetry','Mean pLDDT','Minimum pLDDT','Mean low-confidence fraction','Maximum low-confidence fraction']
METRICS=[('standardized_mean_difference','Target − control\n/ pooled matched SD'),('selection_shift_in_baseline_sd','Retained − original targets\n/ original-target SD'),('retention_shift_in_metadata_matched_sd','Retained − pre-screen matched targets\n/ pre-screen target SD')]
SCREENS=['n30_c50','n30_c70','n30_c90','n50_c50','n50_c70','n50_c90']
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();p=json.loads(a.plan.read_text());ph=sha(a.plan)
 def verify():
  assert sha(a.plan)==ph
  for path,h in p['pins'].items():assert sha(path)==h
 verify();dep=p['audit']
 while psutil.pid_exists(dep['pid']):
  try:
   proc=psutil.Process(dep['pid'])
   if abs(proc.create_time()-dep['created'])>.01 or proc.status()==psutil.STATUS_ZOMBIE:break
   assert proc.cmdline()==dep['cmdline']
  except psutil.NoSuchProcess:break
  time.sleep(30)
 state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',dep['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines());assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
 root=Path(p['source']);rp=root/'receipt.json';r=json.loads(rp.read_text());proof=json.loads(Path(p['proof']).read_text());assert proof['status']=='passed_full_screened_match_balance_readback' and proof['source_receipt_sha256']==sha(rp)
 bindings={str(rp):sha(rp),p['proof']:sha(p['proof'])}
 for name,h in r['artifacts'].items():assert sha(root/name)==h;bindings[str(root/name)]=h
 data={};coverage={}
 for filename,index,extra in [('balance.tsv',data,['feature']),('coverage.tsv',coverage,[])]:
  for row in csv.DictReader((root/filename).open(),delimiter='\t'):
   if row['scenario_id'] in ['S45','S46'] and row['policy']=='alignment_evalue' and row['mask']=='both_masks':
    key=tuple(row[k] for k in ['scenario_id','guide','screen']+extra);assert key not in index;index[key]=row
 assert set(data)=={(s,g,c,f) for s in ['S45','S46'] for g in ['profile','mafft'] for c in SCREENS for f in FEATURES}
 assert set(coverage)=={k[:3] for k in data}
 points=[]
 for (s,g,c,f),row in data.items():
  for field,_ in METRICS:
   assert row[field] and math.isfinite(float(row[field]))
   points.append(dict(scenario=s,guide=g,screen=c,feature=f,metric=field,value=row[field],pairs=row['pairs'],original_baseline_targets=row['baseline_targets'],matched_baseline_targets=row['metadata_matched_baseline_targets'],retained_matches=coverage[s,g,c]['retained_matches']))
 out=Path(p['output']);out.mkdir(exist_ok=False);table=out/'plotted_values.tsv'
 with table.open('w') as f:
  w=csv.DictWriter(f,list(points[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(points)
 # Re-read the exported values and compare every field to the source-derived records.
 assert list(csv.DictReader(table.open(),delimiter='\t'))==[{k:str(v) for k,v in row.items()} for row in points]
 plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'svg.fonttype':'none'})
 bounds={field:max(abs(float(row['value'])) for row in points if row['metric']==field) for field,_ in METRICS}
 outputs=[table]
 for scenario,description in [('S45','Any background taxon'),('S46','Focal taxon required')]:
  fig,axes=plt.subplots(2,3,figsize=(17,9),layout='constrained')
  for rownum,guide in enumerate(['profile','mafft']):
   for col,(field,title) in enumerate(METRICS):
    ax=axes[rownum,col];matrix=np.array([[float(data[scenario,guide,c,f][field]) for c in SCREENS] for f in FEATURES]);bound=bounds[field] or 1
    im=ax.imshow(matrix,cmap='RdBu_r',vmin=-bound,vmax=bound,aspect='auto');assert np.array_equal(np.asarray(im.get_array()),matrix)
    labels=[c.replace('_',' / ')+'\nN='+format(int(coverage[scenario,guide,c]['retained_matches']),',') for c in SCREENS]
    ax.set_xticks(range(6),labels,rotation=45,ha='right');ax.set_yticks(range(8),LABELS if col==0 else ['']*8);ax.set_title(guide+' · '+title,fontsize=10)
    for i in range(8):
     for j in range(6):ax.text(j,i,f'{matrix[i,j]:.2f}',ha='center',va='center',fontsize=7,color='white' if abs(matrix[i,j])>bound*.6 else '#222222')
    fig.colorbar(im,ax=ax,shrink=.7,label='Standardized difference')
  fig.suptitle(f'{scenario}: {description} — balance after structural screening\nBoth full and confidence-masked alignments must pass; alignment-E-value annotation policy',fontsize=13)
  fig.supxlabel('n = minimum aligned residues; c = minimum coverage (%) of each original protein; N = retained matches\nDescriptive diagnostics, not structural effects. Guides are dependent sensitivity alternatives. Color scales differ by metric, fixed across scenarios.\nPositive-log distance excludes zero-distance pairs; exact per-feature denominators are in plotted_values.tsv.',fontsize=9)
  for ext in ['png','pdf','svg']:
   path=out/(scenario+'.'+ext);fig.savefig(path,dpi=170);outputs.append(path)
  plt.close(fig)
 verify()
 for path,h in bindings.items():assert sha(path)==h
 receipt=dict(status='complete_screened_match_balance_figures_pending_visual_review',plan_sha256=ph,source_bindings=bindings,audit_terminal_state=state,plotted_values=len(points),color_bounds=bounds,artifacts={str(path):sha(path) for path in outputs},scope='Illustrative S45/S46, both guides, all six pre-existing screens; both masks required. No selection of primary test, effect estimate, confidence interval, phylogenetic correction or independence claim.')
 (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
