#!/usr/bin/env python3
"""Audit a complete real bootstrap shard, then reject cost and identity corruption."""
import csv,gzip,json,shutil,tempfile,time
from pathlib import Path
from ecology_bootstrap_flow_readback import audit_tree,sha
p=json.loads(Path('metadata/ecology_bootstrap_edges_plan_20260927.json').read_text());spec=p['trees'][0]
with open(p['evidence']) as f:ev=list(csv.DictReader(f,delimiter='\t'))
states={x['taxon_id']:1 if x['state']=='ectomycorrhizal' else 0 for x in ev if x['state'] in ['ectomycorrhizal','saprotrophic','asymbiotic']}
with open(p['manifest']) as f:taxa=[x['taxon_id'] for x in csv.DictReader(f,delimiter='\t')]
line=Path(spec['bootstrap']).read_text().splitlines()[0];name=spec['label']+'-0001.json.gz';path=Path(p['output'])/name
start=time.monotonic();proof=audit_tree((spec['label'],1,line,states,taxa,p['output'],sha(path)));elapsed=time.monotonic()-start
with gzip.open(path,'rt') as f:original=json.load(f)
with tempfile.TemporaryDirectory() as d:
 for mode in ['wrong_cost','duplicate_edge']:
  data=json.loads(json.dumps(original))
  if mode=='wrong_cost':
   values=data[0]['edges'][0][4];i=next(i for i,v in enumerate(values) if v is not None);values[i]+=1
  else:data[0]['edges'].insert(1,data[0]['edges'][0])
  q=Path(d)/name
  with gzip.open(q,'wt') as f:json.dump(data,f)
  try:audit_tree((spec['label'],1,line,states,taxa,d,sha(q)))
  except AssertionError:pass
  else:raise AssertionError('Corruption accepted: '+mode)
r=dict(status='passed_real_bootstrap_shard_and_corruption_checks',constrained_costs=proof['constrained_costs'],elapsed_seconds=elapsed,source_shard_sha256=sha(path),helper_sha256=sha('scripts/ecology_bootstrap_flow_readback.py'),fixture_sha256=sha(__file__),rejected=['wrong_cost_with_updated_artifact_hash','duplicate_edge_with_updated_artifact_hash'])
Path('metadata/ecology_bootstrap_flow_fixture_20260927.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
