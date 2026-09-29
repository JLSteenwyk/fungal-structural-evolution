"""Check every flagged background mapping with the analytic two-point RMSD."""
import csv,json,hashlib,math
from pathlib import Path

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
completion=Path('metadata/background_rmsd_diagnostic_completed_20260928.json')
c=json.loads(completion.read_text());receipt=Path(c['receipt'])
assert sha(receipt)==c['receipt_sha256']
r=json.loads(receipt.read_text())
plan=json.loads(Path('metadata/background_alignment_rmsd_diagnostic_plan_20260928.json').read_text())
source=json.loads(Path(plan['source_plan']).read_text());root=Path(source['output'])
assert sha(root/'receipt.json')==r['producer_receipt_sha256']
pr=json.loads((root/'receipt.json').read_text());manifest=root/'checkpoint_manifest.tsv'
assert sha(manifest)==pr['artifacts']['checkpoint_manifest.tsv']
proofs={v['path']:v['sha256'] for v in csv.DictReader(manifest.open(),delimiter='\t')}
checked=[]
for row in r['rmsd_discrepancies']:
 pair,mask,order=row['pair_key'],row['mask'],row['order']
 rel=f'pairs/{pair[:2]}/{pair}-{mask}-{order}.json';path=root/rel
 assert sha(path)==proofs[rel]
 record=json.loads(path.read_text());assert record['status']=='aligned'
 assert (record['pair_key'],record['mask'],record['order'])==(pair,mask,order)
 strings=[record['metrics']['alignment_left'],record['metrics']['alignment_right']]
 assert len(strings[0])==len(strings[1])
 paired=[i for i,(a,b) in enumerate(zip(*strings)) if a!='-' and b!='-']
 assert len(paired)==row['aligned_length']==2
 segments=[]
 for inp,text in zip(record['inputs'],strings):
  assert sha(inp['path'])==inp['sha256']
  atoms=[l for l in Path(inp['path']).read_text().splitlines() if l.startswith('ATOM  ')]
  assert all(l[12:16].strip()=='CA' and l[21]=='A' for l in atoms)
  coords=[[float(l[a:b]) for a,b in [(30,38),(38,46),(46,54)]] for l in atoms]
  assert len(coords)==len(text.replace('-',''))
  indexes=[sum(v!='-' for v in text[:i+1])-1 for i in paired]
  segments.append(math.dist(coords[indexes[0]],coords[indexes[1]]))
 analytic=abs(segments[0]-segments[1])/2
 assert math.isclose(analytic,row['rmsd_recomputed'],rel_tol=0,abs_tol=1e-10)
 assert abs(analytic-row['rmsd_native'])>.00501
 checked.append(dict(pair_key=pair,mask=mask,order=order,checkpoint_sha256=proofs[rel],segment_lengths=segments,analytic_rmsd=analytic,reconstructed_rmsd=row['rmsd_recomputed'],native_rmsd=row['rmsd_native']))
assert len(checked)==len(c['discrepancies'])
result=dict(status='verified_all_background_discrepancies_by_two_point_distance_formula',checker_sha256=sha(__file__),diagnostic_completion_sha256=sha(completion),diagnostic_receipt_sha256=sha(receipt),checked=checked,scope='Every flagged background mapping contains exactly two paired CA atoms. Optimal proper-rotation RMSD is half the absolute difference in segment lengths, independently confirming the reconstructed RMSD. Rotation is nonunique. Native rounding discrepancies remain flagged; no cause within native software or biological eligibility is inferred.')
path=Path('metadata/background_two_point_discrepancy_review_20260928.json')
with path.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print('All',len(checked),'flagged mappings confirmed by analytic two-point formula.')
