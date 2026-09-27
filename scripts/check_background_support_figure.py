#!/usr/bin/env python3
"""Read back every plotted point independently from verified taxon rows."""
import argparse,csv,json,math
from collections import defaultdict
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--receipt',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=json.loads(a.receipt.read_text())
for path,h in (r['source_bindings']|r['artifacts']).items():assert sha(path)==h
source=next(p for p in r['source_bindings'] if p.endswith('taxon_support.tsv'));groups=defaultdict(list)
for row in csv.DictReader(open(source),delimiter='\t'):
 if row['background_set']=='both_guides_unreported_parents':groups[(row['guide'],row['policy'])].append(row)
coverage=next(p for p in r['artifacts'] if p.endswith('_coverage.tsv'));seen=set()
for row in csv.DictReader(open(coverage),delimiter='\t'):
 key=(row['guide'],row['policy'],row['metric']);assert key not in seen;seen.add(key);rows=groups[key[:2]]
 n=sum(int(x[key[2]]) for x in rows);d=sum(int(x['targets']) for x in rows)
 assert int(row['numerator'])==n and int(row['denominator'])==d and math.isclose(float(row['percent']),100*n/d,abs_tol=1e-12)
assert seen=={(*key,metric) for key in groups for metric in ['supported_1_5','focal_supported_1_5']}
curves=next(p for p in r['artifacts'] if p.endswith('_concentration.tsv'));seen=set()
for row in csv.DictReader(open(curves),delimiter='\t'):
 key=row['guide'],int(row['rank']);assert key not in seen;seen.add(key)
 rows=sorted(groups[(key[0],'alignment_evalue')],key=lambda x:(-int(x['supported_1_5']),x['taxon_id']));rank=key[1];n=sum(int(x['supported_1_5']) for x in rows[:rank]);d=sum(int(x['supported_1_5']) for x in rows)
 assert row['taxon_id']==rows[rank-1]['taxon_id'] and int(row['supported_targets'])==int(rows[rank-1]['supported_1_5']) and int(row['cumulative_supported'])==n and int(row['total_supported'])==d and math.isclose(float(row['percent']),100*n/d,abs_tol=1e-12)
assert seen=={(guide,i+1) for (guide,policy),rows in groups.items() if policy=='alignment_evalue' for i in range(len(rows))}
a.output.write_text(json.dumps(dict(status='passed_all_background_support_figure_points',figure_receipt_sha256=sha(a.receipt),coverage_points=r['coverage_rows'],concentration_points=len(seen),scope='Every plotted numerator, denominator, percentage, ranked taxon and cumulative count independently reconstructed. Export hashes verified; no biological effect inference.'),indent=2)+'\n')
print(a.output.read_text())
