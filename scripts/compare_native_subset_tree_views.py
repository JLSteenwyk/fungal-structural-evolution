#!/usr/bin/env python3
"""All 560 actual subset-fit comparisons with explicit unresolved consensus states."""
import argparse
import csv
import gzip
import itertools
import json
from pathlib import Path

from native_subset_tree_comparison_sources import (load,canonical,supported,pair_grid,comparison_scope,
    METRICS,PAIR_FIELDS,PRESENCE_FIELDS,CONFLICT_FIELDS,BOUNDARY_FIELDS)
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def write_table(path,rows,fields):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'xt') as f:
        writer=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)


def export(groups,root):
    root=Path(root);root.mkdir(parents=True,exist_ok=False)
    pairs=[];presence=[];boundaries=[];summaries=[];conflicts_count=0
    with gzip.open(root/'conflicts.tsv.gz','xt',compresslevel=6) as f:
        writer=csv.DictWriter(f,CONFLICT_FIELDS,delimiter='\t',lineterminator='\n');writer.writeheader()
        for cohort,g in sorted(groups.items()):
            taxa=g['taxa'];views=g['views'];tips=sorted(taxa);positions={t:i for i,t in enumerate(tips)};full=(1<<len(tips))-1
            union=sorted(set().union(*(set(v['splits']) for v in views.values())))
            masks={s:sum(1<<positions[t] for t in s) for s in union}
            for name,v in sorted(views.items()):
                for split in union:
                    metric=v['splits'].get(split)
                    presence.append(dict(cohort=cohort,view=name,family=v['family'],support_rule=v['support_rule'],
                        split_taxa_json=json.dumps(split),present=metric is not None,
                        **{k:metric[k] if metric else None for k in METRICS},support_rule_met=supported(metric,v['support_rule']) if metric else None))
                boundary=canonical(g['outgroups'],taxa);metric=v['splits'].get(boundary)
                boundaries.append(dict(cohort=cohort,view=name,ingroup=g['roles']['ingroup'],outgroup=g['roles']['outgroup'],
                    boundary_taxa_json=json.dumps(boundary),boundary_present=metric is not None,support_rule=v['support_rule'],
                    **{k:metric[k] if metric else None for k in METRICS},support_rule_met=supported(metric,v['support_rule']) if metric else None))
            for a,b in pair_grid(views):
                left,right=views[a]['splits'],views[b]['splits'];common=set(left)&set(right);count=high=0
                for x,y in itertools.product(sorted(set(left)-common),sorted(set(right)-common)):
                    xm,ym=masks[x],masks[y];cells=[xm&ym,xm&(~ym&full),ym&(~xm&full),full&~(xm|ym)]
                    if not all(cells):continue
                    sa,sb=supported(left[x],views[a]['support_rule']),supported(right[y],views[b]['support_rule'])
                    witness=';'.join(tips[(cell&-cell).bit_length()-1] for cell in cells)
                    writer.writerow(dict(cohort=cohort,view_a=a,view_b=b,split_a_taxa_json=json.dumps(x),split_b_taxa_json=json.dumps(y),
                        support_rule_a=views[a]['support_rule'],support_rule_b=views[b]['support_rule'],
                        support_rule_met_a=sa,support_rule_met_b=sb,both_support_rules_met=sa and sb,witness_quartet=witness))
                    count+=1;high+=sa and sb
                distance=len(left)+len(right)-2*len(common);total=len(left)+len(right);possible=len(taxa)-3
                pairs.append(dict(cohort=cohort,view_a=a,view_b=b,comparison_scope=comparison_scope(views[a],views[b]),
                    support_rule_a=views[a]['support_rule'],support_rule_b=views[b]['support_rule'],splits_a=len(left),splits_b=len(right),
                    possible_internal_splits=possible,shared_splits=len(common),unique_splits_a=len(left)-len(common),unique_splits_b=len(right)-len(common),
                    rf_distance=distance,normalized_rf=distance/(2*possible),normalized_rf_by_observed_splits=distance/total if total else None,
                    incompatible_pairs=count,both_support_rules_met_pairs=high))
                conflicts_count+=count
            categories={family:[set(v['splits']) for v in views.values() if v['family']==family] for family in ['native_pmsf','coalescent','concatenated']}
            summaries.append(dict(cohort=cohort,taxa=len(taxa),roles=g['roles'],native_views=len(categories['native_pmsf']),
                coalescent_views=len(categories['coalescent']),projected_reference_views=len(categories['concatenated']),
                unresolved_internal_slots_by_view={n:len(taxa)-3-len(v['splits']) for n,v in sorted(views.items())},
                shared_all_native=len(set.intersection(*categories['native_pmsf'])),
                shared_all_views=len(set.intersection(*(s for values in categories.values() for s in values)))))
            print('full_native_subset_cohort_compared',cohort,len(pair_grid(views)),flush=True)
    for name,rows,fields in [('comparisons.tsv',pairs,PAIR_FIELDS),('split_presence.tsv.gz',presence,PRESENCE_FIELDS),('role_boundaries.tsv',boundaries,BOUNDARY_FIELDS)]:write_table(root/name,rows,fields)
    return dict(cohorts=len(groups),tree_views=sum(len(g['views']) for g in groups.values()),comparison_rows=len(pairs),
        split_presence_rows=len(presence),incompatible_pairs=conflicts_count,role_boundary_rows=len(boundaries),
        boundary_views_with_role_split=sum(r['boundary_present'] for r in boundaries),cohort_summaries=summaries)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);args=parser.parse_args()
    plan=json.loads(args.plan.read_text());groups,universe,bindings=load(plan,args.plan);root=Path(plan['output']);summary=export(groups,root)
    assert summary['cohorts']==4 and summary['tree_views']==summary['role_boundary_rows']==88 and summary['comparison_rows']==560
    verify(bindings)
    result=dict(status='complete_full_native_subset_tree_comparisons_pending_independent_readback',plan_sha256=sha(args.plan),
        **summary,source_hashes=bindings,artifacts={p.name:sha(p) for p in root.iterdir() if p.is_file()},scientific_eligibility=False,scope=plan['scope'])
    with (root/'receipt.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if not isinstance(v,list)}),flush=True)


if __name__=='__main__':main()
