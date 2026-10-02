#!/usr/bin/env python3
"""Compare every coalescent candidate to same-taxon concatenated references."""
import argparse
import csv
import gzip
import itertools
import json
from pathlib import Path

from coalescent_tree_comparison_sources import (load,canonical,rule_met,pair_grid,comparison_scope,
    METRICS,PAIR_FIELDS,PRESENCE_FIELDS,CONFLICT_FIELDS,BOUNDARY_FIELDS)
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def table(path,rows,fields):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'wt') as handle:
        writer=csv.DictWriter(handle,fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)


def export(groups,output):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    pairs=[];presence=[];boundaries=[];summaries=[];conflicts_count=0
    with gzip.open(output/'conflicts.tsv.gz','wt') as handle:
        writer=csv.DictWriter(handle,CONFLICT_FIELDS,delimiter='\t',lineterminator='\n');writer.writeheader()
        for cohort,group in sorted(groups.items()):
            taxa=group['taxa'];views=group['views'];tips=sorted(taxa);index={t:i for i,t in enumerate(tips)}
            full_mask=(1<<len(tips))-1
            union=sorted(set().union(*(set(v['splits']) for v in views.values())))
            masks={s:sum(1<<index[t] for t in s) for s in union}
            for name,view in sorted(views.items()):
                for split in union:
                    row=view['splits'].get(split)
                    presence.append(dict(cohort=cohort,view=name,family=view['family'],support_rule=view['support_rule'],
                        split_taxa_json=json.dumps(split),present=row is not None,
                        **{field:row[field] if row else None for field in METRICS},
                        support_rule_met=rule_met(row,view['support_rule']) if row else None))
                boundary=canonical(group['outgroups'],taxa);row=view['splits'].get(boundary)
                boundaries.append(dict(cohort=cohort,view=name,ingroup=group['roles']['ingroup'],outgroup=group['roles']['outgroup'],
                    boundary_taxa_json=json.dumps(boundary),boundary_present=row is not None,support_rule=view['support_rule'],
                    **{field:row[field] if row else None for field in METRICS},
                    support_rule_met=rule_met(row,view['support_rule']) if row else None))
            for a,b in pair_grid(views):
                left,right=views[a]['splits'],views[b]['splits'];common=set(left)&set(right);count=high=0
                for x,y in itertools.product(sorted(set(left)-common),sorted(set(right)-common)):
                    am,bm=masks[x],masks[y]
                    cells=[am&bm,am&(~bm&full_mask),bm&(~am&full_mask),full_mask&~(am|bm)]
                    if not all(cells):continue
                    sa=rule_met(left[x],views[a]['support_rule']);sb=rule_met(right[y],views[b]['support_rule'])
                    witness=';'.join(tips[(cell&-cell).bit_length()-1] for cell in cells)
                    writer.writerow(dict(cohort=cohort,view_a=a,view_b=b,split_a_taxa_json=json.dumps(x),split_b_taxa_json=json.dumps(y),
                        support_rule_a=views[a]['support_rule'],support_rule_b=views[b]['support_rule'],
                        support_rule_met_a=sa,support_rule_met_b=sb,both_support_rules_met=sa and sb,witness_quartet=witness))
                    count+=1;high+=sa and sb
                distance=len(left)+len(right)-2*len(common)
                pairs.append(dict(cohort=cohort,view_a=a,view_b=b,comparison_scope=comparison_scope(views[a],views[b]),
                    support_rule_a=views[a]['support_rule'],support_rule_b=views[b]['support_rule'],
                    shared_splits=len(common),unique_splits_a=len(left)-len(common),unique_splits_b=len(right)-len(common),
                    rf_distance=distance,normalized_rf=distance/(len(left)+len(right)),incompatible_pairs=count,both_support_rules_met_pairs=high))
                conflicts_count+=count
            coal=[set(v['splits']) for v in views.values() if v['family']=='coalescent']
            refs=[set(v['splits']) for v in views.values() if v['family']=='concatenated']
            summaries.append(dict(cohort=cohort,taxa=len(taxa),roles=group['roles'],coalescent_views=len(coal),concatenated_views=len(refs),
                shared_all_coalescent=len(set.intersection(*coal)),shared_all_references=len(set.intersection(*refs)),
                shared_all_views=len(set.intersection(*(coal+refs)))))
            print('complete_coalescent_reference_cohort_compared',cohort,len(pair_grid(views)),flush=True)
    for name,rows,fields in [('comparisons.tsv',pairs,PAIR_FIELDS),('split_presence.tsv.gz',presence,PRESENCE_FIELDS),
                             ('role_boundaries.tsv',boundaries,BOUNDARY_FIELDS)]:table(output/name,rows,fields)
    return dict(cohorts=len(groups),tree_views=sum(len(g['views']) for g in groups.values()),comparison_rows=len(pairs),
        split_presence_rows=len(presence),incompatible_pairs=conflicts_count,role_boundary_rows=len(boundaries),
        boundary_views_with_role_split=sum(r['boundary_present'] for r in boundaries),cohort_summaries=summaries)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);args=parser.parse_args()
    plan=json.loads(args.plan.read_text());groups,universe,bindings=load(plan,args.plan)
    root=Path(plan['output']);summary=export(groups,root)
    if plan['mode']=='full_batch':assert summary['cohorts']==5 and summary['tree_views']==70 and summary['comparison_rows']==315
    else:assert summary['cohorts']==1 and summary['tree_views']==9 and summary['comparison_rows']==8
    verify(bindings)
    result=dict(status='complete_coalescent_reference_tree_comparisons_pending_independent_readback',
        plan_sha256=sha(args.plan),**summary,source_hashes=bindings,artifacts={p.name:sha(p) for p in root.iterdir() if p.is_file()},
        scientific_eligibility=False,scope=plan['scope'])
    with (root/'receipt.json').open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if not isinstance(v,list)}),flush=True)


if __name__=='__main__':main()
