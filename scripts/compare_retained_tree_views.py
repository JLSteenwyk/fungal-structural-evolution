#!/usr/bin/env python3
"""Compare all 32 matched-taxon baseline views, retaining every split conflict."""
import argparse
import csv
import gzip
import itertools
import json
from pathlib import Path
from retained_tree_comparison_sources import load,PAIR_FIELDS,PRESENCE_FIELDS,CONFLICT_FIELDS,RULE
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def export(groups,output):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    comparison=[];presence=[];policy_summaries=[];conflict_count=0
    with gzip.open(output/'conflicts.tsv.gz','wt') as handle:
        writer=csv.DictWriter(handle,CONFLICT_FIELDS,delimiter='\t',lineterminator='\n');writer.writeheader()
        for policy,group in groups.items():
            views=group['views'];universe=group['taxa'];union=set().union(*(set(v['splits']) for v in views.values()))
            for split in sorted(union):
                for name,view in views.items():
                    row=view['splits'].get(split)
                    presence.append(dict(policy=policy,view=name,split_taxa_json=json.dumps(split),present=row is not None,
                        projected_branch_length_sum=row['projected_branch_length_sum'] if row else None,
                        empirical_projected_ufboot_percent=row['empirical_projected_ufboot_percent'] if row else None,sh_alrt_percent=None))
            pairs=[(a,b) for a,b in itertools.combinations(views,2) if views[a]['kind']==views[b]['kind'] or views[a]['run']==views[b]['run']]
            assert len(pairs)==16
            for a,b in pairs:
                left,right=views[a]['splits'],views[b]['splits'];common=set(left)&set(right);conflicts=high=0
                for x,y in itertools.product(sorted(set(left)-common),sorted(set(right)-common)):
                    xs,ys=set(x),set(y);cells=[xs&ys,xs-ys,ys-xs,universe-(xs|ys)]
                    if not all(cells):continue
                    conflicts+=1;sa=float(left[x]['empirical_projected_ufboot_percent']);sb=float(right[y]['empirical_projected_ufboot_percent']);supported=sa>=95 and sb>=95;high+=supported
                    writer.writerow(dict(policy=policy,view_a=a,view_b=b,split_a_taxa_json=json.dumps(x),split_b_taxa_json=json.dumps(y),
                        empirical_projected_ufboot_percent_a=sa,empirical_projected_ufboot_percent_b=sb,
                        both_UFB95=supported,witness_quartet=';'.join(min(cell) for cell in cells)))
                conflict_count+=conflicts;distance=len(left)+len(right)-2*len(common)
                comparison.append(dict(policy=policy,view_a=a,view_b=b,
                    comparison_scope='within_baseline_ML_vs_consensus' if views[a]['run']==views[b]['run'] else 'cross_baseline_'+views[a]['kind'],
                    support_rule=RULE,shared_internal_splits=len(common),unique_internal_splits_a=len(left)-len(common),unique_internal_splits_b=len(right)-len(common),
                    rf_distance=distance,normalized_rf=distance/(len(left)+len(right)),incompatible_split_pairs=conflicts,both_UFB95_incompatible_pairs=high))
            policy_summaries.append(dict(policy=policy,retained_taxa=len(universe),roles=group['roles'],views=8,comparisons=16,
                shared_all_ml=len(set.intersection(*(set(v['splits']) for v in views.values() if v['kind']=='ml'))),
                shared_all_consensus=len(set.intersection(*(set(v['splits']) for v in views.values() if v['kind']=='consensus'))),
                shared_all_eight=len(set.intersection(*(set(v['splits']) for v in views.values())))))
            print('compared_every_retained_baseline_view',policy,flush=True)
    for name,rows,fields in [('comparisons.tsv',comparison,PAIR_FIELDS),('split_presence.tsv',presence,PRESENCE_FIELDS)]:
        with (output/name).open('x') as handle:
            writer=csv.DictWriter(handle,fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
    return dict(policies=len(groups),tree_views=sum(len(g['views']) for g in groups.values()),comparison_rows=len(comparison),
                split_presence_rows=len(presence),incompatible_pairs=conflict_count,policy_summaries=policy_summaries)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);args=parser.parse_args()
    plan=json.loads(args.plan.read_text());groups,universe,bindings=load(plan,args.plan);output=Path(plan['output'])
    summary=export(groups,output);assert summary['policies']==4 and summary['tree_views']==32 and summary['comparison_rows']==64
    verify(bindings)
    receipt=dict(status='complete_full_retained_baseline_tree_comparisons_pending_readback',plan_sha256=sha(args.plan),
        **summary,source_hashes=bindings,artifacts={name:sha(output/name) for name in ['comparisons.tsv','split_presence.tsv','conflicts.tsv.gz']},
        scientific_eligibility=False,scope=plan['scope'])
    with (output/'receipt.json').open('x') as handle:handle.write(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='policy_summaries'}),flush=True)


if __name__=='__main__':main()
