#!/usr/bin/env python3
"""Complete 32-view/64-pair synthetic comparison and counterfeit export checks."""
import argparse
import csv
import gzip
import io
import json
from pathlib import Path
import shutil
from Bio import Phylo
from project_pmsf_retained_taxa import raw_edges,project
from compare_retained_tree_views import export
from readback_retained_tree_comparisons import inspect
from retained_tree_comparison_sources import PAIR_FIELDS,PRESENCE_FIELDS,CONFLICT_FIELDS
from run_ortholog_pair_guide_comparison import sha


def rewrite(path,fields,rows,compressed=False):
    opener=gzip.open if compressed else open
    with opener(path,'wt') as handle:
        writer=csv.DictWriter(handle,fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    root=args.output;root.mkdir(parents=True,exist_ok=False)
    a='(((A:0.1,B:0.1):0.2,X:0.1):0.3,(C:0.1,D:0.1):0.4,(E:0.1,Y:0.1):0.5);'
    b='(((A:0.1,X:0.1):0.2,B:0.1):0.3,(C:0.1,E:0.1):0.4,(D:0.1,Y:0.1):0.5);'
    universe=set('ABXCDEY');groups={}
    policies=[('drop_X','ABCDEY'),('drop_X_Y','ABCDE'),('drop_X_A','BCDEY'),('drop_X_B_Y','ACDE')]
    source_runs={}
    for label,ml,consensus in [('pp',a,a),('pm',b,a),('mp',a,b),('mm',b,b)]:
        folder=root/'synthetic-baselines'/label;folder.mkdir(parents=True)
        (folder/'pmsf.treefile').write_text(ml+'\n');(folder/'pmsf.contree').write_text(consensus+'\n');source_runs[label]=folder
    for policy,tips in policies:
        retained=set(tips);outgroups=retained & set('EY');views={}
        for i,(label,folder) in enumerate(source_runs.items()):
            for kind,file in [('ml','pmsf.treefile'),('consensus','pmsf.contree')]:
                projected=project(raw_edges(Phylo.read(folder/file,'newick'),universe),retained)
                # Thresholds exercise exact95, below95 and above95 support screens.
                support=['95.0','94.9','99.0','95.0'][i]
                splits={key:dict(projected_branch_length_sum=str(length),empirical_projected_ufboot_percent=support,sh_alrt_percent='')
                    for key,(length,n) in projected.items() if len(key)>1}
                views[label+':'+kind]=dict(run=label,kind=kind,splits=splits,source=dict(run=str(folder)))
        groups[policy]=dict(taxa=retained,outgroups=outgroups,roles=dict(ingroup=len(retained-outgroups),outgroup=len(outgroups)),views=views)
    valid=root/'valid';summary=export(groups,valid);proof=inspect(groups,universe,valid)
    assert summary==proof and summary['comparison_rows']==64 and summary['tree_views']==32 and summary['incompatible_pairs']>0
    mutations=[('wrong_RF','comparisons.tsv',PAIR_FIELDS,'rf_distance','999'),
        ('wrong_normalized_RF','comparisons.tsv',PAIR_FIELDS,'normalized_rf','1.5'),
        ('wrong_shared_count','comparisons.tsv',PAIR_FIELDS,'shared_internal_splits','999'),
        ('wrong_unique_count','comparisons.tsv',PAIR_FIELDS,'unique_internal_splits_a','999'),
        ('wrong_conflict_count','comparisons.tsv',PAIR_FIELDS,'incompatible_split_pairs','999'),
        ('wrong_supported_conflict_count','comparisons.tsv',PAIR_FIELDS,'both_UFB95_incompatible_pairs','999'),
        ('invented_SH_rule','comparisons.tsv',PAIR_FIELDS,'support_rule','SH80_and_UFB95'),
        ('wrong_scope','comparisons.tsv',PAIR_FIELDS,'comparison_scope','native_refit'),
        ('wrong_presence','split_presence.tsv',PRESENCE_FIELDS,'present','False'),
        ('wrong_presence_path','split_presence.tsv',PRESENCE_FIELDS,'projected_branch_length_sum','999'),
        ('wrong_presence_support','split_presence.tsv',PRESENCE_FIELDS,'empirical_projected_ufboot_percent','3'),
        ('invented_presence_SH','split_presence.tsv',PRESENCE_FIELDS,'sh_alrt_percent','100'),
        ('wrong_conflict_support','conflicts.tsv.gz',CONFLICT_FIELDS,'empirical_projected_ufboot_percent_a','3'),
        ('invalid_quartet','conflicts.tsv.gz',CONFLICT_FIELDS,'witness_quartet','A;A;A;A')]
    rejected=[]
    for name,file,fields,key,value in mutations:
        target=root/name;shutil.copytree(valid,target);compressed=file.endswith('.gz');opener=gzip.open if compressed else open
        with opener(target/file,'rt') as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
        if key=='present':
            index=next(i for i,r in enumerate(rows) if r['present']=='True')
        elif file=='split_presence.tsv':index=next(i for i,r in enumerate(rows) if r['present']=='True')
        else:index=0
        rows[index][key]=value;rewrite(target/file,fields,rows,compressed)
        (target/'rehashed_artifacts.json').write_text(json.dumps({p.name:sha(p) for p in target.iterdir() if p.is_file()})+'\n')
        try:inspect(groups,universe,target)
        except (AssertionError,ValueError,KeyError) as e:rejected.append(dict(case=name,reason=str(e) or type(e).__name__))
        else:raise AssertionError('False comparison accepted: '+name)
    report=dict(status='passed_full_retained_baseline_comparison_software_cases',summary=summary,
        rehashed_false_exports_rejected=rejected,source_hashes={p:sha(p) for p in
        ['scripts/retained_tree_comparison_sources.py','scripts/compare_retained_tree_views.py','scripts/readback_retained_tree_comparisons.py',str(Path(__file__))]},
        scope='Complete4retainedsets/32view/64pair syntheticcomparison grid, rawDendroPy deletion/RF/bipartitioncompatibility/pathchecks andallpresence/quartet/supportthresholds verified. Fourteenrehashedfalseexports rejected. Treesandprojected supports synthetic; nooriginalbootstrap/closure/journal/productionproof orbiologicalpilot claimed.')
    (root/'receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],false_exports_rejected=len(rejected))),flush=True)


if __name__=='__main__':main()
