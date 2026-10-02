#!/usr/bin/env python3
"""Exercise all88 views/560 comparisons with unresolved native consensus trees."""
import argparse
import copy
import csv
import gzip
import json
from pathlib import Path

import dendropy

from check_coalescent_tree_comparison_cases import fixture
from compare_native_subset_tree_views import export
from readback_native_subset_tree_views import inspect
from run_ortholog_pair_guide_comparison import sha


def tables(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as f:
        reader=csv.DictReader(f,delimiter='\t');return reader.fieldnames,list(reader)


def write_table(path,fields,rows):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'wt') as f:
        writer=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    out=args.output;out.mkdir(parents=True,exist_ok=False);src=out/'sources';src.mkdir();groups,universe=fixture(src);groups.pop('full_primary')
    coal=['coal:'+a+':'+s for a in ['profile','mafft'] for s in ['uncontracted','sh_alrt10','sh_alrt80']]
    refs=['concat:'+r+':'+k for r in ['profile_profile','profile_mafft','mafft_profile','mafft_mafft'] for k in ['ml','consensus']]
    native=['native:'+r+':'+k for r in ['profile_profile','profile_mafft','mafft_profile','mafft_mafft'] for k in ['ml','consensus']]
    templates=['((A,B),(C,D),((E,F),(G,H)));','((A,C),(B,D),((E,G),(F,H)));']
    for cohort,g in groups.items():
        old=g['views'];g['views']={new:old[('coal:' if i<6 else 'concat:')+str(i).zfill(2)] for i,new in enumerate(coal+refs)}
        taxa=g['taxa']
        for i,name in enumerate(native):
            kind='ml' if i%2==0 else 'consensus'
            text='('+','.join(sorted(taxa))+');' if i in [3,7] else templates[(i//2)%2]
            tree=dendropy.Tree.get(data=text,schema='newick',rooting='force-unrooted',preserve_underscores=True)
            tree.retain_taxa_with_labels(sorted(taxa));tree.deroot();metrics={}
            for node in tree.preorder_node_iter():
                if node is not tree.seed_node:node.edge.length=.02
            for j,node in enumerate(n for n in tree.preorder_node_iter() if n is not tree.seed_node and not n.is_leaf()):
                side={n.taxon.label for n in node.leaf_iter()};a,b=tuple(sorted(side)),tuple(sorted(taxa-side));key=min(a,b,key=lambda s:(len(s),s))
                sh=[79.999,80.,80.001][j%3] if kind=='ml' else None;boot=[94.999,95.,95.001][j%3]
                node.label=str(sh)+'/'+str(boot) if kind=='ml' else str(boot)
                metrics[key]=dict(branch_length=node.edge.length,branch_length_unit='AA_substitutions_per_site',local_posterior=None,
                    effective_genes=None,empirical_ufb_percent=boot,sh_alrt_percent=sh)
            path=src/(cohort+'-'+name.replace(':','-')+'.tree');path.write_text(tree.as_string(schema='newick'))
            g['views'][name]=dict(family='native_pmsf',kind=kind,source_tree=str(path),prune=False,splits=metrics,
                support_rule='native_SH80_and_empirical_UFB95' if kind=='ml' else 'native_consensus_empirical_UFB95_SH_unavailable')
    root=out/'full-grid';produced=export(groups,root);checked=inspect(groups,universe,root);assert produced==checked
    assert produced['cohorts']==4 and produced['tree_views']==88 and produced['comparison_rows']==560
    _,rows=tables(root/'comparisons.tsv');empty=[r for r in rows if r['normalized_rf_by_observed_splits']=='']
    assert len(empty)==4 and all(r['rf_distance']=='0' and r['normalized_rf']=='0.0' for r in empty)
    rejected=0
    changes=[('comparisons.tsv','splits_a'),('comparisons.tsv','possible_internal_splits'),('comparisons.tsv','rf_distance'),
        ('comparisons.tsv','normalized_rf'),('comparisons.tsv','normalized_rf_by_observed_splits'),('comparisons.tsv','comparison_scope'),
        ('split_presence.tsv.gz','present'),('split_presence.tsv.gz','support_rule_met'),('split_presence.tsv.gz','sh_alrt_percent'),
        ('split_presence.tsv.gz','branch_length_unit'),('conflicts.tsv.gz','witness_quartet'),('conflicts.tsv.gz','both_support_rules_met'),
        ('role_boundaries.tsv','outgroup'),('role_boundaries.tsv','boundary_present')]
    for name,field in changes:
        path=root/name;original=path.read_bytes();fields,rows=tables(path);row=next((r for r in rows if r[field]!=''),rows[0]);value=row[field]
        if value in ['True','False']:row[field]=str(value!='True')
        else:
            try:row[field]=str(float(value)+.125)
            except ValueError:row[field]='altered'
        write_table(path,fields,rows)
        try:inspect(groups,universe,root)
        except (AssertionError,ValueError,KeyError,TypeError):rejected+=1
        else:raise AssertionError('Altered native subset export accepted '+field)
        finally:path.write_bytes(original)
    # A fabricated value in the genuine zero-observed-splits denominator must reject.
    path=root/'comparisons.tsv';original=path.read_bytes();fields,rows=tables(path)
    next(r for r in rows if r['normalized_rf_by_observed_splits']=='')['normalized_rf_by_observed_splits']='0'
    write_table(path,fields,rows)
    try:inspect(groups,universe,root)
    except (AssertionError,ValueError,KeyError,TypeError):rejected+=1
    else:raise AssertionError('Unavailable star/star denominator changed to zero')
    finally:path.write_bytes(original)
    paths=['scripts/native_subset_tree_comparison_sources.py','scripts/compare_native_subset_tree_views.py',
           'scripts/readback_native_subset_tree_views.py','scripts/check_coalescent_tree_comparison_cases.py',str(Path(__file__))]
    result=dict(status='passed_full_native_subset_88_view_560_pair_software_contracts',**produced,
        zero_observed_split_denominators=4,altered_exports_rejected=rejected,source_hashes={p:sha(p) for p in paths},
        scientific_eligibility=False,scope='Full four-cohort/88-view/560-pair synthetic software grid:32native ML/consensus views,24coalescent candidates,32projected references. Includes unresolved/all-star native consensus states, exact/below/above support cutoffs, separate units and unavailable0/0observed-split normalization. Complete independent raw-tree pruning/RF/incompatibility/witness/boundary readback; changed exports reject. No actual subset inference, production closure or biological pilot claim.')
    with (out/'receipt.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if isinstance(v,(str,int,bool))}),flush=True)


if __name__=='__main__':main()
