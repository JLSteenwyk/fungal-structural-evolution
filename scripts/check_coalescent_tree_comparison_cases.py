#!/usr/bin/env python3
"""Full 70-view/315-pair software grid with support boundaries and false exports."""
import argparse
import csv
import gzip
import json
from pathlib import Path

import dendropy

from compare_coalescent_species_trees import export
from readback_coalescent_tree_comparisons import inspect
from run_ortholog_pair_guide_comparison import sha


def read_table(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as f:
        reader=csv.DictReader(f,delimiter='\t');return reader.fieldnames,list(reader)


def write_table(path,fields,rows):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'wt') as f:
        writer=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)


def fixture(root):
    universe=set('ABCDEFGH');groups={}
    templates=['((A,B),(C,D),((E,F),(G,H)));','((A,C),(B,D),((E,G),(F,H)));']
    cohorts={'full_primary':universe,'exclude_sparse_boundary_taxon':universe-{'A'},
        'exclude_below10_both_alignments':universe-{'B','G'},'exclude_curated_hybrids':universe-{'C'},
        'exclude_hybrids_and_uncertain_labels':universe-{'A','E'}}
    for cohort,taxa in cohorts.items():
        outgroups=taxa&set('GH');views={};groups[cohort]=dict(taxa=taxa,outgroups=outgroups,
            roles=dict(ingroup=len(taxa)-len(outgroups),outgroup=len(outgroups)),views=views)
        for i in range(14):
            coal=i<6;kind='estimated' if coal else 'ml' if i%2==0 else 'consensus'
            name=('coal:' if coal else 'concat:')+str(i).zfill(2)
            raw=dendropy.Tree.get(data=templates[i%2],schema='newick',rooting='force-unrooted',preserve_underscores=True)
            for node in raw.preorder_node_iter():
                if node is not raw.seed_node:node.edge.length=.01
            prune=cohort!='full_primary' and not coal
            if coal:raw.retain_taxa_with_labels(sorted(taxa));raw.deroot()
            if coal:
                for j,node in enumerate(n for n in raw.preorder_node_iter() if n is not raw.seed_node and not n.is_leaf()):
                    value=[.949999,.95,.950001][j%3]
                    node.label='[pp1='+str(value)+';EN=25.0]'
            source=root/(cohort+'-'+name.replace(':','-')+'.tree')
            source.write_text(raw.as_string(schema='newick'))
            metric_tree=raw.clone(depth=2)
            if prune:metric_tree.retain_taxa_with_labels(sorted(taxa))
            metric_tree.deroot();metrics={}
            for j,node in enumerate(n for n in metric_tree.preorder_node_iter() if n is not metric_tree.seed_node and not n.is_leaf()):
                side={n.taxon.label for n in node.leaf_iter()};a,b=tuple(sorted(side)),tuple(sorted(taxa-side))
                key=a if (len(a),a)<=(len(b),b) else b
                pp=float(node.label[5:].split(';')[0]) if coal else None
                metrics[key]=dict(branch_length=node.edge.length,branch_length_unit='MAP_coalescent_units' if coal else
                    'original_path_sum_AA_substitutions_per_site' if prune else 'AA_substitutions_per_site',
                    local_posterior=pp,effective_genes=25.0 if coal else None,
                    empirical_ufb_percent=None if coal else [94.999,95.,95.001][j%3],
                    sh_alrt_percent=None if coal or prune or kind=='consensus' else [79.999,80.,80.001][j%3])
            assert len(metrics)==len(taxa)-3
            rule=('local_PP95' if coal else 'projected_empirical_UFB95_SH_unavailable' if prune else
                  'original_SH80_and_empirical_UFB95' if kind=='ml' else 'original_consensus_empirical_UFB95_SH_unavailable')
            views[name]=dict(family='coalescent' if coal else 'concatenated',kind=kind,source_tree=str(source),prune=prune,support_rule=rule,splits=metrics)
    return groups,universe


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path);args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False);source=args.output/'sources';source.mkdir();groups,universe=fixture(source)
    root=args.output/'full-grid';produced=export(groups,root);checked=inspect(groups,universe,root);assert produced==checked
    assert produced['cohorts']==5 and produced['tree_views']==70 and produced['comparison_rows']==315
    rejected=0
    mutations=[('comparisons.tsv','rf_distance'),('comparisons.tsv','normalized_rf'),('comparisons.tsv','both_support_rules_met_pairs'),
        ('comparisons.tsv','comparison_scope'),('split_presence.tsv.gz','present'),('split_presence.tsv.gz','support_rule_met'),
        ('split_presence.tsv.gz','branch_length_unit'),('split_presence.tsv.gz','local_posterior'),
        ('split_presence.tsv.gz','sh_alrt_percent'),('conflicts.tsv.gz','both_support_rules_met'),
        ('conflicts.tsv.gz','witness_quartet'),('role_boundaries.tsv','boundary_present'),('role_boundaries.tsv','outgroup')]
    for name,field in mutations:
        path=root/name;original=path.read_bytes();fields,rows=read_table(path)
        index=next((i for i,r in enumerate(rows) if r[field]!=''),0)
        value=rows[index][field]
        if value in ['True','False']:changed=str(value!='True')
        else:
            try:changed=str(float(value)+.125)
            except ValueError:changed='altered'
        rows[index][field]=changed;write_table(path,fields,rows)
        try:
            inspect(groups,universe,root)
        except (AssertionError,ValueError,KeyError,TypeError):rejected+=1
        else:raise AssertionError('False export accepted: '+field)
        finally:path.write_bytes(original)
    result=dict(status='passed_full_coalescent_reference_comparison_software_grid',**produced,
        altered_exports_rejected=rejected,source_hashes={p:sha(p) for p in ['scripts/coalescent_tree_comparison_sources.py',
            'scripts/compare_coalescent_species_trees.py','scripts/readback_coalescent_tree_comparisons.py',str(Path(__file__))]},
        scientific_eligibility=False,scope='Synthetic 5-cohort/70-view/315-pair topology grid and exact/below/above PP .95, UFB95 and SH80 support thresholds. Branch units and unavailable SH remain explicit. Independent raw-tree pruning/RF/bipartition witness readback rejects changed exports. No production source/journal completion or biological pilot claim.')
    (args.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if isinstance(v,(str,int,bool))}),flush=True)


if __name__=='__main__':main()
