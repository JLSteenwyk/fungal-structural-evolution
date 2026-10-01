#!/usr/bin/env python3
"""Synthetic path-collapse/frequency/corruption checks, not a biological pilot."""
import argparse
from collections import Counter
import csv
import io
import json
from pathlib import Path
import shutil
from Bio import Phylo
from project_pmsf_retained_taxa import raw_edges,project
from readback_pmsf_retained_taxa import inspect_tables
from retained_taxon_projection_sources import EDGE_FIELDS,BOOT_FIELDS,BOUNDARY_FIELDS
from run_ortholog_pair_guide_comparison import sha


def export(path,fields,rows):
    with Path(path).open('w') as handle:
        writer=csv.DictWriter(handle,fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    root=args.output;root.mkdir(parents=True,exist_ok=False);run=root/'synthetic-run';run.mkdir();valid=root/'valid';valid.mkdir()
    a='(((A:0.1,B:0.1):0.2,X:0.1):0.3,(C:0.1,D:0.1):0.4,(E:0.1,Y:0.1):0.5);'
    b='(((A:0.1,X:0.1):0.2,B:0.1):0.3,(C:0.1,E:0.1):0.4,(D:0.1,Y:0.1):0.5);'
    (run/'pmsf.ufboot').write_text((a+'\n')*700+(b+'\n')*300)
    (run/'pmsf.treefile').write_text(a+'\n');(run/'pmsf.contree').write_text(b+'\n')
    universe=set('ABXCDEY');policies={}
    for name,tips in [('remove_X','ABCDEY'),('remove_X_Y','ABCDE'),('remove_X_A','BCDEY'),('remove_X_B_Y','ACDE')]:
        retained=set(tips);outgroups=retained & set('EY')
        policies[name]=dict(taxa=retained,outgroups=outgroups,roles=dict(ingroup=len(retained-outgroups),outgroup=len(outgroups)))
    sources={'synthetic':dict(spec=dict(run=str(run)))}
    counts={p:Counter() for p in policies}
    for text,n in [(a,700),(b,300)]:
        edges=raw_edges(Phylo.read(io.StringIO(text),'newick'),universe)
        for policy,spec in policies.items():
            for key in project(edges,spec['taxa']):counts[policy][key]+=n
    assert counts['remove_X_Y'][('A','B')]==1000,'Merged original AB/ABX splits must count once per replicate'
    assert counts['remove_X_Y'][('C','D')]==700 and counts['remove_X_Y'][('C','E')]==300
    edges=[];boots=[];boundary=[]
    for policy,spec in policies.items():
        retained=spec['taxa'];side=spec['outgroups'];key=min(tuple(sorted(side)),tuple(sorted(retained-side)),key=lambda s:(len(s),s))
        for split,n in sorted(counts[policy].items()):
            boots.append(dict(baseline_run='synthetic',policy=policy,split_taxa_json=json.dumps(split),
                replicates_with_projected_split=n,empirical_projected_ufboot_percent=n/10))
        for kind,text in [('ml',a),('consensus',b)]:
            p=project(raw_edges(Phylo.read(io.StringIO(text),'newick'),universe),retained)
            for split,(length,n) in sorted(p.items()):
                edges.append(dict(baseline_run='synthetic',policy=policy,tree_type=kind,split_taxa_json=json.dumps(split),
                    projected_branch_length_sum=length,original_edge_components=n,
                    empirical_projected_ufboot_percent=counts[policy][split]/10,sh_alrt_percent=None))
            boundary.append(dict(baseline_run='synthetic',policy=policy,tree_type=kind,retained_taxa=len(retained),
                retained_ingroup=spec['roles']['ingroup'],retained_outgroup=spec['roles']['outgroup'],boundary_taxa_json=json.dumps(key),
                boundary_present=key in p,empirical_projected_ufboot_percent=counts[policy][key]/10,sh_alrt_percent=None))
    export(valid/'projected_tree_edges.tsv',EDGE_FIELDS,edges);export(valid/'projected_bootstrap_splits.tsv',BOOT_FIELDS,boots)
    export(valid/'role_boundary.tsv',BOUNDARY_FIELDS,boundary)
    summary=inspect_tables(sources,universe,policies,valid);assert summary['tree_views']==8 and summary['projected_bootstrap_tree_states']==4000
    mutations=[('inflated_merged_support','projected_bootstrap_splits.tsv','replicates_with_projected_split','1700'),
               ('wrong_bootstrap_percentage','projected_bootstrap_splits.tsv','empirical_projected_ufboot_percent','99.0'),
               ('wrong_path_length','projected_tree_edges.tsv','projected_branch_length_sum','0.987'),
               ('wrong_edge_component_count','projected_tree_edges.tsv','original_edge_components','100'),
               ('inherited_SH_aLRT','projected_tree_edges.tsv','sh_alrt_percent','100'),
               ('wrong_retained_role','role_boundary.tsv','retained_outgroup','99')]
    rejected=[]
    for name,file,field,value in mutations:
        target=root/name;shutil.copytree(valid,target)
        with (target/file).open() as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
        rows[0][field]=value;export(target/file,list(rows[0]),rows)
        # Rehash every false export. The verifier must reject its contents.
        (target/'rehashed_artifacts.json').write_text(json.dumps({p.name:sha(p) for p in target.glob('*.tsv')})+'\n')
        try:inspect_tables(sources,universe,policies,target)
        except (AssertionError,ValueError,KeyError) as e:rejected.append(dict(case=name,reason=str(e) or type(e).__name__))
        else:raise AssertionError('False projection accepted: '+name)
    proof=dict(status='passed_retained_taxon_projection_path_and_frequency_software_cases',valid_summary=summary,
        recomputed_merged_AB_frequency_percent=100,original_AB_split_frequency_percent=70,
        rehashed_false_exports_rejected=rejected,source_hashes={p:sha(p) for p in
        ['scripts/project_pmsf_retained_taxa.py','scripts/readback_pmsf_retained_taxa.py','scripts/retained_taxon_projection_sources.py',str(Path(__file__))]},
        scope='One synthetic7-tip pair of trees,1000synthetic raw replicates andfour retainedtipsets/eightviews. Raw AB support70% becomes100% after droppingX, countedonceper replicate despitecollapsed originalsplits. IndependentDendroPy pruning checksall projectededges/lengths/pathcounts/frequencies/roles; sixrehashed falseexports rejected. No inference/priorjournal/productionproof orbiologicalpilot claimed.')
    (root/'receipt.json').write_text(json.dumps(proof,indent=2)+'\n')
    print(json.dumps(dict(status=proof['status'],false_exports_rejected=len(rejected))),flush=True)


if __name__=='__main__':main()
