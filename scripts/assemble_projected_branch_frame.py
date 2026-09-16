#!/usr/bin/env python3
"""Join audited marker branch estimates to provisional guide mappings without path allocation."""
import argparse
import json
from collections import defaultdict,Counter
from pathlib import Path
from Bio import Phylo
from audit_joint_path_uncertainty import checked,rows,sha
from run_paired_marker_fits import tree_edges
from prepare_paired_phylogenetic_inputs import write_table

LABELS=['aa','3di_af','3di_af_empirical','3di_llm']


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['projection','projection-readback','output']:ap.add_argument('--'+name,type=Path,required=True)
    ap.add_argument('--cohorts',nargs='+',required=True)
    ap.add_argument('--rates',nargs='+',type=Path,required=True)
    ap.add_argument('--audits',nargs='+',type=Path,required=True)
    a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    if not len(a.cohorts)==len(a.rates)==len(a.audits) or len(set(a.cohorts))!=len(a.cohorts):raise ValueError('Distinct cohorts and paired sources required')
    pr=checked(a.projection);pa=json.loads(a.projection_readback.read_text())
    if pa['status']!='passed_full_pruning_projection_readback' or pa['projection_receipt_sha256']!=sha(a.projection/'receipt.json'):raise ValueError('Projection audit mismatch')
    pins={str(a.projection/'receipt.json'):sha(a.projection/'receipt.json'),str(a.projection_readback):sha(a.projection_readback)}
    mappings=defaultdict(dict)
    for row in rows(a.projection/'edge_projection.tsv'):
        key=(row['cohort'],row['marker'],row['marker_split_taxa'])
        if row['guide'] in mappings[key]:raise ValueError('Duplicate guide projection')
        mappings[key][row['guide']]=json.loads(row['full_edge_ids_json'])
    sensitivities=list(rows(a.projection/'guide_sensitivity.tsv'))
    if {x['cohort'] for x in sensitivities}!=set(a.cohorts):raise ValueError('Projection cohort universe mismatch')
    output=[];sources=[]
    for cohort,rate,audit in zip(a.cohorts,a.rates,a.audits):
        rr=json.loads((rate/'receipt.json').read_text());ar=checked(audit)
        if ar['status']!='passed_full_site_rate_output_audit' or ar['rate_receipt_sha256']!=sha(rate/'receipt.json'):raise ValueError('Rate audit lineage mismatch')
        pins[str(rate/'receipt.json')]=sha(rate/'receipt.json');pins[str(audit/'receipt.json')]=sha(audit/'receipt.json')
        fit_summary={(x['marker'],x['fit']):x for x in rows(audit/'fit_summary.tsv')}
        grouped=defaultdict(list)
        for row in sensitivities:
            if row['cohort']==cohort:grouped[row['marker']].append(row)
        if len(grouped)!=ar['markers']:raise ValueError('Marker universe mismatch')
        for marker,items in sorted(grouped.items()):
            estimates={};sites={};warnings={}
            expected={tuple(x['marker_split_taxa'].split(',')) for x in items}
            taxa=None
            for label in LABELS:
                receipt_path=rate/marker/(label+'.receipt.json');tree=rate/marker/(label+'.treefile')
                if sha(receipt_path)!=rr['fit_receipts'][marker+'/'+label]:raise ValueError('Fit receipt mismatch')
                fr=json.loads(receipt_path.read_text())
                if sha(tree)!=fr['artifacts'][tree.name]:raise ValueError('Changed fitted tree')
                tips={t.name for t in Phylo.read(tree,'newick').get_terminals()}
                if taxa is not None and tips!=taxa:raise ValueError('Model tip universe differs')
                taxa=tips;estimates[label]=tree_edges(tree,taxa)
                if set(estimates[label])!=expected:raise ValueError('Fitted topology differs from projection')
                fs=fit_summary[marker,label]
                if fs['rate_heterogeneity']!='G4':raise ValueError('Gamma baseline required')
                sites[label]=int(fs['sites']);warnings[label]=int(fs['warnings'])
                pins[str(tree)]=sha(tree);pins[str(receipt_path)]=sha(receipt_path)
            if len(set(sites.values()))!=1:raise ValueError('Unpaired alignment site counts')
            for item in items:
                split=tuple(item['marker_split_taxa'].split(','));key=(cohort,marker,item['marker_split_taxa']);mapping=mappings[key]
                if len(mapping)!=2:raise ValueError('Both guide dispositions required')
                unique=item['status']=='unique_same_full_split_in_both_guides'
                if unique and (any(len(v)!=1 for v in mapping.values()) or len({v[0] for v in mapping.values()})!=1):raise ValueError('Unique mapping mismatch')
                row=dict(cohort=cohort,marker=marker,marker_taxa=len(taxa),aligned_sites=sites['aa'],marker_split_taxa=item['marker_split_taxa'],marker_terminal=item['marker_terminal'],projection_status=item['status'],unique_full_edge_id=next(iter(mapping.values()))[0] if unique else '',guide_edge_sets_json=json.dumps(mapping,sort_keys=True,separators=(',',':')))
                for label in LABELS:
                    row[label+'_branch_length']=estimates[label][split]
                    row[label+'_fit_warnings']=warnings[label]
                output.append(row)
        sources.append(dict(cohort=cohort,rates=str(rate),audit=str(audit),markers=len(grouped)))
    if len(output)!=pr['marker_edges']:raise ValueError('Incomplete edge frame')
    coverage=defaultdict(list)
    for row in output:
        if row['unique_full_edge_id']:coverage[row['cohort'],row['unique_full_edge_id']].append(row)
    cov=[dict(cohort=k[0],full_edge_id=k[1],markers=len(v),marker_terminal=v[0]['marker_terminal'],marker_ids_json=json.dumps(sorted(x['marker'] for x in v))) for k,v in sorted(coverage.items())]
    a.output.mkdir(parents=True)
    write_table(a.output/'branch_frame.tsv',output);write_table(a.output/'unique_edge_coverage.tsv',cov)
    result=dict(status='complete_projected_gamma_branch_frame',rows=len(output),numerical_branch_estimates=4*len(output),sources=sources,projection_status_counts=dict(Counter(x['projection_status'] for x in output)),source_pins=pins,script_sha256=sha(Path(__file__)),interpretation='Paired expected substitutions per aligned amino-acid or 3Di state under conditional G4 models. Not physical displacement, time rates or an acceleration/selection test. Revised ESM Gamma optimization is used. Unsupported homogeneous guide correspondence and all collapsed/discordant mappings are retained; no branch length is divided among full-guide edges. Unique-edge coverage counts are eligibility, not replication or support. FreeRate, supported guide, phylogenetic dependence and family/coverage sensitivity remain required.',artifacts={p.name:sha(p) for p in a.output.iterdir()})
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['source_pins','artifacts']},indent=2))

if __name__=='__main__':main()
