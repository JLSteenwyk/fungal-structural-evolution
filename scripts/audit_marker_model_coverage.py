#!/usr/bin/env python3
"""Audit model availability across every taxon, separating sequence and model gaps."""
import argparse
import csv
import json
from collections import Counter,defaultdict
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import ROOT,sha
from prepare_paired_phylogenetic_inputs import write_table


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalog',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--inventory-inputs',type=Path,required=True)
    args=parser.parse_args()
    receipt=checked_receipt(args.catalog)
    inventory=checked_receipt(args.inventory_inputs)
    reuse_rows=list(csv.DictReader((args.inventory_inputs/'all_marker_links.tsv').open(),delimiter='\t'))
    reuse={(r['marker'],r['taxon_id']):r for r in reuse_rows}
    if len(reuse)!=len(reuse_rows):raise ValueError('Duplicate inventory marker links')
    manifest_path=ROOT/'metadata/analysis_manifest.tsv'
    manifest=list(csv.DictReader(manifest_path.open(),delimiter='\t'))
    if len({r['taxon_id'] for r in manifest})!=len(manifest):raise ValueError('Duplicate manifest taxa')
    mapping_path=ROOT/'results/phylogeny/markers-full-v1/protein_mapping.tsv'
    if sha(mapping_path)!=receipt['marker_mapping_sha256']:raise ValueError('Changed marker mapping')
    mapping=list(csv.DictReader(mapping_path.open(),delimiter='\t'))
    by_key={(r['marker'],r['taxon_id']):r for r in mapping}
    if len(by_key)!=len(mapping):raise ValueError('Duplicate marker/taxon records')
    if set(reuse)!=set(by_key):raise ValueError('Inventory and marker universes differ')
    for key,row in by_key.items():
        if any(reuse[key][k]!=row[k] for k in ['protein_id','sequence_sha256']):raise ValueError('Inventory sequence identity differs')
    nmarkers=len({r['marker'] for r in mapping})
    known={r['taxon_id'] for r in manifest}
    if any(r['taxon_id'] not in known for r in mapping):raise ValueError('Unknown marker taxon')
    links=list(csv.DictReader((args.catalog/'marker_model_links.tsv').open(),delimiter='\t'))
    if len({(r['marker'],r['taxon_id']) for r in links})!=len(links):raise ValueError('Duplicate model links')
    for r in links:
        original=by_key[r['marker'],r['taxon_id']]
        if any(r[k]!=original[k] for k in ['protein_id','sequence_sha256']):raise ValueError('Model/marker identity differs')
    seq_count=Counter(r['taxon_id'] for r in mapping);model_count=Counter(r['taxon_id'] for r in links)
    linked={(r['marker'],r['taxon_id']) for r in links}
    nominated=Counter();no_candidate=Counter()
    for key,row in reuse.items():
        if key in linked:
            if row['reuse_state']!='verified_reuse_receipt':raise ValueError('Catalog model absent from later inventory evidence')
        elif row['reuse_state'] in ['verified_reuse_receipt','nominated_reuse_pending']:nominated[row['taxon_id']]+=1
        elif row['reuse_state']=='no_candidate_in_completed_inventory':no_candidate[row['taxon_id']]+=1
        else:raise ValueError('Complete inventory required')
    taxa=[];groups=defaultdict(list)
    for r in manifest:
        taxon=r['taxon_id']; sequence=seq_count[taxon]; models=model_count[taxon]
        if not 0<=models<=sequence<=nmarkers:raise ValueError('Invalid marker counts')
        row={'taxon_id':taxon,'species_name':r['species_name'],'study_role':r['study_role'],
            'manifest_lineage_group':r['lineage'].split(';')[0], 'total_markers':nmarkers,
            'recovered_sequence_markers':sequence,'markers_with_source_model':models,
            'reuse_candidate_outside_catalog':nominated[taxon],'no_candidate_in_query_snapshot':no_candidate[taxon],
            'sequence_present_model_missing':sequence-models,'sequence_missing':nmarkers-sequence,
            'model_fraction_of_recovered_markers':models/sequence if sequence else '',
            'model_fraction_of_all_markers':models/nmarkers}
        if models+nominated[taxon]+no_candidate[taxon]!=sequence:raise ValueError('Coverage categories do not partition sequence markers')
        taxa.append(row);groups[r['study_role'],row['manifest_lineage_group']].append(row)
    summaries=[]
    for (role,group),rows in sorted(groups.items()):
        total=nmarkers*len(rows);seq=sum(r['recovered_sequence_markers'] for r in rows);models=sum(r['markers_with_source_model'] for r in rows)
        summaries.append({'study_role':role,'manifest_lineage_group':group,'taxa':len(rows),
            'taxa_with_any_model':sum(r['markers_with_source_model']>0 for r in rows),
            'taxa_without_models':sum(r['markers_with_source_model']==0 for r in rows),
            'marker_slots':total,'sequence_markers':seq,'model_linked_markers':models,
            'model_fraction':models/total,
            'reuse_candidate_outside_catalog_fraction':sum(r['reuse_candidate_outside_catalog'] for r in rows)/total,
            'no_candidate_in_query_snapshot_fraction':sum(r['no_candidate_in_query_snapshot'] for r in rows)/total,
            'missing_sequence_fraction':(total-seq)/total})
    if args.output.exists():raise FileExistsError('Use a new immutable coverage output')
    args.output.mkdir(parents=True)
    write_table(args.output/'taxon_coverage.tsv',taxa);write_table(args.output/'lineage_coverage.tsv',summaries)
    fig,ax=plt.subplots(figsize=(11,10))
    y=np.arange(len(summaries));left=np.zeros(len(y))
    for key,color,label in [('model_fraction','#26817e','Downloaded GDM model in catalog'),('reuse_candidate_outside_catalog_fraction','#e4be74','Reuse candidate, outside this catalog'),('no_candidate_in_query_snapshot_fraction','#bdcccf','No candidate in completed query snapshot'),('missing_sequence_fraction','#eeeeee','Marker sequence unavailable')]:
        values=np.array([r[key] for r in summaries]);ax.barh(y,values,left=left,color=color,label=label,height=.75);left+=values
    ax.set_yticks(y,[r['manifest_lineage_group']+(' [outgroup]' if r['study_role']=='outgroup' else '')+f" (n={r['taxa']})" for r in summaries])
    ax.invert_yaxis();ax.set_xlim(0,1);ax.set_xlabel(f'Fraction of {nmarkers} marker slots per taxon')
    ax.set_title('Sequence recovery and reusable structural-model availability\nFull 526-taxon design; frozen GDM pipeline catalog',fontsize=12)
    ax.legend(loc='upper center',bbox_to_anchor=(.45,-.065),ncol=1,frameon=False,fontsize=9)
    ax.spines[['top','right']].set_visible(False)
    fig.text(.53,.015,'Model availability does not establish confidence-qualified sites or evolutionary independence.\nLineage bins follow the sampling manifest and are not all the same taxonomic rank.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.12,1,1));fig.savefig(args.output/'marker_model_coverage.svg',metadata={'Date':None});fig.savefig(args.output/'marker_model_coverage.png',dpi=160)
    result={'status':'complete_full_sampling_model_coverage_audit','taxa':len(taxa),'markers':nmarkers,
        'taxa_with_models':sum(r['markers_with_source_model']>0 for r in taxa),'taxa_without_models':sum(r['markers_with_source_model']==0 for r in taxa),
        'sequence_marker_links':len(mapping),'model_marker_links':len(links),'lineage_groups':len(summaries),
        'inventory_input_receipt_sha256':sha(args.inventory_inputs/'receipt.json'),
        'catalog_receipt_sha256':sha(args.catalog/'receipt.json'),'sampling_manifest_sha256':sha(manifest_path),'script_sha256':sha(Path(__file__)),
        'interpretation':'Every taxon retained. Distinguishes downloaded source models, reuse candidates outside this catalog, no candidate in the completed query snapshot, and marker non-recovery. Source catalog and query inventory are separately frozen acquisition checkpoints; outside-catalog candidates can include later downloads or other pipelines. None of these missing states proves biological absence. Model availability precedes residue/confidence validation. Group bars pool equal 125-marker denominators per taxon; linked models can be reused across identical sequences and are not independent evolutionary observations.',
        'artifacts':{p.name:sha(p) for p in args.output.iterdir()}}
    (args.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
