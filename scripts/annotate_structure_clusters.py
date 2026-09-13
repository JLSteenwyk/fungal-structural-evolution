#!/usr/bin/env python3
"""Attach frozen taxon/marker provenance and pending validation flags to similarity groups."""
import argparse,json
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
from audit_busco_gene_copies import ROOT,sha,read_table
from assess_pae_sensitivity import checked_receipt
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--clusters',type=Path,required=True);p.add_argument('--edge-validation',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use new immutable annotations')
    checked_receipt(a.clusters);checked_receipt(a.edge_validation);config=json.loads((a.clusters/'config.json').read_text());ec=json.loads((a.edge_validation/'config.json').read_text())
    if ec['cluster_receipt_sha256']!=sha(a.clusters/'receipt.json'):raise ValueError('Edge review lineage differs')
    manifest=ROOT/'metadata/analysis_manifest.tsv';taxa={r['taxon_id']:r for r in read_table(manifest)};links=defaultdict(list)
    for source,path in [('AlphaFold',ROOT/'results/structural_markers/gdm-expanded-v1'),('ESMFold',ROOT/'results/structural_markers/esmfold-partial-v1')]:
        checked_receipt(path)
        if config['source_receipt_sha256'][source]!=sha(path/'receipt.json'):raise ValueError('Snapshot identity differs')
        for r in read_table(path/'marker_structure_links.tsv'):links[source,r['model_id']].append(r)
    edge=defaultdict(Counter)
    for r in read_table(a.edge_validation/'edge_review.tsv'):edge[r['representative_input_id']][r['status']]+=1
    groups=defaultdict(list);mapped=[]
    for r in read_table(a.clusters/'model_cluster_membership.tsv'):
        groups[r['representative_input_id']].append(r);matches=links[r['source'],r['model_id']]
        if not matches:raise ValueError('Missing model marker provenance')
        for m in matches:
            if m['sequence_sha256']!=r['sequence_sha256'] or m['model_sha256']!=r['model_sha256']:raise ValueError('Model identity mismatch')
            t=taxa[m['taxon_id']]
            mapped.append({'representative_input_id':r['representative_input_id'],'cluster_input_id':r['cluster_input_id'],'source':r['source'],'model_id':r['model_id'],'taxon_id':m['taxon_id'],'protein_id':m['protein_id'],'marker':m['marker'],'species_name':t['species_name'],'study_role':t['study_role'],'major_lineage':t['lineage'].split(';')[0]})
    annotations=defaultdict(list)
    for r in mapped:annotations[r['representative_input_id']].append(r)
    rows=[]
    for rep,models in sorted(groups.items()):
        items=annotations[rep];ts={r['taxon_id'] for r in items};markers={r['marker'] for r in items};counts=edge[rep]
        if sum(counts.values())!=len(models)-1:raise ValueError('Edge/model count mismatch')
        rows.append({'representative_input_id':rep,'models':len(models),'unique_sequences':len({r['sequence_sha256'] for r in models}),'taxa':len(ts),'ingroup_taxa':sum(taxa[t]['study_role']=='ingroup' for t in ts),'outgroup_taxa':sum(taxa[t]['study_role']=='outgroup' for t in ts),'sources':';'.join(sorted({r['source'] for r in models})),'markers':len(markers),'marker_ids_json':json.dumps(sorted(markers)),'major_lineages_json':json.dumps(sorted({r['major_lineage'] for r in items})),'minimum_length':min(int(r['length']) for r in models),'maximum_length':max(int(r['length']) for r in models),'model_weighted_median_mean_CA_plddt':float(np.median([float(r['mean_ca_plddt']) for r in models])),'models_mean_CA_plddt_below70':sum(float(r['mean_ca_plddt'])<70 for r in models),'edges_not_passing_both_approximate_directions':counts['one_direction_passes']+counts['neither_direction_passes'],'edge_qualification':'pending_exact_score_review' if len(models)>1 else 'singleton_no_edge_test','orthology_status':'not_established_by_structural_clustering'})
    a.output.mkdir(parents=True);write_table(a.output/'cluster_annotations.tsv',rows);write_table(a.output/'cluster_taxon_marker_links.tsv',mapped)
    r={'status':'complete_exploratory_cluster_provenance_annotation','cluster_receipt_sha256':sha(a.clusters/'receipt.json'),'edge_validation_receipt_sha256':sha(a.edge_validation/'receipt.json'),'manifest_sha256':sha(manifest),'script_sha256':sha(Path(__file__)),'clusters':len(rows),'models':sum(r['models'] for r in rows),'model_taxon_marker_links':len(mapped),'multi_marker_clusters':sum(r['markers']>1 for r in rows),'clusters_with_outgroup_taxa':sum(r['outgroup_taxa']>0 for r in rows),'clusters_with_failed_approximate_edges':sum(r['edges_not_passing_both_approximate_directions']>0 for r in rows),'interpretation':'Provenance annotation of original exploratory groups. Marker labels do not certify orthology; multiple markers can reflect shared domains, annotation differences or false grouping. Taxa and source composition describe uneven availability, not independent events or trait associations. Exact-score review, confidence/domain sensitivity and biological family inference remain pending.','artifacts':{p.name:sha(p) for p in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))

if __name__=='__main__':main()
