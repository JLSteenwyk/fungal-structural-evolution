#!/usr/bin/env python3
"""Quantify actual paired phylogenetic coverage across every sampled lineage."""
import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT, sha, read_table
from prepare_paired_phylogenetic_inputs import write_table

REASONS=['observed','noncanonical_or_missing_sequence','no_structural_mapping',
         'invalid_native_feature','low_feature_plddt','high_feature_pae']


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable coverage output')
    receipt=checked_receipt(a.inputs)
    if receipt['status']!='complete_paired_phylogenetic_input_preparation':raise ValueError('Complete paired inputs required')
    manifest=ROOT/'metadata/analysis_manifest.tsv'
    taxa={r['taxon_id']:r for r in read_table(manifest)}
    matrix=ROOT/receipt['source_receipts']['matrix']['path']
    if sha(matrix/'receipt.json')!=receipt['source_receipts']['matrix']['sha256']:raise ValueError('Changed matrix receipt')
    if json.loads((matrix/'receipt.json').read_text())['manifest_sha256']!=sha(manifest):raise ValueError('Changed taxon manifest')
    markers={r['marker']:r for r in read_table(a.inputs/'marker_summary.tsv')}
    rows=read_table(a.inputs/'taxon_coverage.tsv')
    seen=set();counts=defaultdict(Counter);per_marker=Counter();model_taxa=set()
    for r in rows:
        key=(r['taxon_id'],r['marker'])
        if key in seen or key[0] not in taxa or key[1] not in markers:raise ValueError('Repeated or unknown coverage identity')
        seen.add(key);n=int(r['original_marker_columns']);c=counts[key[0]]
        if sum(int(r[k]) for k in REASONS)!=n:raise ValueError('Site partition does not sum')
        if n!=int(markers[key[1]]['original_marker_columns']):raise ValueError('Marker site denominator differs')
        cutoff=max(50,math.ceil(.3*n))
        eligible=int(r['observed'])>=cutoff
        if int(r['required_observed_columns'])!=cutoff or r['taxon_eligible']!=str(eligible):raise ValueError('Eligibility definition differs')
        c['marker_slots']+=1;c['matrix_columns']+=n
        for name in REASONS:c[name]+=int(r[name])
        c['markers_with_model']+=bool(r['model_name'])
        c['eligible_markers']+=eligible
        c['usable_markers']+=eligible and markers[key[1]]['status']=='ready_for_inference'
        per_marker[key[1]]+=eligible
        if r['model_name']:model_taxa.add(key[0])
    if seen!={(t,m) for t in taxa for m in markers}:raise ValueError('Incomplete full-design coverage grid')
    for marker,r in markers.items():
        if per_marker[marker]!=int(r['eligible_taxa']) or (per_marker[marker]>=4)!=(r['status']=='ready_for_inference'):raise ValueError('Marker readiness differs')
    taxon_rows=[];groups=defaultdict(list)
    for taxon,t in taxa.items():
        c=counts[taxon]
        row={'taxon_id':taxon,'species_name':t['species_name'],'study_role':t['study_role'],
             'lineage_group':t['lineage'].split(';')[0],**{k:c[k] for k in ['marker_slots','markers_with_model','eligible_markers','usable_markers','matrix_columns',*REASONS]},
             'observed_fraction_of_matrix':c['observed']/c['matrix_columns']}
        taxon_rows.append(row);groups[(row['study_role'],row['lineage_group'])].append(row)
    group_rows=[]
    for (role,lineage),members in sorted(groups.items()):
        group_rows.append({'study_role':role,'lineage_group':lineage,'taxa':len(members),
            'taxa_with_any_model':sum(r['markers_with_model']>0 for r in members),
            'taxa_with_any_usable_marker':sum(r['usable_markers']>0 for r in members),
            'taxa_with_at_least10_usable_markers':sum(r['usable_markers']>=10 for r in members),
            'taxa_with_at_least50_usable_markers':sum(r['usable_markers']>=50 for r in members),
            **{k:sum(r[k] for r in members) for k in ['marker_slots','eligible_markers','usable_markers','matrix_columns',*REASONS]}})
    a.output.mkdir(parents=True)
    write_table(a.output/'taxon_coverage.tsv',taxon_rows);write_table(a.output/'lineage_coverage.tsv',group_rows)
    result={'status':'complete_paired_lineage_coverage','source_receipt_sha256':sha(a.inputs/'receipt.json'),
            'manifest_sha256':sha(manifest),'script_sha256':sha(Path(__file__)),
            'taxa':len(taxa),'markers':len(markers),'taxon_marker_rows':len(seen),
            'taxa_with_any_model':len(model_taxa),'taxa_with_any_usable_marker':sum(r['usable_markers']>0 for r in taxon_rows),
            'represented_groups':len(groups),'groups_with_usable_markers':sum(r['taxa_with_any_usable_marker']>0 for r in group_rows),
            'site_partition':{name:sum(r[name] for r in taxon_rows) for name in REASONS},
            'interpretation':'Every manifest entry retained. Usable means eligible in a marker with >=4 eligible taxa, not resolved orthology, supported branches or selection suitability. Site reasons follow sequential masking and are not independent causal effects. Missing sequence combines matrix gaps/noncanonical symbols; missing structure is not biological absence. Lineage labels are manifest groups, not inferred clades.',
            'artifacts':{x.name:sha(x) for x in a.output.iterdir() if x.is_file()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
