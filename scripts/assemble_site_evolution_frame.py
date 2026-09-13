#!/usr/bin/env python3
"""Assemble matched site-level estimates and covariates without inferential testing."""
import argparse,csv,gzip,json
from collections import defaultdict
from pathlib import Path
import numpy as np
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha,read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['inputs','topologies','topology_audit','diagnostic','rates','rate_audit','projection','output']:p.add_argument('--'+name.replace('_','-'),type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use immutable output')
    names=['inputs','topologies','topology_audit','diagnostic','rate_audit','projection']
    r={k:checked_receipt(getattr(a,k)) for k in names};rate=json.loads((a.rates/'receipt.json').read_text());rateconfig=json.loads((a.rates/'config.json').read_text());topconfig=json.loads((a.topologies/'config.json').read_text())
    if r['topology_audit']['assessment_receipt_sha256']!=sha(a.topologies/'receipt.json') or r['rate_audit']['rate_receipt_sha256']!=sha(a.rates/'receipt.json') or rate['config_sha256']!=sha(a.rates/'config.json') or r['topologies']['config_sha256']!=sha(a.topologies/'config.json'):
        raise ValueError('Audit lineage differs')
    if topconfig['source_receipts']['diagnostic']!=sha(a.diagnostic/'receipt.json') or any(v!=sha(a.inputs/'receipt.json') for v in [topconfig['source_receipts']['inputs'],rateconfig['source_receipts']['inputs'],r['diagnostic']['source_receipts']['inputs'],r['projection']['source_receipts']['inputs']]) or topconfig['source_receipts']['fits']!=rateconfig['source_receipts']['fits']:
        raise ValueError('Sources refer to different paired data or fits')
    rows=read_table(a.topologies/'site_topology_sensitivity.tsv');indexed={(x['marker'],x['paired_column_1based']):x for x in rows}
    if len(indexed)!=len(rows):raise ValueError('Duplicate base sites')
    markers={x['marker']:x for x in read_table(a.inputs/'marker_summary.tsv') if x['status']=='ready_for_inference'}
    confidence=defaultdict(list);taxa=defaultdict(set)
    with gzip.open(a.projection/'paired_site_accessibility.tsv.gz','rt') as f:
        for x in csv.DictReader(f,delimiter='\t'):
            key=x['marker'],x['paired_column_1based']
            if key not in indexed or x['taxon_id'] in taxa[key]:raise ValueError('Unexpected or duplicate observed site')
            taxa[key].add(x['taxon_id']);confidence[key].append(float(x['ca_plddt']))
    fitstats={(x['marker'],x['fit']):x for x in read_table(a.rate_audit/'fit_summary.tsv')};labels=['aa','3di_af','3di_af_empirical','3di_llm'];seen=set()
    for x in read_table(a.rate_audit/'site_rates.tsv'):
        key=x['marker'],x['paired_column_1based'];label=x['fit']
        if key not in indexed or label not in labels or (key,label) in seen:raise ValueError('Unexpected rate cell')
        seen.add((key,label));row=indexed[key];fit=fitstats[x['marker'],label]
        row[label+'_posterior_mean_relative_rate']=x['posterior_mean_relative_rate'];row[label+'_modal_rate_category']=x['modal_rate_category'];row[label+'_modal_category_relative_rate']=x['modal_category_relative_rate'];row[label+'_gamma_alpha']=fit['gamma_alpha'];row[label+'_fit_warning_count']=fit['warnings']
    if len(seen)!=len(rows)*4 or len(seen)!=r['rate_audit']['site_rate_rows']:raise ValueError('Rate grid incomplete')
    summaries=[]
    for key,row in indexed.items():
        values=confidence.pop(key);observed=len(taxa.pop(key));total=int(markers[key[0]]['eligible_taxa'])
        if observed!=int(row['observed_taxa']) or observed!=len(values):raise ValueError('Observed site count differs')
        row['tree_taxa']=total;row['observed_taxon_fraction']=observed/total;row['median_focal_ca_plddt']=float(np.median(values));row['minimum_focal_ca_plddt']=min(values);row['focal_ca_plddt_q25']=float(np.quantile(values,.25));row['focal_ca_plddt_q75']=float(np.quantile(values,.75))
    if confidence or taxa or set(markers)!={x['marker'] for x in rows}:raise ValueError('Incomplete base universe')
    for marker in sorted(markers):
        selected=[x for x in rows if x['marker']==marker];summaries.append({'marker':marker,'sites':len(selected),'taxon_site_observations':sum(int(x['observed_taxa']) for x in selected),'tree_taxa':int(markers[marker]['eligible_taxa']),'minimum_observed_taxon_fraction':min(x['observed_taxon_fraction'] for x in selected),'median_observed_taxon_fraction':float(np.median([x['observed_taxon_fraction'] for x in selected]))})
    a.output.mkdir(parents=True);write_table(a.output/'site_evolution_frame.tsv',rows);write_table(a.output/'marker_summary.tsv',summaries)
    out={'status':'complete_matched_site_evolution_frame','markers':len(markers),'sites':len(rows),'rate_cells':len(seen),'taxon_site_observations':sum(x['taxon_site_observations'] for x in summaries),'source_receipts':{k:sha(getattr(a,k)/'receipt.json') for k in names+['rates']},'script_sha256':sha(Path(__file__)),'interpretation':'One row per original paired site, preserving parsimony/exposure/topology diagnostics and all four conditional empirical-Bayes site rates. Coverage and focal confidence are observed covariates; gamma parameters and fit warnings remain explicit. This assembly performs no regression, significance test, causal adjustment, selection test or ancestral reconstruction. Sites and estimated quantities remain dependent.','artifacts':{f.name:sha(f) for f in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
