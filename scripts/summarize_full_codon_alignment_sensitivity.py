#!/usr/bin/env python3
"""Summarize verified correspondence sensitivity, retaining all cases and pairs."""
import csv,json,hashlib,statistics
from pathlib import Path
from collections import defaultdict

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    root=Path('results/cds/full-group-codon-realignment-projection-20260927-v1');rp=root/'receipt.json';proof=Path('metadata/full_codon_realignment_projection_readback_20260927.json');r=json.loads(rp.read_text());pr=json.loads(proof.read_text());assert pr['status']=='passed_full_codon_projection_and_correspondence_readback' and pr['producer_receipt_sha256']==sha(rp)
    for name in ['cases.tsv','pair_correspondence.tsv']:assert sha(root/name)==r['artifacts'][name]
    with (root/'cases.tsv').open() as f:cases=list(csv.DictReader(f,delimiter='\t'))
    groups=defaultdict(list)
    with (root/'pair_correspondence.tsv').open() as f:
        for row in csv.DictReader(f,delimiter='\t'):groups[row['case_id']].append(row)
    output=[];totals=defaultdict(int)
    for case in cases:
        cid=case['case_id'];pairs=groups[cid];assert len(pairs)==int(case['taxa'])*(int(case['taxa'])-1)//2
        sums={k:sum(int(x[k]) for x in pairs) for k in ['original_pairs','original_preserved_retained','original_absent_from_raw','original_preserved_but_filtered','local_retained_not_original']}
        for k,v in sums.items():totals[k]+=v
        n=sums['original_pairs'];assert n==sum(sums[k] for k in ['original_preserved_retained','original_absent_from_raw','original_preserved_but_filtered'])
        fractions=[float(x['original_fraction_retained']) for x in pairs if x['original_fraction_retained']!='']
        output.append(dict(case_id=cid,taxon_pairs=len(pairs),original_columns=case['original_columns'],local_retained_columns=case['retained_codon_columns'],coverage_disposition=case['coverage_disposition'],median_pair_retention=statistics.median(fractions) if fractions else '',minimum_pair_retention=min(fractions) if fractions else '',**sums,original_pair_weighted_retention=sums['original_preserved_retained']/n if n else ''))
    out=Path('results/cds/full-group-alignment-sensitivity-summary-20260927-v1');out.mkdir(parents=True,exist_ok=False)
    with (out/'cases.tsv').open('w') as f:
        w=csv.DictWriter(f,list(output[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(output)
    valid=[x['median_pair_retention'] for x in output if x['median_pair_retention']!=''];den=totals['original_pairs']
    result=dict(status='complete_verified_codon_alignment_sensitivity_summary',source_receipt_sha256=sha(rp),source_readback_sha256=sha(proof),script_sha256=sha(__file__),cases=len(output),taxon_pairs=sum(x['taxon_pairs'] for x in output),original_columns=sum(int(x['original_columns']) for x in cases),local_retained_columns=sum(int(x['retained_codon_columns']) for x in cases),median_of_case_median_pair_retention=statistics.median(valid),minimum_case_median_pair_retention=min(valid),cases_with_any_original_pair_absent_from_raw=sum(x['original_absent_from_raw']>0 for x in output),cases_with_any_original_pair_lost_only_to_filter=sum(x['original_preserved_but_filtered']>0 for x in output),pair_observation_totals=dict(totals),original_pair_observation_fractions={k:totals[k]/den for k in ['original_preserved_retained','original_absent_from_raw','original_preserved_but_filtered']},artifacts={'cases.tsv':sha(out/'cases.tsv')},scope='Descriptive all-case summary. Pair observations repeat residues/taxa/families and are not independent samples; fractions are explicitly pair-observation weighted, while median-of-case-medians is case weighted. Local group alignment and masking differ from global alignment; neither extra columns nor retained correspondence establishes correct homology, statistical power or selection eligibility.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
