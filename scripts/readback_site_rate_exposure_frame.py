#!/usr/bin/env python3
"""Independently verify every joined rate/exposure field and alignment covariate."""
import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha


def table(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def indexed(rows, columns):
    result={tuple(r[k] for k in columns):r for r in rows}
    if len(result)!=len(rows):raise ValueError('Duplicate table keys')
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['frame','rates','rate-readback','exposure','exposure-readback','inputs','output']:
        ap.add_argument('--'+name,type=Path,required=True)
    a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    receipts={k:checked_receipt(getattr(a,k)) for k in ['frame','rates','exposure','inputs']}
    frame=receipts['frame']
    if frame['status']!='complete_site_rate_exposure_analysis_frame':raise ValueError('Incomplete frame')
    for k in ['rates','exposure','inputs']:
        if frame['source_receipts'][k]!=sha(getattr(a,k)/'receipt.json'):raise ValueError('Frame lineage differs')
    rr=json.loads(a.rate_readback.read_text());er=json.loads(a.exposure_readback.read_text())
    if (rr['status']!='passed_full_rate_comparison_readback' or rr['comparison_receipt_sha256']!=sha(a.rates/'receipt.json') or frame['rate_readback_sha256']!=sha(a.rate_readback)):
        raise ValueError('Rate readback mismatch')
    if er['status']!='passed_full_site_parsimony_exposure_readback' or er['source_receipt_sha256']!=sha(a.exposure/'receipt.json'):raise ValueError('Exposure readback mismatch')
    for k in ['rates','exposure']:
        if receipts[k]['source_receipts']['inputs']!=sha(a.inputs/'receipt.json'):raise ValueError('Input cohort differs')
    keys=['marker','paired_column_1based'];exposures=indexed(table(a.exposure/'site_parsimony_exposure.tsv'),keys)
    rate=indexed(table(a.rates/'sites.tsv'),keys+['fit'])
    fits=['aa','3di_af','3di_af_empirical','3di_llm'];alphabet='ACDEFGHIKLMNPQRSTVWY'
    if set(rate)!={(*k,f) for k in exposures for f in fits}:raise ValueError('Rate/exposure universe differs')
    ready={r['marker']:r for r in table(a.inputs/'marker_summary.tsv') if r['status']=='ready_for_inference'}
    if len(ready)!=receipts['inputs']['ready_markers']:raise ValueError('Marker count differs')
    output=indexed(table(a.frame/'site_rate_exposure.tsv'),keys)
    if set(output)!=set(exposures):raise ValueError('Frame universe differs')
    checked=set();observations=0;rate_values=0
    for marker,m in ready.items():
        sequences=list(SeqIO.parse(a.inputs/marker/'aa.faa','fasta'));n=len(sequences)
        if n!=int(m['eligible_taxa']) or len({x.id for x in sequences})!=n:raise ValueError('Taxon grid differs')
        columns=table(a.inputs/marker/'columns.tsv');width=len(columns)
        if width!=int(m['retained_columns']) or any(len(s.seq)!=width for s in sequences):raise ValueError('Width differs')
        for j,col in enumerate(columns):
            key=marker,str(j+1);row=output[key];source=exposures[key];checked.add(key)
            if source['matrix_column_1based']!=col['matrix_column_1based']:raise ValueError('Column mapping differs')
            expected_fields=set(source)|{'marker_taxa','marker_columns','observed_fraction','aa_entropy_nats'}|{'aa_count_'+c for c in alphabet}|{f+'_'+v for f in fits for v in ['gamma_rate','freerate_rate']}
            if set(row)!=expected_fields:raise ValueError('Frame schema differs')
            for field,value in source.items():
                if row[field]!=value:raise ValueError(('Exposure value differs',key,field))
            chars=[str(s.seq[j]) for s in sequences if s.seq[j]!='?'];total=len(chars)
            if not total or set(chars)-set(alphabet):raise ValueError('Invalid observed state')
            counts={c:chars.count(c) for c in alphabet}
            if int(row['marker_taxa'])!=n or int(row['marker_columns'])!=width or int(row['observed_taxa'])!=total:raise ValueError('Counts differ')
            for c,v in counts.items():
                if int(row['aa_count_'+c])!=v:raise ValueError('Composition differs')
            entropy=math.log(total)-math.fsum(v*math.log(v) for v in counts.values() if v)/total
            if not math.isclose(float(row['aa_entropy_nats']),entropy,rel_tol=1e-12,abs_tol=1e-12) or not math.isclose(float(row['observed_fraction']),total/n,rel_tol=1e-12,abs_tol=1e-12):raise ValueError('Derived covariate differs')
            for fit in fits:
                for field in ['gamma_rate','freerate_rate']:
                    if row[fit+'_'+field]!=rate[(*key,fit)][field]:raise ValueError('Rate value differs')
                    rate_values+=1
            observations+=total
    if checked!=set(output) or (len(ready),len(checked),observations,rate_values)!=(frame['markers'],frame['sites'],frame['taxon_site_observations'],frame['site_rate_values']):raise ValueError('Incomplete frame')
    if indexed(table(a.frame/'fit_diagnostics.tsv'),['marker','fit'])!=indexed(table(a.rates/'fit_summary.tsv'),['marker','fit']):raise ValueError('Fit diagnostics differ')
    result=dict(status='passed_full_site_rate_exposure_frame_readback',markers=len(ready),sites=len(checked),taxon_site_observations=observations,rate_values=rate_values,frame_receipt_sha256=sha(a.frame/'receipt.json'),rate_readback_sha256=sha(a.rate_readback),exposure_readback_sha256=sha(a.exposure_readback),script_sha256=sha(Path(__file__)),scope='Every copied exposure and rate value, diagnostics row, exact site/fit grid, amino-acid composition, coverage and entropy independently read back. No new rate estimation or biological inference.')
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
