#!/usr/bin/env python3
"""Summarize audited alignment sensitivity, using cases as descriptive units."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    root = Path('results/cds/codon-alignment-divergence-comparison-20260927-v1')
    proof_path = Path('results/cds/codon-divergence-readback-handoff-20260927-v1/readback.json')
    proof = json.loads(proof_path.read_text())
    assert proof['status'] == 'passed_full_normalized_codon_divergence_independent_readback'
    assert proof['source_receipt_sha256'] == sha(root/'receipt.json')
    for name, digest in json.loads((root/'receipt.json').read_text())['artifacts'].items(): assert sha(root/name) == digest
    cases = pd.read_csv(root/'cases.tsv',sep='\t',keep_default_na=False)
    pairs = pd.read_csv(root/'pairs.tsv',sep='\t')
    out = cases[cases.comparison_status == 'matched'].copy().set_index('case_id')
    assert len(out) == proof['dispositions']['matched'] == 1625
    assert set(pairs.case_id) == set(out.index)
    out['historical_flags_present'] = out.historical_review_flags.ne('')
    topology_root = Path('results/cds/codon-alignment-tree-comparison-20260927-v1')
    tp = Path('metadata/codon_tree_comparison_completed_readback_20260927.json')
    assert json.loads(tp.read_text())['source_receipt_sha256'] == sha(topology_root/'receipt.json')
    tr = json.loads((topology_root/'receipt.json').read_text())
    assert tr['artifacts']['cases.tsv'] == sha(topology_root/'cases.tsv')
    topology = pd.read_csv(topology_root/'cases.tsv',sep='\t',keep_default_na=False).set_index('case_id')
    assert set(topology.index[topology.comparison_status=='matched']) == set(out.index)
    out['topology_changed'] = topology.loc[out.index,'rf_distance'].astype(int).gt(0)
    summaries = []
    for metric in ('global_omega','tree_ds_equal_alternative','tree_dn_equal_alternative'):
        a = out['original_'+metric].astype(float)
        b = out['local_'+metric].astype(float)
        assert np.isfinite(a).all() and np.isfinite(b).all() and (a>=0).all() and (b>=0).all()
        good = (a>0)&(b>0)
        out[metric+'_log2_local_over_original'] = np.nan
        out.loc[good,metric+'_log2_local_over_original'] = np.log2(b[good])-np.log2(a[good])
        out[metric+'_ratio_disposition'] = np.where(good,'positive_both',np.where((a==0)&(b==0),'zero_both',np.where(a==0,'zero_original','zero_local')))
    for metric in ('ds_equal_alternative','dn_equal_alternative'):
        a,b=pairs['original_'+metric],pairs['local_'+metric]
        valid=(a>0)&(b>0)
        values=np.log2(b[valid])-np.log2(a[valid])
        frame=pd.DataFrame({'case_id':pairs.loc[valid,'case_id'],'absolute_log2_fold':values.abs()})
        out['pair_'+metric+'_positive_both']=frame.groupby('case_id').size().reindex(out.index,fill_value=0)
        out['pair_'+metric+'_not_positive_both']=pairs.groupby('case_id').size().reindex(out.index)-out['pair_'+metric+'_positive_both']
        out['pair_'+metric+'_median_absolute_log2_fold']=frame.groupby('case_id').absolute_log2_fold.median().reindex(out.index)
    metrics=[c for c in out if c.endswith('_log2_local_over_original') or c.endswith('_median_absolute_log2_fold')]
    for group,subset in [('all',out),('topology_changed',out[out.topology_changed]),('topology_unchanged',out[~out.topology_changed])]:
        for metric in metrics:
            values=subset[metric].dropna()
            r=dict(group=group,metric=metric,cases=len(subset),defined_cases=len(values),undefined_cases=len(subset)-len(values))
            for q,value in values.quantile([0,.25,.5,.75,1]).items():r['quantile_'+str(q)]=value
            summaries.append(r)
    target=Path('results/cds/codon-divergence-sensitivity-summary-20260927-v1')
    target.mkdir(parents=True,exist_ok=False)
    out.to_csv(target/'matched_cases.tsv',sep='\t')
    pd.DataFrame(summaries).to_csv(target/'case_quantiles.tsv',sep='\t',index=False)
    result=dict(status='complete_codon_divergence_sensitivity_summary_pending_readback',matched_cases=len(out),comparison_receipt_sha256=sha(root/'receipt.json'),comparison_readback_sha256=sha(proof_path),topology_receipt_sha256=sha(topology_root/'receipt.json'),script_sha256=sha(__file__),ratio_dispositions={m:out[m+'_ratio_disposition'].value_counts().to_dict() for m in ('global_omega','tree_ds_equal_alternative','tree_dn_equal_alternative')},artifacts={p.name:sha(p) for p in target.iterdir()},scope='Descriptive case-level sensitivity only. Signed log2(local/original) requires both values positive; zeros remain explicit and receive no pseudocount. Pair summaries are first aggregated within case. Cases remain phylogenetically/family correlated; no tests or confidence intervals. Changed alignment coverage and topology both contribute. Historical warnings persist; neither alignment is established as correct.')
    (target/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__=='__main__':main()
