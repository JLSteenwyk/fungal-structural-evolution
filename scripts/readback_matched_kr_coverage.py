"""Summarize completed simulation audit and retain every optimization-review case."""
from collections import Counter
import json
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
from ancestral_chain_attempt import sha,write_json


def main():
    roots=[Path('results/model_validation/matched-kr-simulations-20260928-v1'),
           Path('results/model_validation/matched-kr-simulation-audit-20260928-v1'),
           Path('results/figures/matched-kr-coverage-20260928-v1')]
    receipts=[];bindings={};states={}
    for suffix in ['simulation-audit','coverage-report']:
        unit='fungal-matched-kr-'+suffix+'-20260928.service'
        state=dict(l.split('=',1) for l in subprocess.check_output(['systemctl','--user','show',unit,
            '-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
        assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0'),state
        states[unit]=state
    for root in roots:
        receipt=json.loads((root/'receipt.json').read_text());receipts.append(receipt)
        bindings[str(root/'receipt.json')]=sha(root/'receipt.json')
        for name,digest in receipt['artifacts'].items():
            assert sha(root/name)==digest;bindings[str(root/name)]=digest
    assert receipts[1]['source_receipt_sha256']==sha(roots[0]/'receipt.json')
    assert receipts[2]['source_audit_receipt_sha256']==sha(roots[1]/'receipt.json')
    assert receipts[0]['attempted']==receipts[1]['attempted']==15984
    manifest=json.loads((roots[0]/'manifest.json').read_text());reviews=[];causes=Counter()
    assert len(manifest)==15984
    assert dict(Counter(r['refit_status'] for r in manifest))==receipts[1]['status_counts']
    for row in manifest:
        if row['refit_status']=='refit_numerically_checked':continue
        assert sha(row['path'])==row['sha256']
        saved=json.loads(Path(row['path']).read_text());f=saved['refit']['fit']
        failed=[k for k,v in f['checks'].items() if (v if k=='upper_bound_contact' else not v)]
        causes.update(failed)
        full=[c['objective'] for c in f['candidates'] if all(c['active']) and c['start']!='original_unoptimized_reference']
        reviews.append(dict(case=row['case'],replicate=row['replicate'],path=row['path'],sha256=row['sha256'],
            failed_checks=';'.join(failed),full_face_objective_spread=float(np.ptp(full)),
            maximum_projected_gradient=float(np.max(abs(np.array(f['projected_gradient']))))))
    frame=pd.read_csv(roots[2]/'all_coefficient_coverage.tsv',sep='\t')
    assert len(frame)==112 and not frame.duplicated(['case','method','coefficient']).any()
    diagnostic={}
    for method,group in frame.groupby('method'):
        diagnostic[method]=dict(rows=len(group),marginal_upper_below_095=int((group.upper<.95).sum()),
            marginal_lower_above_095=int((group.lower>.95).sum()),
            unresolved_per_coefficient_range=[int(group.unresolved.min()),int(group.unresolved.max())])
    out=Path('results/model_validation/matched-kr-coverage-readback-20260928-v1');out.mkdir(parents=True,exist_ok=False)
    pd.DataFrame(reviews).to_csv(out/'optimization_review_cases.tsv',sep='\t',index=False)
    result=dict(status='completed_coverage_report_readback_and_review_inventory',attempted=15984,
        numerically_qualified_refits=receipts[1]['status_counts']['refit_numerically_checked'],
        review_refits=len(reviews),review_causes=dict(causes),method_descriptive_counts=diagnostic,
        full_face_objective_spread_quantiles={str(q):float(np.quantile([r['full_face_objective_spread'] for r in reviews],q)) for q in [0,.5,.9,1]},
        terminal_states=states,source_bindings=bindings,script_sha256=sha(__file__),
        artifacts={p.name:sha(p) for p in out.iterdir()},
        scope='All published coverage rows and184 source-bound review cases inventoried. Marginal '
              'interval exclusions are descriptive across56 dependent coefficient/design comparisons; '
              'not multiplicity-adjusted rejections. No global coverage certification or optimizer '
              'repair claimed. Full response/fit replay is bound through the completed audit.')
    write_json(out/'receipt.json',result)
    write_json(Path('metadata/matched_kr_coverage_completed_20260928.json'),result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_bindings','terminal_states','artifacts']},indent=2))


if __name__=='__main__':main()
