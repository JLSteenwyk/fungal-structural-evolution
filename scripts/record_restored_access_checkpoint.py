#!/usr/bin/env python3
"""Verify newly terminal analysis evidence after execution access restoration."""
import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess
from ancestral_chain_attempt import sha,write_json


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();assert not args.output.exists()
    checked={}
    def verify(p,h):
        p=Path(p)
        if str(p) not in checked:checked[str(p)]=sha(p)
        assert checked[str(p)]==h,str(p)
    units=['fungal-full-matched-working-models-20260927','fungal-full-matched-working-model-output-audit-20260927','fungal-full-working-model-analytic-stationarity-20260927','fungal-full-analytic-stationarity-readback-20260927','fungal-full-model-grid-export-20260927','fungal-full-fastml-variants-20260928','fungal-full-fastml-roundoff-readback-20260928']
    states={}
    for unit in units:
        s=dict(l.split('=',1) for l in subprocess.check_output(['systemctl','--user','show',unit+'.service','-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
        assert s==dict(ActiveState='inactive',Result='success',ExecMainStatus='0');states[unit]=s
    root=Path('results/structural_comparisons/full-matched-working-models-20260927-v1')
    audit=Path('metadata/full_matched_working_model_output_audit_20260927.json');a=json.loads(audit.read_text());assert a['tree_fit_dispositions']==144040;verify(root/'receipt.json',a['source_receipt_sha256'])
    analytic=Path('results/model_validation/full-working-model-analytic-stationarity-20260927-v1')
    readback=Path('results/model_validation/full-working-model-analytic-stationarity-readback-20260927-v1')
    export=Path('results/structural_comparisons/full-working-model-grid-export-20260927-v1')
    roots=[root,analytic,readback,export];reports={}
    for folder in roots:
        p=folder/'receipt.json';r=json.loads(p.read_text());reports[str(folder)]=r
        for name,h in r.get('artifacts',{}).items():verify(folder/name,h)
        for name,h in r.get('source_bindings',{}).items():verify(name,h)
    ar=reports[str(analytic)];rr=reports[str(readback)];er=reports[str(export)]
    verify(analytic/'receipt.json',rr['source_receipt_sha256']);verify(root/'receipt.json',ar['production_receipt_sha256']);verify(audit,ar['production_audit_sha256'])
    assert ar['dispositions']==rr['dispositions']==er['unique_fits']==144040
    assert rr['remaining_review_cases']==er['remaining_review_unique_fits']==658
    fp=Path('results/ancestral/full-fastml-roundoff-readback-20260928-v1/receipt.json');fr=json.loads(fp.read_text());assert len(fr['readbacks'])==312 and not fr['unresolved']
    rows=[]
    for binding in fr['readbacks'].values():
        p=Path(binding['path']);verify(p,binding['sha256']);r=json.loads(p.read_text());rows.append(r)
        for name,h in r.get('evidence',{}).items():
            verify(name,h)
            ep=Path(name)
            if ep.name=='receipt.json':
                ev=json.loads(ep.read_text())
                for artifact,expected in ev.get('artifacts',{}).items():verify(ep.parent/artifact,expected)
    statuses=dict(Counter(r['status'] for r in rows));assert statuses=={'integrity_and_independent_replay_complete_not_fit_qualification':306,'no_coded_characters':6}
    summaries={}
    for variant in ['precision-only','precision-cache-refresh']:
        selected=[r for r in rows if r.get('variant')==variant and 'likelihood_difference' in r];assert len(selected)==153
        summaries[variant]=dict(fits=153,probability_rows=sum(r['probability_rows'] for r in selected),maximum_absolute_likelihood_difference=max(abs(r['likelihood_difference']) for r in selected),likelihood_difference_exceeds_1e_minus5=sum(abs(r['likelihood_difference'])>1e-5 for r in selected),maximum_probability_difference=max(r['maximum_probability_difference'] for r in selected),raw_boundary_excursions=sum(r['raw_boundary_excursions'] for r in selected),maximum_raw_boundary_error=max(r['maximum_raw_boundary_error'] for r in selected))
    benchmark=Path('metadata/ancestral_categorical_runtime_benchmark_20260928.json');b=json.loads(benchmark.read_text());assert b['patterns']==120
    for name,h in b['pins'].items():verify(name,h)
    roots_evidence=[audit,fp,benchmark,*[p/'receipt.json' for p in roots]]
    write_json(args.output,dict(status='new_terminal_stages_verified_after_access_restoration',services=states,source_roots={str(p):sha(p) for p in roots_evidence},checked_artifact_files=len(checked),matched_models=dict(fits=144040,original_status_counts=a['status_counts'],analytic_gradient_transitions=rr['gradient_threshold_transitions'],remaining_review_fits=658,review_reasons=rr['review_reason_counts'],expanded_setting_rows=er['full_setting_fits']),fastml=dict(status_counts=statuses,variants=summaries),benchmark=dict(patterns=120,elapsed_seconds=b['elapsed_seconds']),script_sha256=sha(__file__),scope='Terminal states and bound artifacts checked; no independent refitting or derivative recomputation. Optimization review, FastML numerical/model qualification, uncertainty and all eight scientific aims remain open.'))
    print(json.dumps(dict(checked_artifact_files=len(checked),matched_review_fits=658,fastml=summaries)))


if __name__=='__main__':main()
