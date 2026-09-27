#!/usr/bin/env python3
"""Export matched ordinary-ML gains with all fit and nesting-review flags retained."""
import argparse,csv,json,math,subprocess,time
from collections import Counter
from pathlib import Path
import pandas as pd
import psutil
from screen_duplication_domain_alignment_coverage import sha

LABELS=['linear','quadratic','cubic']
CONTRASTS=[('linear','quadratic'),('quadratic','cubic'),('linear','cubic')]
PASS='ml_candidate_passed_numerical_optimization_checks'


def compare(payloads):
    row={};valid={}
    for label,p in zip(LABELS,payloads):
        row[label+'_status']=p['status']
        error=p['status']=='fit_error_requires_review'
        row[label+'_negative_log_likelihood']='' if error else p['negative_profiled_ml']
        row[label+'_fixed_coefficients']='' if error else len(p['beta'])
        row[label+'_records']='' if error else p['records']
        if not error:
            assert p['likelihood']=='ordinary_gaussian_ml'
            assert math.isfinite(p['negative_profiled_ml'])
            valid[label]=p
    assert len({p['records'] for p in valid.values()})<=1
    violations=0
    for reduced,full in CONTRASTS:
        name=reduced+'_to_'+full
        if reduced not in valid or full not in valid:
            row[name+'_log_likelihood_gain']='';row[name+'_nesting_status']='missing_fit'
            continue
        a,b=valid[reduced],valid[full]
        assert a['active_covariates']==b['active_covariates'][:len(a['active_covariates'])]
        gain=a['negative_profiled_ml']-b['negative_profiled_ml']
        tolerance=1e-7+1e-9*max(abs(a['negative_profiled_ml']),abs(b['negative_profiled_ml']))
        bad=gain < -tolerance
        row[name+'_log_likelihood_gain']=gain
        row[name+'_nesting_status']='worse_full_model_requires_review' if bad else 'nondecreasing_within_numeric_tolerance'
        violations+=bad
    row['comparison_status']=('missing_fit_requires_review' if len(valid)<3 else
        'nested_likelihood_violation_requires_review' if violations else
        'optimization_review_required' if any(p['status']!=PASS for p in payloads) else
        'candidate_passed_checks_not_inferential_acceptance')
    return row


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);args=ap.parse_args()
    plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        assert sha(args.plan)==ph
        for path,digest in plan['pins'].items():assert sha(path)==digest,path
    verify();dep=plan['producer']
    while True:
        try:
            p=psutil.Process(dep['pid'])
            if p.create_time()!=dep['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==dep['cmdline']
        except psutil.NoSuchProcess:break
        print('waiting_for_full_polynomial_ml_audit',dep['pid'],flush=True);time.sleep(30)
    terminal=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',dep['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert terminal==dict(ActiveState='inactive',Result='success',ExecMainStatus='0'),terminal
    verify();root=Path(plan['fits']);receipt=json.loads((root/'receipt.json').read_text());rh=sha(root/'receipt.json')
    audit=json.loads(Path(plan['audit']).read_text());ah=sha(plan['audit'])
    assert audit['status']=='passed_full_polynomial_ml_output_integrity_audit' and audit['source_receipt_sha256']==rh
    assert receipt['plan_sha256']==sha(plan['fit_plan']) and receipt['tree_fit_dispositions']==432120
    for name,h in receipt['artifacts'].items():assert sha(root/name)==h
    manifest={}
    for line in (root/'fit_manifest.jsonl').open():
        row=json.loads(line);key=row['fit_input_id'],row['tree'];assert key not in manifest;manifest[key]=row
    assert len(manifest)==432120
    links=Path(plan['links']);completion=json.loads(Path(plan['links_completion']).read_text())
    assert completion['status']=='complete_verified_identical_observation_polynomial_input_links'
    for name,h in completion['artifacts'].items():assert sha(links/name)==h
    triples=pd.read_csv(links/'unique_comparison_sets.tsv',sep='\t');settings=pd.read_csv(links/'full_setting_links.tsv',sep='\t')
    assert len(triples)==28808 and len(settings)==82944
    trees=sorted({tree for _,tree in manifest});assert len(trees)==5
    output=[];used=set()
    for ix,ids in enumerate(triples.to_dict('records'),1):
        for tree in trees:
            payloads=[]
            for degree,label in enumerate(LABELS,1):
                key=ids[label+'_input_id'],tree;assert key not in used;used.add(key)
                entry=manifest[key];path=Path(entry['path']);assert sha(path)==entry['sha256']
                saved=json.loads(path.read_text())
                assert (saved['fit_input_id'],saved['tree'],saved['polynomial_degree'])==(*key,degree)
                payloads.append(saved['payload'])
            output.append(dict(**ids,tree=tree,**compare(payloads)))
        if ix%1000==0:print('Compared ML triplets',ix,'/ 28808',flush=True)
    assert used==set(manifest)
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    unique=pd.DataFrame(output);unique.to_csv(out/'unique_comparisons.tsv',sep='\t',index=False)
    # Read every emitted scalar back, including explicit blank failed-fit metrics.
    with (out/'unique_comparisons.tsv').open() as f:
        exported=list(csv.DictReader(f,delimiter='\t'))
    assert len(exported)==len(output)==144040
    for actual,expected in zip(exported,output):
        assert set(actual)==set(expected)
        for field,value in expected.items():
            if isinstance(value,(float,int)):assert math.isfinite(float(actual[field])) and float(actual[field])==value
            else:assert actual[field]==value
    keys=[label+'_input_id' for label in LABELS]
    full=settings.merge(unique,on=keys,validate='many_to_many',sort=False)
    assert len(full)==414720
    for label in LABELS:
        present=full[label+'_records'].ne('')
        assert (full.loc[present,label+'_records'].astype(int)==full.loc[present,'records']).all()
    # Each original setting must appear once per tree; all linked identifiers retained.
    identity=[c for c in settings.columns if c!='records']
    assert not full.duplicated(identity+['tree']).any()
    assert full.groupby(identity,dropna=False).size().eq(5).all()
    full.to_csv(out/'full_setting_comparisons.tsv',sep='\t',index=False)
    pd.testing.assert_frame_equal(pd.read_csv(out/'full_setting_comparisons.tsv',sep='\t',float_precision='round_trip').fillna(''),full.fillna(''),check_dtype=False)
    verify();assert sha(root/'receipt.json')==rh and sha(plan['audit'])==ah
    result=dict(status='complete_polynomial_ml_comparison_export_with_full_serialized_readback',plan_sha256=ph,prerequisite_terminal_state=terminal,unique_triplet_tree_rows=len(unique),full_setting_tree_rows=len(full),fit_dispositions_consumed=len(used),comparison_status_counts=unique.comparison_status.value_counts().to_dict(),source_receipt_sha256=rh,source_audit_sha256=ah,artifacts={p.name:sha(p) for p in out.iterdir()},scope='Same-observation ordinary ML gains for all three nested degree contrasts, with raw negative gains retained and numerical nesting violations flagged at 1e-7 + 1e-9 max absolute objective. Failed and review fits remain visible. Every exported scalar and full setting map read back. No chi-square reference, p-values, significance, predictive preference, global-optimum guarantee or calibrated uncertainty.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
