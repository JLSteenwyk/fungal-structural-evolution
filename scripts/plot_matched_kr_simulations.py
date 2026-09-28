"""Report all synthetic coverage outcomes only after complete replay audit."""
import argparse
import json
from pathlib import Path
import subprocess
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psutil
from ancestral_chain_attempt import sha,write_json
from matched_calibration_intervals import exact_coverage_bounds


def coverage_table(plan, kr, baseline):
    expected={c['id'] for c in plan['cases']}
    assert set(kr)==set(baseline)==expected
    rows=[]
    for case in plan['cases']:
        for method, summary, key in [('KR candidate',kr[case['id']],'marginal_coverage_intervals'),
                                      ('Conditional t',baseline[case['id']],'marginal_monte_carlo_intervals')]:
            values=summary[key];assert len(values)==case['coefficients']
            for coefficient,value in enumerate(values):
                assert value['attempted']==plan['replicates']
                recomputed=exact_coverage_bounds(value['covered'],value['unresolved'],value['attempted'],value['alpha'])
                assert value==recomputed
                rows.append(dict(case=case['id'],records=case['records'],method=method,
                    coefficient=coefficient,background_ratio=case['ratios'][0],family_ratio=case['ratios'][1],
                    species_ratio=case['ratios'][2],**value))
    frame=pd.DataFrame(rows)
    assert not frame.duplicated(['case','method','coefficient']).any()
    assert len(frame)==2*sum(c['coefficients'] for c in plan['cases'])==112
    return frame


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--plan',type=Path,required=True)
    args=parser.parse_args();config=json.loads(args.plan.read_text());ch=sha(args.plan)
    def verify():
        assert sha(args.plan)==ch
        for path,digest in config['pins'].items():assert sha(path)==digest,path
    verify();launch=json.loads(Path(config['audit_launch']).read_text())
    while True:
        state=dict(l.split('=',1) for l in subprocess.check_output(['systemctl','--user','show',launch['unit'],
            '-p','ActiveState','-p','MainPID','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
        if state['ActiveState'] in ['inactive','failed']:break
        assert state['ActiveState']=='active' and int(state['MainPID'])==launch['pid']
        try:
            p=psutil.Process(launch['pid']);assert p.create_time()==launch['created'] and p.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:time.sleep(1);continue
        time.sleep(30)
    assert state['ActiveState']=='inactive' and state['Result']=='success' and state['ExecMainStatus']=='0',state
    audit_plan=json.loads(Path(config['audit_plan']).read_text())
    assert launch['plan_sha256']==sha(config['audit_plan'])
    root=Path(audit_plan['output']);receipt=json.loads((root/'receipt.json').read_text())
    assert receipt['status']=='completed_all_synthetic_response_fit_interval_replays' and receipt['attempted']==15984
    assert receipt['plan_sha256']==sha(config['audit_plan'])
    for name,digest in receipt['artifacts'].items():assert sha(root/name)==digest
    plan=json.loads(Path(audit_plan['producer_plan']).read_text())
    source=Path(plan['output'])/'receipt.json'
    assert receipt['source_receipt_sha256']==sha(source)
    kr=json.loads((root/'kr_coverage.json').read_text());baseline=json.loads((root/'conditional_t_coverage.json').read_text())
    frame=coverage_table(plan,kr,baseline)
    out=Path(config['output']);out.mkdir(parents=True,exist_ok=False)
    frame.to_csv(out/'all_coefficient_coverage.tsv',sep='\t',index=False)
    pd.testing.assert_frame_equal(pd.read_csv(out/'all_coefficient_coverage.tsv',sep='\t'),frame,check_exact=False,rtol=1e-12,atol=1e-14)
    focal=frame[frame.coefficient==0]
    fig,ax=plt.subplots(figsize=(10,9))
    labels=[]
    for y,case in enumerate(plan['cases']):
        labels.append(f"n={case['records']}; ratios="+','.join(f'{v:g}' for v in case['ratios']))
        for method,color,offset in [('KR candidate','#167d9a',-.15),('Conditional t','#a74c39',.15)]:
            row=focal[(focal.case==case['id'])&(focal.method==method)].iloc[0]
            ax.plot([row.lower,row.upper],[y+offset]*2,color=color,lw=1.3,
                    label=method if y==0 else None)
            ax.plot([row.observed_fraction_lower,row.observed_fraction_upper],[y+offset]*2,color=color,lw=5)
            ax.plot(row.observed_fraction_lower,y+offset,'o',color=color,ms=3)
    ax.axvline(.95,color='black',ls='--',lw=1,label='Nominal 95%')
    ax.set_yticks(range(len(labels)),labels);ax.invert_yaxis();ax.set_xlim(0,1)
    ax.set_xlabel('Coverage of generating intercept')
    ax.set_title('Synthetic Gaussian designs: interval coverage after variance refitting\n999 attempted responses per configuration')
    ax.legend(loc='lower left');ax.spines[['top','right']].set_visible(False)
    fig.text(.5,.02,'Thin segments: marginal 95% Monte Carlo envelopes. Thick segments: unresolved-outcome range.\nAll attempts retained; ratios ordered background, family, species. No simultaneous or real-grid coverage claim.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.07,1,1));fig.savefig(out/'intercept_coverage.png',dpi=160);fig.savefig(out/'intercept_coverage.pdf');plt.close(fig)
    verify()
    write_json(out/'receipt.json',dict(status='complete_audited_synthetic_coverage_table_and_figure',
        plan_sha256=ch,source_audit_receipt_sha256=sha(root/'receipt.json'),table_rows=len(frame),
        figure_intercept_rows=len(focal),attempted_per_case=plan['replicates'],
        artifacts={p.name:sha(p) for p in out.iterdir()},
        scope='All coefficients and both methods tabulated; figure shows intercept only. '
              'Marginal exact Monte Carlo envelopes retain unresolved attempts. Synthetic '
              'Gaussian cases do not establish full-real-grid or misspecification coverage. '
              'Rendered figure still requires visual inspection.'))


if __name__=='__main__':main()
