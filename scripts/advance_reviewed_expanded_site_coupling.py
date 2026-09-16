#!/usr/bin/env python3
"""Run pinned coupling and marker-resampling stages after an audited frame handoff."""
import argparse,fcntl,json,os,shutil,subprocess,sys,time
from pathlib import Path
import psutil
from audit_joint_path_uncertainty import sha,checked


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',type=Path,required=True);ap.add_argument('--check-config',action='store_true')
    a=ap.parse_args();c=json.loads(a.config.read_text());config_sha=sha(a.config)
    def verify():
        if sha(a.config)!=config_sha:raise ValueError('Changed controller configuration')
        for p,digest in c['pins'].items():
            if sha(Path(p))!=digest:raise ValueError('Changed source: '+p)
    def live():
        try:
            p=psutil.Process(c['producer']['pid'])
            if p.create_time()!=c['producer']['created'] or p.status()==psutil.STATUS_ZOMBIE:return False
            if p.cmdline()!=c['producer']['command']:raise ValueError('Producer command changed')
            return True
        except psutil.NoSuchProcess:return False
    verify()
    if a.check_config:
        print('Configuration pins checked; producer_live=',live());return
    root=Path(c['controller_output']);root.mkdir(parents=True,exist_ok=True)
    lock=(root/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    if (root/'state.json').exists():raise FileExistsError('Review existing controller state before restart')
    state={'status':'waiting_for_frame_producer','pid':os.getpid(),'created':psutil.Process().create_time(),'config_sha256':config_sha,'started_at':time.time()}
    def save():
        tmp=root/'state.partial';tmp.write_text(json.dumps(state,indent=2)+'\n');tmp.replace(root/'state.json')
    save();os.sched_setaffinity(0,c['cpu_affinity'])
    try:
        while live():time.sleep(30)
        verify()
        handoff=Path(c['producer_receipt']);hr=json.loads(handoff.read_text())
        if hr['status']!='complete_site_rate_comparison_and_frame_handoff' or hr['config_sha256']!=c['producer_config_sha256']:raise ValueError('Upstream controller incomplete or mismatched')
        frame=Path(c['frame']);access=Path(c['accessibility']);fr=checked(frame);checked(access)
        if fr['status']!='complete_site_rate_exposure_analysis_frame' or fr['markers']!=c['markers'] or fr['sites']!=c['sites']:raise ValueError('Wrong completed frame')
        matches=[s for s in hr['stages'] if s['receipt']==str(frame/'receipt.json')]
        if len(matches)!=1 or matches[0]['receipt_sha256']!=sha(frame/'receipt.json'):raise ValueError('Frame not bound to upstream handoff')
        state['producer_receipt_sha256']=sha(handoff)
        for output in [c['fit_output'],c['resampling_output']]:
            if Path(output).exists():raise FileExistsError('Fresh outputs required: '+output)
        def run(stage,command):
            verify()
            if psutil.virtual_memory().available<c['minimum_available_memory_gib']*2**30 or shutil.disk_usage(root).free<c['minimum_free_disk_gib']*2**30:raise RuntimeError('Insufficient resource headroom')
            state['status']='running_'+stage;save()
            with (root/(stage+'.log')).open('w') as log:
                subprocess.run([sys.executable]+command,stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1'),check=True)
        plan=json.loads(Path(c['fit_template']).read_text())
        plan.update(sites=c['sites'],markers=c['markers'],analysis_sites=c['analysis_sites'],analysis_markers=c['analysis_markers'],marker_review_policy=c['marker_review_policy'],source_frame_sha256=sha(frame/'receipt.json'),source_accessibility_sha256=sha(access/'receipt.json'),memory_allowance_gb=16,output_allowance_gb=4,planning_hours=c['fit_planning_hours'],uncertainty=f"CR1 marker-cluster covariance, finite-sample correction and t({c['analysis_markers']-1}) reference; conditional on estimated rates and fixed gene trees. BH across 72 focal tests.",interpretation=c['interpretation'])
        pp=root/'fit_plan.json';pp.write_text(json.dumps(plan,indent=2)+'\n')
        run('coupling',['scripts/fit_reviewed_conditional_site_coupling.py','--frame',str(frame),'--accessibility',str(access),'--plan',str(pp),'--output',c['fit_output']])
        fits=Path(c['fit_output']);fit=checked(fits)
        if fit['status']!='complete_exploratory_conditional_site_coupling_fits' or fit['fits']!=24 or fit['markers']!=c['analysis_markers'] or fit['sites']!=c['analysis_sites']:raise ValueError('Incomplete coupling fit grid')
        plan=json.loads(Path(c['resampling_template']).read_text())
        plan.update(source_fit_receipt_sha256=sha(fits/'receipt.json'),markers=c['analysis_markers'],leave_one_marker_out_fits=24*c['analysis_markers'],memory_allowance_gb=16,output_allowance_gb=4,planning_hours=c['resampling_planning_hours'],method=f"Within-marker centering; paired multinomial resampling of {c['analysis_markers']} whole markers, 2000 draws across 24 specifications, and omission of every marker. Expanded-row verification retained.")
        pp=root/'resampling_plan.json';pp.write_text(json.dumps(plan,indent=2)+'\n')
        run('resampling',['scripts/resample_conditional_site_coupling.py','--fits',str(fits),'--plan',str(pp),'--output',c['resampling_output']])
        rs=checked(Path(c['resampling_output']))
        if rs['status']!='complete_conditional_marker_resampling' or rs['models']!=24 or rs['markers']!=c['analysis_markers'] or rs['bootstrap_fits']!=48000 or rs['leave_one_marker_out_fits']!=24*c['analysis_markers']:
            raise ValueError('Incomplete marker-resampling grid')
        state.update(status='complete_coupling_and_marker_resampling',fit_receipt_sha256=sha(fits/'receipt.json'),resampling_receipt_sha256=sha(Path(c['resampling_output'])/'receipt.json'),resampling_status=rs['status'])
    except Exception as exc:
        state.update(status='failed_requires_review',error=repr(exc));raise
    finally:
        state['updated_at']=time.time();save()

if __name__=='__main__':main()
