#!/usr/bin/env python3
"""Finish remaining annotation shards and catalog them after pinned producers exit."""
import argparse,json,os,subprocess,sys,time,shutil
from pathlib import Path
import psutil
from prepare_pfam import digest


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--config',type=Path,required=True);ap.add_argument('--check-config',action='store_true');a=ap.parse_args();c=json.loads(a.config.read_text());config_sha=digest(a.config)
    def verify():
        if digest(a.config)!=config_sha:raise ValueError('Changed controller configuration')
        for path,sha in c['pins'].items():
            if digest(Path(path))!=sha:raise ValueError('Changed source '+path)
    def live(spec):
        try:
            p=psutil.Process(spec['pid'])
            if p.create_time()!=spec['created'] or p.status()==psutil.STATUS_ZOMBIE:return False
            if p.cmdline()!=spec['command']:raise ValueError('Producer command changed')
            return True
        except psutil.NoSuchProcess:return False
    verify()
    if a.check_config:print('Pins checked; producer states',[live(p) for p in c['producers']]);return
    root=Path(c['controller_output']);root.mkdir(parents=True,exist_ok=False);state={'status':'waiting','pid':os.getpid(),'created':psutil.Process().create_time(),'config_sha256':config_sha,'started_at':time.time(),'completed_stages':[]}
    def save():
        p=root/'state.partial';p.write_text(json.dumps(state,indent=2)+'\n');p.replace(root/'state.json')
    save();os.sched_setaffinity(0,c['cpu_affinity'])
    try:
        while any(live(p) for p in c['producers']):time.sleep(30)
        verify();template=json.loads(Path(c['plan_template']).read_text());search=Path(template['search']);raw=json.loads((search/'receipt.json').read_text());initial=Path(c['initial_snapshot']);initial_readback=Path(c['initial_readback']);ir=json.loads((initial/'receipt.json').read_text());ar=json.loads((initial_readback/'receipt.json').read_text())
        if raw['status']!='complete_raw_domain_search' or ar['status']!='passed_full_pfam_annotation_raw_field_readback' or ar['annotation_receipt_sha256']!=digest(initial/'receipt.json'):raise ValueError('Prerequisite completion failure')
        done={x['chunk'] for x in ir['shards']};allchunks={x['chunk'] for x in raw['chunks']}
        if len(allchunks)!=64 or len(done)!=52 or not done<allchunks:raise ValueError('Unexpected chunk universe')
        template['chunks']=[{'chunk':key,'receipt_sha256':digest(search/(key+'.receipt.json')),'rows':next(x['domain_hit_rows'] for x in raw['chunks'] if x['chunk']==key)} for key in sorted(allchunks-done)]
        template['output']=c['remaining_snapshot'];template['created_at']=time.time();template['basis']='Remaining 12 complete raw shards after initial 52-shard audited snapshot; same parser and query validation.';template['pins'][str(search/'receipt.json')]=digest(search/'receipt.json')
        pp=root/'remaining_annotation_plan.json';pp.write_text(json.dumps(template,indent=2)+'\n')
        def run(name,command):
            verify()
            if psutil.virtual_memory().available<32*2**30 or shutil.disk_usage(root).free<50*2**30:raise RuntimeError('Insufficient resource headroom')
            state['status']='running_'+name;save()
            with (root/(name+'.log')).open('w') as log:subprocess.run([sys.executable]+command,stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'),check=True)
            state['completed_stages'].append(name);save()
        run('remaining_annotation',['scripts/annotate_completed_full_pfam_chunks.py','--plan',str(pp)])
        run('remaining_readback',['scripts/readback_full_pfam_annotation_chunks.py','--annotations',c['remaining_snapshot'],'--plan',str(pp),'--output',c['remaining_readback']])
        run('catalog',['scripts/combine_full_pfam_annotation_snapshots.py','--search',str(search),'--pfam-receipt',template['pfam_receipt'],'--snapshots',str(initial),c['remaining_snapshot'],'--readbacks',str(initial_readback),c['remaining_readback'],'--output',c['catalog_output']])
        catalog=Path(c['catalog_output'])/'receipt.json';rr=json.loads(catalog.read_text())
        if rr['status']!='complete_audited_additional_full_proteome_annotation_catalog' or rr['chunks']!=64:raise ValueError('Incomplete catalog')
        state.update(status='complete_full_pfam_annotation_catalog_handoff',catalog_receipt_sha256=digest(catalog))
    except Exception as exc:state.update(status='failed_requires_review',error=repr(exc));raise
    finally:state['updated_at']=time.time();save()

if __name__=='__main__':main()
