#!/usr/bin/env python3
"""Advance the expanded confidence pipeline after the live prefetch completes."""
import argparse,csv,fcntl,json,subprocess,sys,time
from pathlib import Path
import numpy as np
from audit_busco_gene_copies import ROOT,sha,read_table


def process_identity(pid):
    path=Path('/proc')/str(pid)
    try:
        stat=(path/'stat').read_text().split(') ',1)[1].split()
        if stat[0]=='Z':return None
        return {'start_ticks':stat[19],'command':(path/'cmdline').read_bytes().replace(b'\0',b' ').decode().strip()}
    except FileNotFoundError:return None


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--prefetch-pid',type=int,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    catalog=ROOT/'results/structural_markers/gdm-prefetch-catalog-v1';mapping=ROOT/'results/structural_markers/gdm-expanded-v1';coordinates=ROOT/'results/structural_alphabet/coordinate-gdm-expanded-v1';prefetch=ROOT/'results/structural_pae/gdm-prefetch-v1'
    bound=ROOT/'results/structural_pae/gdm-expanded-v1';qualified=ROOT/'results/structural_alphabet/audited-gdm-expanded-v1'
    a.output.mkdir(parents=True,exist_ok=True);lock=(a.output/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    scripts=['verify_marker_catalog_mapping.py','retrieve_marker_pae.py','qualify_native_pae.py','retrieve_matched_models.py','assess_pae_sensitivity.py','compare_marker_structures.py']
    current=process_identity(a.prefetch_pid);cp=a.output/'config.json'
    if current and ('scripts/retrieve_marker_pae.py' not in current['command'] or 'gdm-prefetch-catalog-v1' not in current['command'] or 'gdm-prefetch-v1' not in current['command']):raise ValueError('PID does not identify the expected prefetch')
    dependency=json.loads(cp.read_text())['dependency_identity'] if cp.exists() else current
    config={'prefetch_pid':a.prefetch_pid,'dependency_identity':dependency,'catalog_receipt_sha256':sha(catalog/'receipt.json'),'mapping_receipt_sha256':sha(mapping/'receipt.json'),'coordinate_receipt_sha256':sha(coordinates/'receipt.json'),'producer_scripts':{name:sha(ROOT/'scripts'/name) for name in scripts},'script_sha256':sha(Path(__file__)),'mapping_bound_output':str(bound.relative_to(ROOT)),'qualified_output':str(qualified.relative_to(ROOT))}
    if cp.exists() and json.loads(cp.read_text())!=config:raise ValueError('Changed confidence pipeline configuration')
    cp.write_text(json.dumps(config,indent=2)+'\n')
    if current and dependency!=current:raise ValueError('Prefetch PID identity changed')
    last=0
    while current:
        if time.monotonic()-last>=60:print('Waiting for verified live prefetch PID',a.prefetch_pid,flush=True);last=time.monotonic()
        time.sleep(10);current=process_identity(a.prefetch_pid)
        if current and current!=dependency:raise ValueError('Prefetch PID reused before completion audit')
    pr=json.loads((prefetch/'receipt.json').read_text());models=json.loads((catalog/'model_provenance.json').read_text())
    if pr['mapping_receipt_sha256']!=config['catalog_receipt_sha256'] or pr['models_failed'] or pr['models_requested']!=len(models) or pr['models_verified']!=len(models):raise ValueError('Prefetch did not complete the full catalog')
    if sha(prefetch/'pae_manifest.json')!=pr['artifacts']['pae_manifest.json']:raise ValueError('Changed prefetch manifest')
    def stage(script,args,log):
        for name,h in config['producer_scripts'].items():
            if sha(ROOT/'scripts'/name)!=h:raise ValueError('Changed stage code '+name)
        print('Starting',script,flush=True)
        with (a.output/log).open('a') as f:subprocess.run([sys.executable,str(ROOT/'scripts'/script),*map(str,args)],cwd=ROOT,stdout=f,stderr=f,check=True)
    agreement=a.output/'catalog_mapping_agreement.json'
    if not agreement.exists():stage('verify_marker_catalog_mapping.py',['--catalog',catalog,'--mapping',mapping,'--output',agreement],'agreement.log')
    ag=json.loads(agreement.read_text())
    if ag['catalog_receipt_sha256']!=config['catalog_receipt_sha256'] or ag['mapping_receipt_sha256']!=config['mapping_receipt_sha256']:raise ValueError('Changed catalog/mapping agreement')
    if not (bound/'receipt.json').exists():stage('retrieve_marker_pae.py',['--snapshot',mapping,'--output',bound],'mapping_bound_pae.log')
    br=json.loads((bound/'receipt.json').read_text())
    if br['mapping_receipt_sha256']!=config['mapping_receipt_sha256'] or br['models_failed'] or br['models_verified']!=len(models):raise ValueError('Incomplete final mapping confidence')
    if not (qualified/'receipt.json').exists():stage('qualify_native_pae.py',['--coordinates',coordinates,'--pae',bound,'--output',qualified],'qualification.log')
    qr=json.loads((qualified/'receipt.json').read_text())
    if qr['status']!='complete_native_3di_feature_audit' or qr['pae_receipt_sha256']!=sha(bound/'receipt.json') or qr['coordinate_audit_receipt_sha256']!=config['coordinate_receipt_sha256']:raise ValueError('Qualified source linkage differs')
    for name,h in qr['artifacts'].items():
        if sha(qualified/name)!=h:raise ValueError('Changed qualification table')
    rows=read_table(qualified/'model_summary.tsv');seen=set();total=0
    for row in rows:
        key=(row['model_id'],str(row['version']))
        if key in seen:raise ValueError('Repeated qualified model')
        seen.add(key);path=ROOT/row['encoding_path']
        if sha(path)!=row['encoding_sha256']:raise ValueError('Changed qualified encoding')
        with np.load(path,allow_pickle=False) as d:
            valid=d['valid'];pae=d['feature_max_pae']
            if pae.shape!=valid.shape or len(valid)!=int(row['length']) or not np.isfinite(pae[valid]).all() or (pae[valid]<0).any() or not np.isnan(pae[~valid]).all():raise ValueError('Invalid feature PAE array')
            n=int((valid & (d['feature_min_plddt']>=70) & (pae<=10)).sum())
            if n!=int(row['valid_feature_plddt70_pae10']):raise ValueError('Joint confidence recount differs')
            total+=n
    if seen!={(x['model_id'],str(x['version'])) for x in models} or total!=qr['totals']['valid_feature_plddt70_pae10']:raise ValueError('Final qualified universe/totals differ')
    result={'status':'complete_expanded_mapping_bound_structural_confidence','config_sha256':sha(cp),'prefetch_receipt_sha256':sha(prefetch/'receipt.json'),'mapping_bound_pae_receipt_sha256':sha(bound/'receipt.json'),'qualified_receipt_sha256':sha(qualified/'receipt.json'),'models':len(seen),'joint_plddt70_pae10_valid_states':total,'interpretation':'Complete prefetch, exact catalog/final mapping agreement, cache revalidation bound to final mapping, native feature confidence qualification and per-model output/count readback. Not yet expanded structural evolutionary inference or calibration of prediction accuracy.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
