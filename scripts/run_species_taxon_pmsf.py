#!/usr/bin/env python3
"""Run all native taxon/alignment/guide PMSF sensitivities serially with memory gates.

Completed old homogeneous guides are frozen conditioning inputs, checked against
identical matrix bytes and raw reports; their reuse is not a new supported tree
or retrospective original-process completion claim. Missing guides are inferred
from the full corresponding sensitivity matrix. Every PMSF mixture profile and
supported tree is refitted, followed by complete native-output readback.
"""
import argparse
import fcntl
import io
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
import psutil
from Bio import Phylo, SeqIO
from audit_species_guide import edges
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def write(path, value):
    with Path(path).open('x') as handle: handle.write(json.dumps(value,indent=2)+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);args=parser.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    bindings=dict(plan['pins']);bind(bindings,args.plan)
    closure=json.loads(Path(plan['input_completion']).read_text());assert closure['status']=='complete_verified_native_species_taxon_refit_inputs' and len(closure['services'])==2 and closure['summary']['matrices']==8
    for path,digest in closure['source_hashes'].items():bind(bindings,path,digest)
    verify(bindings);inputs=Path(plan['inputs']);source=json.loads((inputs/'receipt.json').read_text());assert str(inputs/'receipt.json') in bindings
    matrices={(r['alignment'],r['policy']):r for r in source['matrices']};assert len(matrices)==8
    out=Path(plan['output']);out.mkdir(exist_ok=True);lock=(out/'.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    config=out/'batch_config.json';expected=dict(plan_sha256=ph,source_hashes=bindings)
    if config.exists():assert json.loads(config.read_text())==expected
    else:write(config,expected)
    assert not (out/'receipt.json').exists(),'Batch already complete; inspect original completion'
    native=Path(plan['executable']);resources=plan['resources'];env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')

    def ready(memory):
        while psutil.virtual_memory().available<memory*2**30:
            print('waiting_for_available_memory_gib',memory,'observed',round(psutil.virtual_memory().available/2**30,2),flush=True);time.sleep(30);verify(bindings)
        assert shutil.disk_usage(out).free>=resources['minimum_free_disk_gib']*2**30

    def execute(command, folder, label):
        with (folder/(label+'.stdout.log')).open('a') as handle:
            process=subprocess.Popen(command,stdout=handle,stderr=subprocess.STDOUT,env=env);p=psutil.Process(process.pid)
            lp=folder/(label+'.native_launch.json');number=1
            while lp.exists():number+=1;lp=folder/(label+f'.native_launch_{number}.json')
            write(lp,dict(pid=p.pid,created=p.create_time(),cmdline=p.cmdline(),command=command,plan_sha256=ph))
            print('native_stage_started',label,p.pid,command,flush=True);code=process.wait();assert code==0,(label,code)

    def guide(label, policy):
        spec=matrices[label,policy];matrix=Path(spec['path'])/'matrix.faa';folder=out/'guides'/(label+'-'+policy);folder.mkdir(parents=True,exist_ok=True);rp=folder/'receipt.json'
        if rp.exists():
            r=json.loads(rp.read_text());assert r['plan_sha256']==ph and r['matrix_sha256']==sha(matrix)
            for name,digest in r['artifacts'].items():assert sha(folder/name)==digest
            return folder/'guide.treefile'
        cached=plan['cached_guides'].get(label+'-'+policy);provenance={}
        if cached:
            root=Path(cached['path']);r=json.loads((root/'receipt.json').read_text());assert r['status']=='complete_sensitivity_guide' and r['input_sha256']==sha(matrix) and r['taxa']==spec['taxa']
            assert r['config_sha256']==sha(cached['config'])
            command=r['command'];assert command[command.index('-m')+1]=='LG+F+G4' and command[0]==str(native)
            assert sha(command[command.index('-s')+1])==r['input_sha256']
            for name,digest in r['artifacts'].items():assert sha(root/name)==digest
            for name in ['guide.treefile','guide.iqtree','guide.log']:shutil.copyfile(root/name,folder/name)
            provenance=dict(kind='frozen_completed_homogeneous_guide_conditioning_input',source=str(root),source_receipt_sha256=sha(root/'receipt.json'),original_process_completion_revalidated=False)
        else:
            ready(resources['guide_minimum_available_memory_gib']);command=[str(native),'-s',str(matrix.resolve()),'-st','AA','-m','LG+F+G4','-T',str(resources['threads']),'--mem','32G','--seed',str(plan['seed']),'--prefix',str((folder/'guide').resolve())];execute(command,folder,'guide');provenance=dict(kind='new_full_sensitivity_homogeneous_guide',command=command)
        taxa={r.id for r in SeqIO.parse(matrix,'fasta')};assert len(taxa)==spec['taxa'];observed=edges(Phylo.read(folder/'guide.treefile','newick'),taxa)
        report=(folder/'guide.iqtree').read_text();reported=[line.strip() for line in report.splitlines() if line.startswith('(') and line.rstrip().endswith(';')]
        assert len(reported)==1 and edges(Phylo.read(io.StringIO(reported[0]),'newick'),taxa)==observed
        assert f"Input data: {spec['taxa']} sequences with {spec['columns']} amino-acid sites" in report
        assert re.search(r'Model of substitution: (\S+)',report)[1]=='LG+F+G4'
        assert math.isfinite(float(re.search(r'Log-likelihood of the tree: ([-\d.]+)',report)[1]))
        assert abs(math.fsum(observed.values())-float(re.search(r'Total tree length \(sum of branch lengths\): ([\d.]+)',report)[1]))<.000051
        r=dict(status='checked_full_taxon_sensitivity_conditioning_guide',plan_sha256=ph,matrix_sha256=sha(matrix),taxa=len(taxa),guide_provenance=provenance,artifacts={p.name:sha(p) for p in folder.iterdir() if p.is_file()},scientific_eligibility=False);write(rp,r);return folder/'guide.treefile'

    completed=[]
    for job in plan['jobs']:
        label,policy,gl=job['alignment'],job['policy'],job['guide_alignment'];spec=matrices[label,policy];matrix=Path(spec['path'])/'matrix.faa';g=guide(gl,policy);folder=out/job['label'];folder.mkdir(exist_ok=True);rp=folder/'receipt.json';audit=folder/'audit'
        if rp.exists():
            r=json.loads(rp.read_text());assert r['status']=='complete_pmsf_execution_pending_full_audit' and r['returncode']==0
            for name,digest in r['artifacts'].items():assert sha(folder/name)==digest
            c=json.loads((folder/'config.json').read_text());assert c['batch_plan_sha256']==ph and r['config_sha256']==sha(folder/'config.json')
        else:
            ready(resources['mixture_minimum_available_memory_gib'])
            shared=(Path('results/phylogeny')/'.species_pmsf.lock').open('a')
            while True:
                try:fcntl.flock(shared,fcntl.LOCK_EX|fcntl.LOCK_NB);break
                except BlockingIOError:print('waiting_for_existing_species_mixture_lock',flush=True);time.sleep(30)
            ready(resources['mixture_minimum_available_memory_gib']);verify(bindings)
            frozen=folder/'input_guide.treefile';shutil.copyfile(g,frozen);taxa={r.id for r in SeqIO.parse(matrix,'fasta')};assert {t.name for t in Phylo.read(frozen,'newick').get_terminals()}==taxa
            command=[str(native),'-s',str(matrix.resolve()),'-st','AA','-m','LG+C20+F+G4','--tree-freq',str(frozen.resolve()),'-T',str(resources['threads']),'--mem','600G','--seed',str(plan['seed']),'--alrt','1000','-B','1000','--bnni','--boot-trees','--prefix',str((folder/'pmsf').resolve())]
            config_path=folder/'config.json';gr=json.loads((g.parent/'receipt.json').read_text());run_pins={**bindings,str(matrix):sha(matrix),str(g):sha(g),str(frozen):sha(frozen),str(g.parent/'receipt.json'):sha(g.parent/'receipt.json')}
            for name,digest in gr['artifacts'].items():run_pins[str(g.parent/name)]=digest
            c=dict(command=command,pinned_files=run_pins,batch_plan_sha256=ph,resource_plan=dict(model='LG+C20+F+G4',threads=resources['threads']),alignment=label,policy=policy,guide_alignment=gl)
            if config_path.exists():assert json.loads(config_path.read_text())==c
            else:write(config_path,c)
            started=time.monotonic();execute(command,folder,'pmsf');verify(run_pins);shared.close()
            required=['pmsf.treefile','pmsf.contree','pmsf.ufboot','pmsf.sitefreq','pmsf.iqtree','pmsf.log','pmsf.splits.nex'];assert all((folder/name).is_file() for name in required)
            r=dict(status='complete_pmsf_execution_pending_full_audit',returncode=0,taxa=spec['taxa'],columns=spec['columns'],config_sha256=sha(config_path),elapsed_seconds=time.monotonic()-started,artifacts={p.name:sha(p) for p in folder.iterdir() if p.is_file()},scientific_eligibility=False);write(rp,r)
        if not (audit/'receipt.json').exists():subprocess.run([sys.executable,'scripts/audit_species_taxon_pmsf.py','--run',str(folder),'--output',str(audit)],check=True,env=env)
        a=json.loads((audit/'receipt.json').read_text());assert a['status']=='passed_taxon_sensitivity_pmsf_profile_tree_and_bootstrap_readback' and a['source_receipt_sha256']==sha(rp) and a['taxa']==spec['taxa'] and a['bootstrap_trees']==1000 and a['sites']==spec['columns']
        for name,digest in a['artifacts'].items():assert sha(audit/name)==digest
        completed.append(dict(**job,run=str(folder),run_receipt_sha256=sha(rp),audit=str(audit),audit_receipt_sha256=sha(audit/'receipt.json'),taxa=spec['taxa'],columns=spec['columns']));checkpoint=out/(job['label']+'.completed.json')
        if checkpoint.exists():assert json.loads(checkpoint.read_text())==completed[-1]
        else:write(checkpoint,completed[-1])
        print('completed_supported_taxon_sensitivity',job['label'],len(completed),'/',len(plan['jobs']),flush=True)
    assert len(completed)==16;verify(bindings);write(out/'receipt.json',dict(status='complete_native_taxon_pmsf_batch_pending_full_independent_collection_readback',plan_sha256=ph,runs=completed,source_hashes=bindings,scientific_eligibility=False,scope=plan['scope']))


if __name__=='__main__':main()
