#!/usr/bin/env python3
"""Resume a reviewed PMSF checkpoint after a verified host reboot."""
import argparse
import fcntl
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
import psutil
from Bio import Phylo, SeqIO
from run_species_pmsf import validate_tree
from readback_whole_proteome_family_coverage import sha


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());plan_sha=sha(args.plan)
    if Path('/proc/sys/kernel/random/boot_id').read_text().strip()!=plan['recovery_boot_id']:raise ValueError('Recovery boot changed')
    folder=Path(plan['original_output']);config_path=folder/'config.json';config=json.loads(config_path.read_text())
    if (folder/'receipt.json').exists():raise FileExistsError('Original job already has a completion/review receipt')
    lock=Path('results/phylogeny/.species_pmsf.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    def verify():
        if sha(args.plan)!=plan_sha:raise ValueError('Recovery plan changed')
        for path,digest in {**plan['pins'],**config['pinned_files']}.items():
            if sha(path)!=digest:raise ValueError('Changed dependency: '+path)
    verify()
    for process in psutil.process_iter(['pid','cmdline']):
        if process.pid==os.getpid():continue
        command=process.info['cmdline'] or []
        if command==config['command']:raise RuntimeError('Native PMSF already running')
    ckp=folder/'pmsf.ckp.gz'
    if sha(ckp)!=plan['checkpoint_sha256']:raise ValueError('Checkpoint changed before recovery')
    with gzip.open(ckp,'rb') as handle:
        if not handle.readline().startswith(b'--- # IQ-TREE Checkpoint'):raise ValueError('Unexpected checkpoint')
        for _ in handle:pass
    if psutil.virtual_memory().available<plan['resources']['minimum_available_memory_gib']*2**30:raise ValueError('Insufficient memory')
    if shutil.disk_usage(folder).free<plan['resources']['minimum_free_disk_gib']*2**30:raise ValueError('Insufficient disk')
    control=Path(plan['control']);control.mkdir(parents=True,exist_ok=False)
    shutil.copy2(ckp,control/'original_checkpoint.ckp.gz')
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='')
    started=time.monotonic()
    with (control/'resume.stdout.log').open('w') as log:
        subprocess.run(config['command'],stdout=log,stderr=subprocess.STDOUT,env=env,check=True)
    verify()
    matrix_path=Path(config['command'][config['command'].index('-s')+1])
    records=list(SeqIO.parse(matrix_path,'fasta'));taxa={r.id for r in records}
    if len(taxa)!=526 or len(records)!=526:raise ValueError('Matrix taxa differ')
    validate_tree(folder/'pmsf.treefile',taxa);count=0
    for tree in Phylo.parse(folder/'pmsf.ufboot','newick'):
        names=[n.name for n in tree.get_terminals()]
        if len(names)!=526 or set(names)!=taxa:raise ValueError('Bootstrap taxon grid differs')
        count+=1
    if count!=1000 or not (folder/'pmsf.sitefreq').is_file():raise ValueError('Incomplete support/profile output')
    receipt=dict(status='complete_pmsf_execution_pending_full_audit',returncode=0,config_sha256=sha(config_path),
                 elapsed_seconds=time.monotonic()-started,elapsed_scope='This recovery process only; excludes pre-reboot work',
                 taxa=526,columns=len(records[0]),bootstrap_trees=count,
                 recovery_plan_sha256=plan_sha,original_checkpoint_sha256=plan['checkpoint_sha256'],
                 artifacts={p.name:sha(p) for p in folder.iterdir() if p.is_file()},
                 interpretation='Checkpoint-resumed execution and full bootstrap tip grids verified. Full model/profile/support audit remains required.')
    (folder/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    (control/'receipt.json').write_text(json.dumps(dict(status='complete_checkpoint_recovery_pending_full_pmsf_audit',original_receipt_sha256=sha(folder/'receipt.json'),plan_sha256=plan_sha),indent=2)+'\n')


if __name__=='__main__':main()
