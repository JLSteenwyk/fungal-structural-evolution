#!/usr/bin/env python3
"""Run supported nucleotide gene-tree diagnostics for the full eligible genus set."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha, read_table
from run_paired_marker_fits import tree_edges, verify_bootstrap_tips


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['inputs','readback','information','resources','output']:
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();r=checked_receipt(a.inputs);rb=json.loads(a.readback.read_text());plan=json.loads(a.resources.read_text())
    source_hash=sha(a.inputs/'receipt.json')
    if rb['status']!='passed_full_genus_codon_diagnostic_readback' or rb['source_receipt_sha256']!=source_hash or plan['input_receipt_sha256']!=source_hash or plan['information_table_sha256']!=sha(a.information):
        raise ValueError('Input/readback/resource provenance differs')
    info=read_table(a.information);cases={x['case_id']:x for x in read_table(a.inputs/'case_summary.tsv') if x['status']=='ready_for_tree_and_divergence_diagnostics'}
    if len(info)!=len(cases) or {x['case_id'] for x in info}!=set(cases):raise ValueError('Incomplete information grid')
    ready=[x for x in info if x['status']=='ready_for_supported_tree_diagnostic']
    if len(ready)!=plan['ready_cases']:raise ValueError('Resource case count differs')
    executable=shutil.which('iqtree3')
    if executable is None or sha(Path(executable))!='40424ccdb1d79c304641f910cb6c172ebb670214e50958352018e3ff9906ab8f':raise ValueError('Pinned IQ-TREE build required')
    available=int(next(l.split()[1] for l in Path('/proc/meminfo').read_text().splitlines() if l.startswith('MemAvailable:')))*1024
    if available<16*2**30 or shutil.disk_usage(a.output.parent).free<plan['output_allowance_gb']*10**9:raise RuntimeError('Insufficient headroom')
    a.output.mkdir(parents=True,exist_ok=True);lock=(a.output/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    config={'input_receipt_sha256':source_hash,'readback_sha256':sha(a.readback),'information_sha256':sha(a.information),'resource_plan_sha256':sha(a.resources),'script_sha256':sha(Path(__file__)),'tree_helper_sha256':sha(Path(__file__).with_name('run_paired_marker_fits.py')),'executable':executable,'executable_sha256':sha(Path(executable)),'workers':4,'model':'GTR+F+G4','support_replicates':1000,'ready_cases':len(ready),'interpretation':'Unrooted within-genus/code nucleotide topology diagnostics; all codon positions under one GTR+F+G4 model. Genetic codes remain case metadata for later codon models. No dN/dS, selection or saturation inference.'}
    cp=a.output/'config.json'
    if cp.exists() and json.loads(cp.read_text())!=config:raise ValueError('Changed run configuration')
    cp.write_text(json.dumps(config,indent=2)+'\n')
    shutil.copyfile(a.information,a.output/'information.tsv')
    def run(row):
        case=row['case_id'];folder=a.output/case;folder.mkdir(exist_ok=True)
        alignment=a.inputs/case/'codons.fna';records=list(SeqIO.parse(alignment,'fasta'));taxa={x.id for x in records}
        if len(records)!=len(taxa) or len(taxa)!=int(row['taxa']) or any(len(x.seq)!=int(row['nucleotide_columns']) for x in records):raise ValueError('Alignment grid differs')
        seed=int(hashlib.sha256(case.encode()).hexdigest()[:8],16)%2147483646+1
        prefix=folder/'tree';command=[executable,'-s',str(alignment.resolve()),'-st','DNA','-m','GTR+F+G4','-T','1','--mem','2G','--seed',str(seed),'-keep-ident','--alrt','1000','-B','1000','--bnni','--boot-trees','--prefix',str(prefix.resolve())]
        rc={'command':command,'alignment_sha256':sha(alignment),'parent_config_sha256':sha(cp),'marker_copy_caveat':row['marker_copy_caveat'],'translation_table':cases[case]['translation_table']};rp=folder/'config.json'
        if rp.exists() and json.loads(rp.read_text())!=rc:raise ValueError('Changed case configuration')
        rp.write_text(json.dumps(rc,indent=2)+'\n');done=folder/'receipt.json'
        if done.exists():
            old=json.loads(done.read_text())
            if old['config_sha256']!=sha(rp):raise ValueError('Changed completed case configuration')
            for name,digest in old['artifacts'].items():
                if sha(folder/name)!=digest:raise ValueError('Changed completed output')
            return {'case_id':case,'receipt_sha256':sha(done)}
        start=time.monotonic();env=os.environ.copy();env['OPENBLAS_NUM_THREADS']=env['OMP_NUM_THREADS']='1'
        with (folder/'stdout.log').open('a') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=env,check=True)
        edges=tree_edges(folder/'tree.treefile',taxa);boots=verify_bootstrap_tips(folder/'tree.ufboot',taxa,1000)
        if sha(alignment)!=rc['alignment_sha256']:raise ValueError('Alignment changed during inference')
        result={'status':'complete_supported_nucleotide_tree_pending_full_audit','case_id':case,'taxa':len(taxa),'branches':len(edges),'bootstrap_trees':boots,'config_sha256':sha(rp),'elapsed_seconds':time.monotonic()-start,'artifacts':{name:sha(folder/name) for name in ['tree.treefile','tree.ufboot','tree.iqtree','tree.log']},'interpretation':'Execution, exact tip grids and finite nonnegative tree edges checked. Full report/model/support audit, codon divergence, alignment sensitivity and gene-copy review remain required.'}
        done.write_text(json.dumps(result,indent=2)+'\n');print('Completed',case,len(taxa),flush=True)
        return {'case_id':case,'receipt_sha256':sha(done)}
    results=[]
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures=[pool.submit(run,row) for row in ready]
        try:
            for future in as_completed(futures):results.append(future.result())
        except Exception:
            for future in futures:future.cancel()
            raise
    if sha(a.inputs/'receipt.json')!=source_hash or sha(Path(__file__))!=config['script_sha256']:raise ValueError('Inputs or producer changed')
    result={'status':'complete_genus_nucleotide_tree_execution_pending_full_audit','config_sha256':sha(cp),'completed_cases':len(results),'case_receipts':sorted(results,key=lambda r:r['case_id']),'excluded_cases':len(info)-len(ready),'interpretation':'Complete supported topology execution for this information-screened diagnostic set, not a final species tree, selection eligibility or selection result.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print('Completed all',len(results),'cases',flush=True)


if __name__=='__main__':main()
