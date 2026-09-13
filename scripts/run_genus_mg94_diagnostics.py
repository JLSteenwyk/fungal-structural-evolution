#!/usr/bin/env python3
"""Run the full MG94 diagnostic queue as per-case nucleotide trees complete."""
import argparse
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time
from Bio import Phylo, SeqIO
from audit_genus_codon_trees import sha, table, split_map


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ['inputs','trees','information','resources','code-check','hyphy-install','hyphy-source','installed-manifest','output']:
        parser.add_argument('--'+key,required=True,type=Path)
    parser.add_argument('--tree-producer-pid',required=True,type=int)
    args=parser.parse_args()
    plan=json.loads(args.resources.read_text());code_check=json.loads(args.code_check.read_text())
    if sha(args.inputs/'receipt.json')!=plan['source_input_receipt_sha256'] or sha(args.information)!=plan['information_table_sha256'] or sha(args.code_check)!=plan['code_check_sha256']:
        raise ValueError('Changed input or resource provenance')
    if code_check['status']!='passed_executed_hyphy_fungal_code_table_check':raise ValueError('Code validation required')
    for name,digest in json.loads(args.installed_manifest.read_text()).items():
        if sha(args.hyphy_install/name)!=digest:raise ValueError('Changed HyPhy installation')
    exe=args.hyphy_install/'bin/hyphy'
    model=args.hyphy_source/'tests/hbltests/libv3/support/FitMG94.bf'
    if sha(exe)!=code_check['executable_sha256']:raise ValueError('Changed tested executable')
    if shutil.disk_usage(args.output.parent).free<plan['output_allowance_gb']*10**9:raise RuntimeError('Insufficient disk headroom')
    available=int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))*1024
    if available<16*2**30:raise RuntimeError('Insufficient RAM headroom')
    info={r['case_id']:r for r in table(args.information) if r['status']=='ready_for_supported_tree_diagnostic'}
    cases={r['case_id']:r for r in table(args.inputs/'case_summary.tsv')}
    if len(info)!=plan['cases']:raise ValueError('Case count differs')
    tree_config=args.trees/'config.json'
    if json.loads(tree_config.read_text())['input_receipt_sha256']!=plan['source_input_receipt_sha256']:raise ValueError('Tree inputs differ')
    producer=Path('/proc')/str(args.tree_producer_pid)
    if not producer.exists():raise RuntimeError('Specified tree producer is not live at launch')
    producer_start=producer.joinpath('stat').read_text().split()[21]
    producer_cmd=producer.joinpath('cmdline').read_bytes()
    if b'run_genus_codon_trees.py' not in producer_cmd:raise ValueError('Wrong producer process')
    args.output.mkdir(parents=True,exist_ok=True)
    lock=(args.output/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    config={'resource_plan_sha256':sha(args.resources),'code_check_sha256':sha(args.code_check),'input_receipt_sha256':sha(args.inputs/'receipt.json'),'information_sha256':sha(args.information),'tree_config_sha256':sha(tree_config),'installed_manifest_sha256':sha(args.installed_manifest),'executable_sha256':sha(exe),'model_source_sha256':sha(model),'script_sha256':sha(Path(__file__)),'tree_helper_sha256':sha(Path(__file__).with_name('audit_genus_codon_trees.py')),'tree_producer_pid':args.tree_producer_pid,'tree_producer_start_ticks':producer_start,'workers':4,'cases':len(info),'model':'MG94xREV global CF3x4','selection_lrt':False,'interpretation':plan['interpretation']}
    cp=args.output/'config.json'
    if cp.exists() and json.loads(cp.read_text())!=config:raise ValueError('Changed queue configuration')
    cp.write_text(json.dumps(config,indent=2)+'\n')
    def run(case):
        folder=args.output/case;folder.mkdir(exist_ok=True)
        source=args.trees/case;tr=json.loads((source/'receipt.json').read_text());tc=json.loads((source/'config.json').read_text())
        if tr['status']!='complete_supported_nucleotide_tree_pending_full_audit' or tr['config_sha256']!=sha(source/'config.json') or tc['parent_config_sha256']!=sha(tree_config):raise ValueError('Changed tree execution provenance')
        for name,digest in tr['artifacts'].items():
            if sha(source/name)!=digest:raise ValueError('Changed tree output')
        alignment=args.inputs/case/'codons.fna'
        if sha(alignment)!=tc['alignment_sha256']:raise ValueError('Changed alignment')
        records=list(SeqIO.parse(alignment,'fasta'));taxa={r.id for r in records}
        if len(records)!=len(taxa) or len(taxa)!=int(info[case]['taxa']):raise ValueError('Input grid differs')
        tree=Phylo.read(source/'tree.treefile','newick');original_splits=set(split_map(tree,taxa))
        for i,node in enumerate(tree.get_nonterminals()):node.name='InternalNode'+str(i);node.confidence=None
        treefile=folder/'tree.nwk'
        Phylo.write(tree,treefile,'newick',format_branch_length='%.10f')
        if set(split_map(Phylo.read(treefile,'newick'),taxa))!=original_splits:raise ValueError('Topology changed during label removal')
        seed=int(hashlib.sha256(case.encode()).hexdigest()[:8],16)%2147483646+1
        code=code_check['mapping'][cases[case]['translation_table']]
        command=[str(exe.resolve()),'CPU=1','ENV=RANDOM_SEED='+str(seed)+';',str(model.resolve()),'--code',code,'--alignment',str(alignment.resolve()),'--tree',str(treefile.resolve()),'--type','global','--frequencies','CF3x4','--lrt','No','--output',str((folder/'fit.json').resolve()),'--save-fit',str((folder/'fit.bf').resolve())]
        rc={'command':command,'parent_config_sha256':sha(cp),'source_tree_receipt_sha256':sha(source/'receipt.json'),'alignment_sha256':sha(alignment),'tree_sha256':sha(treefile),'translation_table':cases[case]['translation_table'],'hyphy_code':code,'marker_copy_caveat':cases[case]['marker_copy_caveat'],'independent_tree_audit':'required_before_interpretation'}
        rp=folder/'config.json'
        if rp.exists() and json.loads(rp.read_text())!=rc:raise ValueError('Changed case configuration')
        rp.write_text(json.dumps(rc,indent=2)+'\n')
        done=folder/'receipt.json'
        if done.exists():
            old=json.loads(done.read_text())
            if old['config_sha256']!=sha(rp):raise ValueError('Changed completed case')
            for name,digest in old['artifacts'].items():
                if sha(folder/name)!=digest:raise ValueError('Changed completed fit')
            return {'case_id':case,'receipt_sha256':sha(done)}
        start=time.monotonic();env=os.environ.copy();env['OMP_NUM_THREADS']=env['OPENBLAS_NUM_THREADS']='1'
        with (folder/'fit.log').open('w') as stream:subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,env=env,check=True)
        result=json.loads((folder/'fit.json').read_text());fit=result['fits']['Standard MG94']
        if result['input']['number of sequences']!=len(taxa) or result['input']['number of sites']!=int(cases[case]['retained_codon_columns']):raise ValueError('Fit input dimensions differ')
        if not math.isfinite(fit['Log Likelihood']):raise ValueError('Nonfinite fit')
        values=result['branch attributes']['0']
        for value in values.values():
            for label in ['Standard MG94','synonymous','nonsynonymous']:
                if not math.isfinite(value[label]) or value[label]<0:raise ValueError('Invalid fitted branch component')
        if sha(alignment)!=rc['alignment_sha256'] or sha(model)!=config['model_source_sha256']:raise ValueError('Source changed during fit')
        out={'status':'complete_mg94_execution_pending_independent_audit','case_id':case,'config_sha256':sha(rp),'elapsed_seconds':time.monotonic()-start,'log_likelihood':fit['Log Likelihood'],'reported_omega':fit['Rate Distributions']['non-synonymous/synonymous rate ratio'],'artifacts':{name:sha(folder/name) for name in ['fit.json','fit.bf','fit.log','tree.nwk']},'interpretation':'Global-omega diagnostic on a fixed nucleotide topology. Independent tree/fit audit, component scaling and numerical/biological review pending. No selection test.'}
        done.write_text(json.dumps(out,indent=2)+'\n');print('Completed',case,flush=True)
        return {'case_id':case,'receipt_sha256':sha(done)}
    pending=set(info);active={};completed=[]
    with ThreadPoolExecutor(max_workers=4) as pool:
        while pending or active:
            ready=[case for case in sorted(pending) if (args.trees/case/'receipt.json').exists()]
            for case in ready[:4-len(active)]:active[pool.submit(run,case)]=case;pending.remove(case)
            if active:
                done,_=wait(active,timeout=10,return_when=FIRST_COMPLETED)
                for future in done:completed.append(future.result());del active[future]
            elif pending:
                if not producer.exists() or producer.joinpath('stat').read_text().split()[21]!=producer_start:raise RuntimeError('Tree producer stopped with unresolved pending cases')
                time.sleep(10)
    if sha(Path(__file__))!=config['script_sha256']:raise ValueError('Producer changed during execution')
    result={'status':'complete_full_genus_mg94_execution_pending_audit','completed_cases':len(completed),'config_sha256':sha(cp),'case_receipts':sorted(completed,key=lambda r:r['case_id']),'interpretation':plan['interpretation']}
    (args.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print('Completed full MG94 queue',len(completed),flush=True)


if __name__=='__main__':main()
