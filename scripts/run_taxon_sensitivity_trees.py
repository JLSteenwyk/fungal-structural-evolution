#!/usr/bin/env python3
"""Fit homogeneous-model guide trees for all frozen taxon-identity sensitivities."""
import argparse,concurrent.futures,fcntl,json,math,shutil,subprocess,time
from pathlib import Path
from Bio import Phylo,SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT,sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    source=checked_receipt(a.inputs);exe=shutil.which('iqtree3')
    if not exe:raise FileNotFoundError('iqtree3')
    version=subprocess.run([exe,'--version'],capture_output=True,text=True,check=True).stdout
    a.output.mkdir(parents=True,exist_ok=True);lock=(a.output/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    config={'source_receipt_sha256':sha(a.inputs/'receipt.json'),'script_sha256':sha(Path(__file__)),'executable_sha256':sha(Path(exe)),'version':version,'workers':2,'threads_per_tree':8,'memory_per_tree':'32G','model':'LG+F+G4','seed':20260913,'support':'None: unpartitioned homogeneous-model guide sensitivity only'}
    cp=a.output/'config.json'
    if cp.exists() and json.loads(cp.read_text())!=config:raise ValueError('Changed configuration')
    cp.write_text(json.dumps(config,indent=2)+'\n')
    def run(item):
        name=item['alignment']+'-'+item['policy'];folder=a.output/name;folder.mkdir(exist_ok=True);rp=folder/'receipt.json';matrix=a.inputs/item['path']
        if sha(matrix)!=item['sha256']:raise ValueError('Changed matrix')
        if rp.exists():
            r=json.loads(rp.read_text())
            if r['status']!='complete_sensitivity_guide' or r['config_sha256']!=sha(cp):raise ValueError('Invalid cached tree')
            for n,h in r['artifacts'].items():
                if sha(folder/n)!=h:raise ValueError('Changed cached artifact')
            return r
        prefix=folder/'guide';command=[exe,'-s',str(matrix.resolve()),'-st','AA','-m','LG+F+G4','-T','8','--mem','32G','--seed','20260913','--prefix',str(prefix.resolve())]
        started=time.monotonic()
        with (folder/'stdout.log').open('a') as f:subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,check=True)
        treepath=folder/'guide.treefile';tree=Phylo.read(treepath,'newick');tips=[t.name for t in tree.get_terminals()];expected={r.id for r in SeqIO.parse(matrix,'fasta')}
        if len(tips)!=len(set(tips)) or set(tips)!=expected or len(tips)!=item['taxa']:raise ValueError('Tree identities differ')
        if any(c.branch_length is not None and (not math.isfinite(c.branch_length) or c.branch_length<0) for c in tree.find_clades()):raise ValueError('Invalid branch length')
        r={'status':'complete_sensitivity_guide','name':name,'config_sha256':sha(cp),'input_sha256':sha(matrix),'command':command,'taxa':len(tips),'elapsed_seconds':time.monotonic()-started,'artifacts':{p.name:sha(p) for p in folder.iterdir() if p.is_file()}}
        rp.write_text(json.dumps(r,indent=2)+'\n');print(name,'completed',flush=True);return r
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(run,source['matrices']))
    receipt={'status':'complete_taxon_identity_guide_sensitivities','config_sha256':sha(cp),'results':results,'interpretation':'Four guide trees fitted under the same homogeneous model as the full-data guides. Taxon/aligner sensitivity only; no support, species delimitation or final species phylogeny claim.'}
    (a.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':main()
