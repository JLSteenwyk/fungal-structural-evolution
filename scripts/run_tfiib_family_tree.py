#!/usr/bin/env python3
"""Infer an unrooted, supported tree from verified paired TFIIB domains."""
import argparse,fcntl,json,math,shutil,subprocess,time
from pathlib import Path
from Bio import Phylo,SeqIO
from audit_busco_gene_copies import ROOT,sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    base=ROOT/'results/phylogeny/tfiib-domain-pairs-v1';r=json.loads((base/'receipt.json').read_text());rp=ROOT/'metadata/tfiib_domain_alignment_readback.json';readback=json.loads(rp.read_text())
    if r['status']!='complete_ordered_tfiib_domain_pair_alignment' or readback['status']!='complete_domain_pair_residue_readback' or readback['source_receipt_sha256']!=sha(base/'receipt.json'):raise ValueError('Verified full domain alignment required')
    for name,h in r['artifacts'].items():
        if sha(base/name)!=h:raise ValueError('Changed domain alignment')
    source=base/'paired_domains.faa';records=list(SeqIO.parse(source,'fasta'));tips={x.id for x in records}
    if len(tips)!=len(records) or len(tips)!=r['aligned_proteins']:raise ValueError('Domain input identities differ')
    exe=shutil.which('iqtree3');version=subprocess.run([exe,'--version'],capture_output=True,text=True,check=True).stdout
    a.output.mkdir(parents=True,exist_ok=True);lock=(a.output/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    prefix=a.output/'tree';command=[exe,'-s',str(source.resolve()),'-st','AA','-m','MFP','-mset','LG,WAG,JTT','-mfreq','F','-mrate','G','--alrt','1000','-B','1000','--bnni','--boot-trees','-T','4','--mem','8G','--seed','20260913','--prefix',str(prefix)]
    config={'command':command,'version':version,'executable_sha256':sha(Path(exe)),'source_receipt_sha256':sha(base/'receipt.json'),'input_sha256':sha(source),'readback_sha256':sha(rp),'tips':len(tips),'columns':r['concatenated_columns'],'script_sha256':sha(Path(__file__)),'support':'SH-aLRT 1000; ultrafast bootstrap 1000 with bootstrap-tree NNI optimization; bootstrap trees retained','interpretation':'Unrooted domain-based candidate-family topology. Identical protein/domain sequences retain distinct input identities; zero-length placements are not resolved duplication events. No species-tree rooting/reconciliation, domain-order validation or selection claim. Model/partition/individual-repeat sensitivity remains required.'}
    cp=a.output/'run_config.json'
    if cp.exists() and json.loads(cp.read_text())!=config:raise ValueError('Changed tree configuration')
    cp.write_text(json.dumps(config,indent=2)+'\n');receipt=a.output/'receipt.json'
    if receipt.exists():
        old=json.loads(receipt.read_text())
        if old['status']=='complete_supported_tfiib_domain_tree' and old['config_sha256']==sha(cp) and sha(prefix.with_suffix('.treefile'))==old['tree_sha256']:return
        raise ValueError('Existing tree receipt requires review')
    start=time.monotonic();print('Starting',len(tips),'tips',r['concatenated_columns'],'columns',flush=True)
    with (a.output/'stdout.log').open('a') as f:subprocess.run(command,stdout=f,stderr=f,check=True)
    tree_path=prefix.with_suffix('.treefile');tree=Phylo.read(tree_path,'newick');observed=[x.name for x in tree.get_terminals()]
    if len(observed)!=len(set(observed)) or set(observed)!=tips:raise ValueError('Completed tree identity universe differs')
    for branch in tree.find_clades():
        if branch.branch_length is not None and (not math.isfinite(branch.branch_length) or branch.branch_length<0):raise ValueError('Invalid tree branch length')
    boots=prefix.with_suffix('.ufboot');count=0
    for bt in Phylo.parse(boots,'newick'):
        names=[x.name for x in bt.get_terminals()]
        if len(names)!=len(set(names)) or set(names)!=tips:raise ValueError('Bootstrap identity universe differs')
        count+=1
    if count!=1000:raise ValueError('Incomplete bootstrap tree output')
    result={'status':'complete_supported_tfiib_domain_tree','config_sha256':sha(cp),'tree_sha256':sha(tree_path),'bootstrap_trees':count,'bootstrap_sha256':sha(boots),'tips':len(tips),'elapsed_seconds':time.monotonic()-start,'interpretation':config['interpretation']}
    receipt.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
