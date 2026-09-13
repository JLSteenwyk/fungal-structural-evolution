#!/usr/bin/env python3
"""Run supported unrooted trees for each repeat with identical gene sampling."""
import argparse,fcntl,json,math,shutil,subprocess,time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from Bio import Phylo,SeqIO
from audit_busco_gene_copies import ROOT,sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    base=ROOT/'results/phylogeny/tfiib-repeat-sensitivity-inputs-v1';r=json.loads((base/'receipt.json').read_text())
    if r['status']!='complete_tfiib_repeat_sensitivity_inputs':raise ValueError('Complete repeat inputs required')
    for name,h in r['artifacts'].items():
        if sha(base/name)!=h:raise ValueError('Changed repeat input')
    a.output.mkdir(parents=True,exist_ok=True);lock=(a.output/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    exe=shutil.which('iqtree3');version=subprocess.run([exe,'--version'],capture_output=True,text=True,check=True).stdout;source_hash=sha(Path(__file__))
    def run(repeat):
        number=repeat['repeat_order'];folder=a.output/('repeat'+str(number));folder.mkdir(exist_ok=True)
        source=base/repeat['path'];tips={x.id for x in SeqIO.parse(source,'fasta')}
        if len(tips)!=repeat['proteins']:raise ValueError('Repeat input identity mismatch')
        prefix=folder/'tree';command=[exe,'-s',str(source.resolve()),'-st','AA','-m','MFP','-mset','LG,WAG,JTT','-mfreq','F','-mrate','G','--alrt','1000','-B','1000','--bnni','--boot-trees','-T','2','--mem','4G','--seed',str(20260913+number),'--prefix',str(prefix)]
        config={'repeat_order':number,'command':command,'executable_sha256':sha(Path(exe)),'version':version,'input_sha256':sha(source),'source_receipt_sha256':sha(base/'receipt.json'),'tips':len(tips),'columns':repeat['columns'],'script_sha256':source_hash,'support':'SH-aLRT 1000; UFBoot 1000 with NNI refinement and retained bootstrap trees','interpretation':'Unrooted positional-repeat topology; comparisons with paired-domain and other-repeat trees remain required. Identical inputs and zero-length branches are not resolved gene duplications.'}
        cp=folder/'run_config.json'
        if cp.exists() and json.loads(cp.read_text())!=config:raise ValueError('Changed repeat tree configuration')
        cp.write_text(json.dumps(config,indent=2)+'\n');rp=folder/'receipt.json'
        if rp.exists():
            old=json.loads(rp.read_text())
            if old['status']=='complete_supported_tfiib_repeat_tree' and old['config_sha256']==sha(cp) and sha(prefix.with_suffix('.treefile'))==old['tree_sha256']:return old
            raise ValueError('Existing repeat output requires review')
        start=time.monotonic();print('Starting repeat',number,len(tips),'tips',flush=True)
        with (folder/'stdout.log').open('a') as f:subprocess.run(command,stdout=f,stderr=f,check=True)
        treefile=prefix.with_suffix('.treefile');tree=Phylo.read(treefile,'newick');names=[x.name for x in tree.get_terminals()]
        if len(names)!=len(set(names)) or set(names)!=tips:raise ValueError('Repeat tree identities differ')
        for c in tree.find_clades():
            if c.branch_length is not None and (not math.isfinite(c.branch_length) or c.branch_length<0):raise ValueError('Invalid repeat branch length')
        boots=prefix.with_suffix('.ufboot');count=0
        for bt in Phylo.parse(boots,'newick'):
            names=[x.name for x in bt.get_terminals()]
            if len(names)!=len(set(names)) or set(names)!=tips:raise ValueError('Repeat bootstrap identities differ')
            count+=1
        if count!=1000:raise ValueError('Incomplete repeat bootstrap set')
        out={'status':'complete_supported_tfiib_repeat_tree','repeat_order':number,'config_sha256':sha(cp),'tree_sha256':sha(treefile),'bootstrap_sha256':sha(boots),'bootstrap_trees':count,'tips':len(tips),'elapsed_seconds':time.monotonic()-start}
        rp.write_text(json.dumps(out,indent=2)+'\n');print('Completed repeat',number,flush=True);return out
    with ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(run,r['repeats']))
    (a.output/'receipt.json').write_text(json.dumps({'status':'complete_tfiib_repeat_tree_runs','source_receipt_sha256':sha(base/'receipt.json'),'repeats':results,'interpretation':'Computational tree completion only; repeat concordance, rooting, reconciliation and model adequacy remain to assess.'},indent=2)+'\n')

if __name__=='__main__':main()
