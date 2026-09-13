#!/usr/bin/env python3
"""Prepare a verified missing family tree from the identical bundled alignment rerun."""
import argparse
import hashlib
import json
import math
import subprocess
import time
from pathlib import Path
from Bio import Phylo, SeqIO

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--diagnostic',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use new isolated recovery output')
    r=json.loads((a.diagnostic/'receipt.json').read_text())
    if r['status']!='complete_isolated_alignment_reruns' or r['diagnostics']['bundled']['returncode']!=0:raise ValueError('Successful bundled rerun required')
    source=ROOT/r['source_path'];alignment=a.diagnostic/'bundled.faa'
    if sha(source)!=r['source_sha256'] or sha(alignment)!=r['diagnostics']['bundled']['alignment_sha256']:raise ValueError('Changed source alignment')
    seq={x.id:str(x.seq) for x in SeqIO.parse(source,'fasta')};records=list(SeqIO.parse(alignment,'fasta'));msa={x.id:str(x.seq) for x in records}
    if len(records)!=len(msa) or set(msa)!=set(seq) or any(msa[k].replace('-','')!=seq[k] for k in seq):raise ValueError('Alignment identity or residue mismatch')
    env=ROOT/'.cache/envs/orthofinder';module=env/'lib/python3.12/site-packages/orthofinder/tools/trim.py';exe=module.parent.parent/'bin/FastTree'
    a.output.mkdir(parents=True);trimmed=a.output/'alignment.faa';tree=a.output/'tree.nwk'
    code='from orthofinder.tools import trim; import sys; trim.main(sys.argv[1],sys.argv[2],10,0.1,500,0.75,False)'
    commands=[[str(env/'bin/python'),'-c',code,str(alignment.resolve()),str(trimmed.resolve())],[str(exe),str(trimmed.resolve())]]
    config={'diagnostic_receipt_sha256':sha(a.diagnostic/'receipt.json'),'source_sha256':sha(source),'alignment_sha256':sha(alignment),'trim_module_sha256':sha(module),'fasttree_sha256':sha(exe),'commands':commands,'script_sha256':sha(Path(__file__)),'interpretation':'Same trim parameters and FastTree invocation as production. Isolated repair; no live output replacement or completed species reconciliation claim.'}
    (a.output/'config.json').write_text(json.dumps(config,indent=2)+'\n');start=time.monotonic()
    with (a.output/'trim.log').open('w') as f:subprocess.run(commands[0],stdout=f,stderr=f,check=True)
    trimmed_records=list(SeqIO.parse(trimmed,'fasta'));retained={x.id:str(x.seq) for x in trimmed_records}
    if len(trimmed_records)!=len(retained) or set(retained)!=set(msa) or len({len(s) for s in retained.values()})!=1:raise ValueError('Trimmed identity or dimensions differ')
    # Every retained alignment column must occur in order in the untrimmed alignment.
    ids=sorted(msa);original_columns=list(zip(*(msa[k] for k in ids)));kept_columns=list(zip(*(retained[k] for k in ids)));cursor=0
    for column in kept_columns:
        while cursor<len(original_columns) and original_columns[cursor]!=column:cursor+=1
        if cursor==len(original_columns):raise ValueError('Trimmed column is not an ordered source column')
        cursor+=1
    with tree.open('w') as out,(a.output/'fasttree.log').open('w') as err:subprocess.run(commands[1],stdout=out,stderr=err,check=True)
    parsed=Phylo.read(tree,'newick');tips=[x.name for x in parsed.get_terminals()]
    if len(tips)!=len(set(tips)) or set(tips)!=set(seq):raise ValueError('Recovered gene tree tip universe differs')
    for node in parsed.find_clades():
        if node.branch_length is not None and (not math.isfinite(node.branch_length) or node.branch_length<0):raise ValueError('Invalid recovered branch length')
    result={'status':'complete_isolated_orthology_family_tree_recovery','orthogroup':r['orthogroup'],'config_sha256':sha(a.output/'config.json'),'tips':len(tips),'untrimmed_columns':len(records[0]),'trimmed_columns':len(trimmed_records[0]),'elapsed_seconds':time.monotonic()-start,'interpretation':config['interpretation'],'artifacts':{n:sha(a.output/n) for n in ['alignment.faa','tree.nwk','trim.log','fasttree.log']}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
