#!/usr/bin/env python3
"""Restore missing family artifacts before downstream OrthoFinder stages begin."""
import argparse
import json
import os
import shutil
import tempfile
from pathlib import Path
from recover_orthology_family_tree import ROOT, sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--recovery',type=Path,required=True);p.add_argument('--production-pid',type=int,required=True);p.add_argument('--receipt',type=Path,required=True)
    a=p.parse_args()
    if a.receipt.exists():raise FileExistsError('Installation receipt already exists')
    r=json.loads((a.recovery/'receipt.json').read_text())
    if r['status']!='complete_isolated_orthology_family_tree_recovery' or r['orthogroup']!='OG0001522':raise ValueError('Verified expected recovery required')
    if sha(a.recovery/'config.json')!=r['config_sha256']:raise ValueError('Changed recovery configuration')
    for n,h in r['artifacts'].items():
        if sha(a.recovery/n)!=h:raise ValueError('Changed recovered artifact')
    control=ROOT/'results/orthology/assignment-control-v2';production=ROOT/'results/orthology/core-v1/Results_full526v2/WorkingDirectory'
    config=json.loads((control/'run_config.json').read_text())
    for n,h in config['core_artifacts'].items():
        if sha(ROOT/n)!=h:raise ValueError('Protected reference core changed')
    diagnostic=json.loads((ROOT/'results/orthology/famsa-OG0001522-diagnostic-v1/receipt.json').read_text())
    if sha(production/'Sequences_ids/OG0001522.fa')!=diagnostic['source_sha256']:raise ValueError('Changed production family input')
    proc=Path('/proc')/str(a.production_pid)
    command=(proc/'cmdline').read_bytes().replace(b'\0',b' ').decode()
    if 'orthofinder' not in command or 'full526v2' not in command or '--assign' not in command:raise ValueError('Expected live assignment process required')
    log=(control/'stdout.log').read_text()
    if 'Starting MSA/Trees' not in log or 'Inferring unrooted species tree' in log or (control/'receipt.json').exists():raise ValueError('Production passed the repairable stage; reassess downstream rerun needs')
    for child in Path('/proc').iterdir():
        if not child.name.isdigit():continue
        try:args=[x.decode() for x in (child/'cmdline').read_bytes().split(b'\0') if x]
        except (FileNotFoundError,ProcessLookupError,PermissionError):continue
        if args and Path(args[0]).name in ['famsa','FastTree','fasttree','astral-pro']:
            if Path(args[0]).name=='astral-pro' or any('OG0001522' in x for x in args):raise ValueError('Active downstream stage or family writer; defer installation')
    pairs=[(a.recovery/'alignment.faa',production/'Alignments_ids/OG0001522.fa'),(a.recovery/'tree.nwk',production/'Trees_ids/OG0001522.txt')]
    if any(dest.exists() for _,dest in pairs):raise FileExistsError('Never overwrite production artifacts')
    installed=[]
    for source,dest in pairs:
        # Independent inode, complete bytes before atomic exclusive link; no partial tree exposure.
        with tempfile.NamedTemporaryFile(dir=dest.parent,prefix='.recovery-',delete=False) as f:
            temporary=Path(f.name)
            with source.open('rb') as src:shutil.copyfileobj(src,f)
            f.flush();os.fsync(f.fileno())
        try:
            if sha(temporary)!=sha(source):raise ValueError('Copy checksum differs')
            os.link(temporary,dest)
        finally:temporary.unlink(missing_ok=True)
        installed.append({'path':str(dest.relative_to(ROOT)),'sha256':sha(dest),'source_path':str(source)})
    result={'status':'complete_missing_family_artifact_installation','orthogroup':r['orthogroup'],'recovery_receipt_sha256':sha(a.recovery/'receipt.json'),'production_config_sha256':sha(control/'run_config.json'),'production_pid':a.production_pid,'production_command':command,'script_sha256':sha(Path(__file__)),'installed':installed,'interpretation':'Only two previously absent family files restored, using verified identical-binary recovery before downstream species-tree/reconciliation stages. No existing production artifacts overwritten. Full assignment and scientific validation remain pending.'}
    a.receipt.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
