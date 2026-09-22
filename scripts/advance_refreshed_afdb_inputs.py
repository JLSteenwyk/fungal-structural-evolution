#!/usr/bin/env python3
"""Prepare independently runnable mapping and PAE stages for a frozen catalog."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
from assess_small_family_output_exposure import sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',required=True,type=Path)
    p.add_argument('--stage',required=True,choices=['mapping','pae'])
    a=p.parse_args();plan=json.loads(a.plan.read_text());pins={str(a.plan):sha(a.plan),**plan['pins']}
    def verify():
        for path,digest in pins.items():
            if sha(path)!=digest:raise ValueError('Changed pinned source: '+path)
    verify()
    out=Path(plan['controller_output'])/a.stage;out.mkdir(parents=True,exist_ok=False)
    started=time.time()
    def state(status,**details):
        temp=out/'state.tmp'
        temp.write_text(json.dumps(dict(status=status,stage=a.stage,elapsed_seconds=time.time()-started,**details),indent=2)+'\n')
        temp.replace(out/'state.json')
    def run(script,args,log):
        command=[sys.executable,'scripts/'+script,*map(str,args)]
        state('running',command=command)
        with (out/log).open('w') as f:subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,check=True)
    try:
        if shutil.disk_usage(out).free<plan['resources']['minimum_free_disk_gib']*2**30:
            raise ValueError('Insufficient disk headroom')
        catalog=Path(plan['catalog']);cr=json.loads((catalog/'receipt.json').read_text())
        audit=json.loads(Path(plan['catalog_readback']).read_text())
        if (audit['status']!='passed_full_frozen_inventory_marker_catalog_selection_readback'
                or audit['catalog_receipt_sha256']!=sha(catalog/'receipt.json')):
            raise ValueError('Catalog selection readback is not bound')
        if a.stage=='mapping':
            target=Path(plan['mapping_output'])
            run('map_marker_structures.py',['--inventory',plan['inventory'],'--output',target],'mapping.log')
            run('verify_marker_catalog_mapping.py',['--catalog',catalog,'--mapping',target,
                '--output',out/'catalog_mapping_agreement.json'],'agreement.log')
            mr=json.loads((target/'receipt.json').read_text())
            if mr['distinct_models']!=cr['distinct_models'] or mr['marker_proteins_linked']!=cr['marker_proteins_linked']:
                raise ValueError('Mapping coverage differs from catalog')
            status='complete_refreshed_mapping_and_catalog_agreement_pending_residue_readback'
        else:
            target=Path(plan['pae_output'])
            run('retrieve_marker_pae.py',['--snapshot',catalog,'--output',target],'pae.log')
            pr=json.loads((target/'receipt.json').read_text())
            if (pr['mapping_receipt_sha256']!=sha(catalog/'receipt.json') or pr['models_failed']
                    or pr['models_requested']!=cr['distinct_models'] or pr['models_verified']!=cr['distinct_models']
                    or sha(target/'pae_manifest.json')!=pr['artifacts']['pae_manifest.json']):
                raise ValueError('PAE prefetch incomplete; preserve successes and inspect failure manifest')
            status='complete_refreshed_catalog_pae_prefetch_pending_mapping_binding'
        verify()
        receipt=dict(status=status,plan_sha256=sha(a.plan),stage=a.stage,
                     result_receipt=str(target/'receipt.json'),result_receipt_sha256=sha(target/'receipt.json'),
                     elapsed_seconds=time.time()-started)
        (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');state(status)
    except Exception as exc:
        state('failed_requires_review',error=repr(exc));raise


if __name__=='__main__':main()
