#!/usr/bin/env python3
"""Run full resolved-tree membership readback after both native executions finish."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
from assess_small_family_output_exposure import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',required=True,type=Path)
    a=p.parse_args();plan=json.loads(a.plan.read_text());pins={str(a.plan):sha(a.plan),**plan['pins']}
    def verify():
        for path,digest in pins.items():
            if sha(path)!=digest:raise ValueError('Changed pin: '+path)
    verify();out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);start=time.time();done=[]
    def state(status,**details):
        temp=out/'state.tmp';temp.write_text(json.dumps(dict(status=status,completed=done,elapsed_seconds=time.time()-start,**details),indent=2)+'\n');temp.replace(out/'state.json')
    try:
        state('waiting_for_exact_native_controller')
        proc=Path('/proc')/str(plan['predecessor_pid'])/'stat'
        while proc.exists():
            try:fields=proc.read_text().rsplit(')',1)[1].split()
            except FileNotFoundError:break
            if fields[19]!=str(plan['predecessor_start_ticks']) or fields[0]=='Z':break
            time.sleep(20)
        execution=json.loads(Path(plan['execution_plan']).read_text());base=Path(execution['output'])
        native=json.loads((base/'state.json').read_text())
        if native['status']!='both_native_executions_complete_pending_full_output_readback' or native['plan_sha256']!=sha(plan['execution_plan']):
            raise ValueError('Native execution did not complete under pinned plan')
        inputs=Path(execution['inputs']);ir=json.loads((inputs/'receipt.json').read_text())
        if sha(inputs/'receipt.json')!=execution['pins'][str(inputs/'receipt.json')]:raise ValueError('Changed input receipt')
        for summary in ir['guides']:
            guide=summary['guide'];manifest=inputs/guide/'copied_files.tsv'
            if sha(manifest)!=summary['manifest_sha256']:raise ValueError('Changed staged manifest')
            rows={r['relative_path']:r for r in csv.DictReader(manifest.open(),delimiter='\t')}
            source=base/guide/'Source/WorkingDirectory'
            for name in ['SpeciesIDs.txt','SequenceIDs.txt','clusters_OrthoFinder.txt_id_pairs.txt']:
                if sha(source/name)!=rows['Source/WorkingDirectory/'+name]['sha256']:raise ValueError('Native source changed')
            stage=next(s for s in native['stages'] if s['guide']==guide)
            result=Path(stage['result']);resolved=result/'Resolved_Gene_Trees/Resolved_Gene_Trees.txt'
            expected_sha=stage['mandatory_artifacts']['Resolved_Gene_Trees/Resolved_Gene_Trees.txt']
            if sha(resolved)!=expected_sha:raise ValueError('Native resolved-tree output changed')
            if shutil.disk_usage(out).free<plan['resources']['minimum_free_disk_gib']*2**30:raise ValueError('Insufficient free disk')
            verify();state('checking_complete_resolved_tree_membership',guide=guide)
            target=out/(guide+'_readback.json')
            command=[sys.executable,'scripts/readback_resolved_tree_memberships.py','--source',str(source),'--results',str(result),'--output',str(target)]
            with (out/(guide+'.log')).open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
            r=json.loads(target.read_text())
            if (r['status']!='passed_complete_resolved_tree_membership_readback' or r['source_families']!=summary['families']
                    or r['source_proteins']!=summary['proteins'] or r['resolved_tree_file_sha256']!=expected_sha):
                raise ValueError('Wrong resolved-tree/source scope')
            done.append(dict(guide=guide,receipt_sha256=sha(target),resolved_trees=r['resolved_trees'],resolved_tips=r['resolved_tips']))
        verify()
        receipt=dict(status='complete_both_guide_resolved_tree_membership_readbacks',guides=done,plan_sha256=sha(a.plan),scope='Complete resolved-tree membership and branch-length checks only. Reconciliation events, HOGs and complete ortholog-table semantics remain separate validation requirements.')
        (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');state(receipt['status'])
    except Exception as exc:
        state('failed_requires_review',error=repr(exc));raise


if __name__=='__main__':main()
