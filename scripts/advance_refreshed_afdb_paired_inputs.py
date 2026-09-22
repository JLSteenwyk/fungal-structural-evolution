#!/usr/bin/env python3
"""Prepare refreshed AlphaFold paired inputs after full confidence readback."""
import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import psutil
from audit_busco_gene_copies import ROOT,sha,read_table
from assess_pae_sensitivity import checked_receipt


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan',required=True,type=Path)
    a=ap.parse_args();plan=json.loads(a.plan.read_text());plan_sha=sha(a.plan)
    def verify():
        if sha(a.plan)!=plan_sha:raise ValueError('Plan changed')
        for name,digest in plan['pins'].items():
            if sha(ROOT/name)!=digest:raise ValueError('Changed input: '+name)
    verify()
    control=ROOT/plan['control'];inputs=ROOT/plan['inputs'];audit_path=ROOT/plan['readback']
    if any(p.exists() for p in [control,inputs,audit_path]):raise FileExistsError('Use immutable new outputs')
    control.mkdir(parents=True)
    state=dict(status='waiting_for_exact_predecessors',started_unix=time.time(),plan_sha256=plan_sha)
    def save():
        p=control/'state.tmp';p.write_text(json.dumps(state,indent=2)+'\n');p.replace(control/'state.json')
    save()
    for dependency in plan['predecessors']:
        while True:
            try:
                p=psutil.Process(dependency['pid'])
                live=p.create_time()==dependency['create_time'] and p.status()!=psutil.STATUS_ZOMBIE
            except psutil.NoSuchProcess:live=False
            if not live:break
            time.sleep(20)
        result=json.loads((ROOT/dependency['receipt']).read_text())
        if any(result.get(k)!=v for k,v in dependency['expected'].items()):raise ValueError('Unsuccessful predecessor')
    verify()
    baseline=checked_receipt(ROOT/plan['baseline'])
    confidence_path=ROOT/plan['predecessors'][0]['receipt']
    confidence=json.loads(confidence_path.read_text())
    if confidence['plan_sha256']!=sha(ROOT/plan['confidence_plan']):
        raise ValueError('Wrong confidence controller plan')
    if confidence['mapping_receipt_sha256']!=sha(ROOT/plan['mapping']/'receipt.json'):
        raise ValueError('Confidence is bound to another mapping')
    for stage in confidence['stages']:
        if sha(ROOT/stage['result'])!=stage['sha256']:
            raise ValueError('Confidence stage receipt changed')
    encoding=checked_receipt(ROOT/plan['encodings'])
    qualified_stage=next(s for s in confidence['stages'] if s['name']=='confidence_qualification')
    if (encoding['models']!=plan['models'] or encoding['status']!='complete_native_3di_feature_audit'
            or sha(ROOT/plan['encodings']/'receipt.json')!=qualified_stage['sha256']):
        raise ValueError('Wrong confidence-qualified encoding set')
    plan['pins'][str(confidence_path.relative_to(ROOT))]=sha(confidence_path)
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='')
    commands=[('preparation',[sys.executable,'scripts/prepare_paired_phylogenetic_inputs.py','--encodings',plan['encodings'],
                '--snapshot',plan['mapping'],'--matrix',plan['matrix'],'--output',plan['inputs']]),
              ('array_readback',[sys.executable,'scripts/readback_paired_inputs_from_encodings.py','--inputs',plan['inputs'],
                                '--output',plan['readback']])]
    for stage,command in commands:
        verify()
        if shutil.disk_usage(ROOT).free<plan['resources']['minimum_free_disk_gib']*2**30:raise ValueError('Insufficient disk')
        state.update(status='running_stage',stage=stage,command=command);save()
        with (control/(stage+'.log')).open('w') as f:
            subprocess.run(command,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
    result=checked_receipt(inputs);audit=json.loads(audit_path.read_text())
    if (result['status']!='complete_paired_phylogenetic_input_preparation'
            or audit['status']!='passed_complete_paired_inputs_from_qualified_arrays_readback'
            or audit['source_receipt_sha256']!=sha(inputs/'receipt.json')
            or result['source_receipts']['encodings']['sha256']!=sha(ROOT/plan['encodings']/'receipt.json')
            or result['mask']!=baseline['mask'] or result['eligibility']!=baseline['eligibility']
            or result['source_receipts']['matrix']['sha256']!=baseline['source_receipts']['matrix']['sha256']):
        raise ValueError('Paired result binding or comparable filtering differs')
    old={r['marker']:r for r in read_table(ROOT/plan['baseline']/'marker_summary.tsv')}
    new={r['marker']:r for r in read_table(inputs/'marker_summary.tsv')}
    if set(old)!=set(new):raise ValueError('Marker universe changed')
    rows=[dict(marker=k,previous_eligible_taxa=old[k]['eligible_taxa'],expanded_eligible_taxa=new[k]['eligible_taxa'],
               previous_retained_columns=old[k]['retained_columns'],expanded_retained_columns=new[k]['retained_columns'],
               previous_status=old[k]['status'],expanded_status=new[k]['status']) for k in sorted(new)]
    path=control/'marker_coverage_change.tsv'
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    verify()
    receipt=dict(status='complete_refreshed_afdb_paired_inputs_and_independent_array_readback',plan_sha256=plan_sha,
                 input_receipt_sha256=sha(inputs/'receipt.json'),readback_sha256=sha(audit_path),
                 previous_ready_markers=baseline['ready_markers'],expanded_ready_markers=result['ready_markers'],
                 observed_paired_cells=audit['observed_paired_cells'],taxa=audit['taxa'],
                 artifacts={'marker_coverage_change.tsv':sha(path)},
                 scope='Complete refreshed AlphaFold model cohort with unchanged confidence and eligibility rules; all paired characters independently reconstructed. Coverage changes are not evolutionary effects. Tree fitting, branch uncertainty, direct-geometry benchmarking and copy/annotation sensitivities remain required.')
    (control/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    state.update(status=receipt['status']);save()


if __name__=='__main__':main()
