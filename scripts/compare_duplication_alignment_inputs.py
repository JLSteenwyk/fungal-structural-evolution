#!/usr/bin/env python3
"""Screen exact old/new alignment inputs; native result reuse remains separate."""
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import time
import psutil
from screen_duplication_alignment_reuse import sha, load


def signature(row):
    # Paths and coordinate-shard boundaries may differ after catalog expansion.
    fields=['model_id','version','mask','source_sha256','status','reason','sequence','original_positions','retained_residues','original_length','sha256']
    value={k:row[k] for k in fields if k in row}
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def main():
    pp=Path('metadata/duplication_exact_input_comparison_plan_20260929.json')
    plan=json.loads(pp.read_text());bindings={str(pp):sha(pp),**plan['pins']}
    launch=json.loads(Path(plan['audit_launch']).read_text())
    while psutil.pid_exists(launch['pid']):
        try:
            p=psutil.Process(launch['pid'])
            if abs(p.create_time()-launch['created'])>.01 or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
    proof=json.loads(Path(plan['audit_proof']).read_text());newroot=Path(plan['new_inputs'])
    assert proof['status']=='passed_full_duplication_alignment_input_readback' and proof['source_receipt_sha256']==sha(newroot/'receipt.json')
    bindings[plan['audit_proof']]=sha(plan['audit_proof'])
    oldmodels,oldpairs,oldbindings=load(Path(plan['old_queue']))
    newmodels,newpairs,newbindings=load(Path(plan['new_queue']))
    bindings.update(oldbindings);bindings.update(newbindings)
    def verify():
        for path,h in bindings.items():assert sha(path)==h,path
    verify();tables=[];hashed_files=0
    for name,models in [('old_inputs',oldmodels),('new_inputs',newmodels)]:
        root=Path(plan[name]);rp=root/'receipt.json';r=json.loads(rp.read_text())
        assert r['status']=='complete_duplication_alignment_input_materialization'
        queue=Path(plan['old_queue' if name=='old_inputs' else 'new_queue'])
        assert r['queue_receipt_sha256']==sha(queue/'receipt.json')
        bindings[str(rp)]=sha(rp)
        bindings.update({str(root/n):h for n,h in r['artifacts'].items()})
        assert sha(root/'inputs.jsonl')==r['artifacts']['inputs.jsonl']
        table={};counts=Counter()
        for line in (root/'inputs.jsonl').open():
            row=json.loads(line);key=row['model_id'],row['version'],row['mask']
            assert key not in table and key[:2] in models and key[2] in ['full','plddt70']
            assert row['source_sha256']==models[key[:2]]['sha256']
            if row['status']=='ready':
                assert sha(row['path'])==row['sha256'];hashed_files+=1
            else:
                assert row['status'] in ['too_few_retained_residues','source_rejected'] and 'path' not in row
            table[key]=(signature(row),row['status'])
            counts[row['mask']+':'+row['status']]+=1
        assert len(table)==2*len(models)==r['input_dispositions'] and dict(counts)==r['counts']
        tables.append(table);print('Checked exact alignment inputs',name,len(table),flush=True)
    old,new=tables;out=Path(plan['output']);out.mkdir(exist_ok=False)
    counts=Counter()
    with (out/'pair_mask_input_comparison.tsv').open('w') as handle:
        writer=csv.writer(handle,delimiter='\t',lineterminator='\n')
        writer.writerow(['pair_key','mask','disposition'])
        for pair,ends in sorted(newpairs.items()):
            if pair in oldpairs:assert ends==oldpairs[pair]
            for mask in ['full','plddt70']:
                keys=[(*end,mask) for end in ends]
                if pair not in oldpairs:status='new_pair_requires_alignment_or_input_disposition'
                elif any(old[k][0]!=new[k][0] for k in keys):status='shared_pair_changed_input_requires_alignment_or_input_disposition'
                elif all(new[k][1]=='ready' for k in keys):status='exact_ready_inputs_pending_native_settings_and_result_checks'
                else:status='exact_unavailable_inputs_pending_prior_disposition_check'
                counts[status]+=1;writer.writerow([pair,mask,status])
    assert sum(counts.values())==2*len(newpairs)==269624
    verify()
    result=dict(status='complete_exact_duplication_input_comparison_not_result_reuse',pairs=len(newpairs),pair_masks=sum(counts.values()),old_input_dispositions=len(old),new_input_dispositions=len(new),pdb_files_rehashed=hashed_files,counts=dict(counts),source_hashes=bindings,script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in out.iterdir()},scope='Full old/new input manifests compared by model/version, source bytes, confidence mask, disposition, sequence, original positions, lengths and actual PDB hash. Matching input pairs remain candidates only; executable/options, both orders, old native checkpoints and numerical discrepancy flags must be checked before result reuse. No alignment results copied. Full comparison readback pending.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'}),flush=True)


if __name__=='__main__':main()
