#!/usr/bin/env python3
"""Check every materialized alignment PDB against independently validated coordinates."""
import argparse
from collections import Counter
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time
import psutil
from Bio.Data.PDBData import protein_letters_3to1
from screen_duplication_alignment_reuse import sha


def check_entry(row, record, mask, blob=None):
    assert row['model_id'] == record['model_id'] and row['version'] == record['version']
    assert row['mask'] == mask and row['source_sha256'] == record['source_sha256']
    if record['status'] != 'validated':
        assert row['status'] == 'source_rejected' and row['reason'] == record['reason']
        assert 'path' not in row and blob is None
        return 0
    positions = [i for i,b in enumerate(record['ca_plddt'],1) if mask == 'full' or b >= 70]
    sequence = ''.join(record['sequence'][i-1] for i in positions)
    assert row['original_positions'] == positions and row['sequence'] == sequence
    assert row['retained_residues'] == len(positions) and row['original_length'] == record['length']
    if len(positions) < 3:
        assert row['status'] == 'too_few_retained_residues' and 'path' not in row and blob is None
        return 0
    assert row['status'] == 'ready' and blob is not None
    assert hashlib.sha256(blob).hexdigest() == row['sha256']
    lines = blob.decode('ascii').splitlines()
    assert lines[-2:] == ['TER','END'] and len(lines) == len(positions)+2
    for serial,(line,pos,aa) in enumerate(zip(lines[:-2],positions,sequence),1):
        assert len(line) == 80 and line[:6] == 'ATOM  ' and int(line[6:11]) == serial
        assert line[12:16].strip() == 'CA' and line[16] == ' ' and line[21] == 'A'
        assert int(line[22:26]) == pos and line[26] == ' ' and protein_letters_3to1[line[17:20]] == aa
        assert float(line[54:60]) == 1. and line[76:78].strip() == 'C'
        for value,(start,end) in zip(record['ca_xyz'][pos-1],[(30,38),(38,46),(46,54)]):
            parsed = float(line[start:end])
            assert math.isfinite(parsed) and abs(parsed-value) <= .000501
        confidence = float(line[60:66])
        assert math.isfinite(confidence) and abs(confidence-record['ca_plddt'][pos-1]) <= .005001
    return len(blob)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--plan',type=Path,required=True)
    args=parser.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    bindings={str(args.plan):ph,**plan['pins']}
    launch=json.loads(Path(plan['launch']).read_text())
    while psutil.pid_exists(launch['pid']):
        try:
            p=psutil.Process(launch['pid'])
            if abs(p.create_time()-launch['created'])>.01 or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
    producer=json.loads(Path(plan['producer_plan']).read_text());assert sha(plan['producer_plan'])==launch['plan_sha256']
    root=Path(producer['output']);rp=root/'receipt.json';receipt=json.loads(rp.read_text())
    assert receipt['status']=='complete_duplication_alignment_input_materialization' and receipt['plan_sha256']==sha(plan['producer_plan'])
    coordinates=Path(producer['coordinates']);crp=coordinates/'receipt.json';cr=json.loads(crp.read_text())
    auditdir=Path(producer['readback']);arp=auditdir/'receipt.json';ar=json.loads(arp.read_text())
    assert sha(crp)==receipt['coordinate_receipt_sha256']==ar['producer_receipt_sha256']
    assert sha(arp)==receipt['readback_receipt_sha256'] and ar['status']=='passed_duplication_exported_ca_readback'
    bindings.update(producer['pins']);bindings.update({str(rp):sha(rp),str(crp):sha(crp),str(arp):sha(arp)})
    bindings.update({str(root/n):h for n,h in receipt['artifacts'].items()})
    def verify():
        for path,h in bindings.items():assert sha(path)==h,path
    verify();seen=set();counts=Counter();bytes_checked=0
    with (root/'inputs.jsonl').open() as manifest:
        for i,shard in enumerate(cr['shards']):
            source=coordinates/shard['output'];assert sha(source)==shard['output_sha256']
            proofpath=auditdir/(source.name+'.readback.json')
            assert sha(proofpath)==ar['proofs'][proofpath.name]
            assert json.loads(proofpath.read_text())['source_sha256']==shard['output_sha256']
            with gzip.open(source,'rt') as handle:
                for line in handle:
                    record=json.loads(line);key=record['model_id'],record['version'];assert key not in seen;seen.add(key)
                    name=hashlib.sha256(json.dumps(key,separators=(',',':')).encode()).hexdigest()
                    for mask in ['full','plddt70']:
                        row=json.loads(next(manifest));assert row['coordinate_shard']==source.name
                        blob=None
                        if row['status']=='ready':
                            path=root/'pdb'/name[:2]/(name+'-'+mask+'.pdb');assert row['path']==str(path);blob=path.read_bytes()
                        bytes_checked+=check_entry(row,record,mask,blob)
                        counts[mask+':'+row['status']]+=1
            assert sha(source)==shard['output_sha256']
            print('Checked duplication alignment-input shards',i+1,'/',len(cr['shards']),flush=True)
        assert not manifest.read()
    assert len(seen)==receipt['models']==cr['models']==ar['models']==276682
    assert sum(counts.values())==receipt['input_dispositions']==553364 and dict(counts)==receipt['counts']
    assert bytes_checked==receipt['pdb_bytes']
    verify()
    proof=dict(status='passed_full_duplication_alignment_input_readback',models=len(seen),input_dispositions=sum(counts.values()),counts=dict(counts),pdb_bytes=bytes_checked,source_receipt_sha256=sha(rp),plan_sha256=ph,script_sha256=sha(__file__),producer_terminal_state=state,scope='Every full and pLDDT70 disposition, original position, sequence and written C-alpha coordinate checked against the fully audited source records. PDB rounding tolerances retained. No pair alignment, PAE qualification or biological effect inference.')
    with Path(plan['proof']).open('x') as handle:json.dump(proof,handle,indent=2);handle.write('\n')
    print(json.dumps(proof),flush=True)


if __name__=='__main__':main()
