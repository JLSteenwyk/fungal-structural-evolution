#!/usr/bin/env python3
"""Realign every verified full-protein codon group with immutable per-case receipts."""
import argparse,csv,fcntl,json,hashlib,os,shutil,subprocess,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def validate(source,target):
    left=list(SeqIO.parse(source,'fasta'));right=list(SeqIO.parse(target,'fasta'))
    assert left and len({r.id for r in left})==len(left)
    assert [r.id for r in right]==[r.id for r in left]
    assert len({len(r) for r in right})==1
    for a,b in zip(left,right):assert str(b.seq).upper().replace('-','')==str(a.seq).upper()
    assert all(any(c!='-' for c in col) for col in zip(*(str(r.seq) for r in right)))
    return len(right),len(right[0])
def run_case(job):
    case,p,ph=job;cid=case['case_id'];source=Path(p['inputs'])/cid/'full_amino_acids.faa';out=Path(p['output'])/cid;out.mkdir(parents=True,exist_ok=True);rp=out/'receipt.json';ih=sha(source)
    if rp.exists():
        r=json.loads(rp.read_text());assert r['plan_sha256']==ph and r['input_sha256']==ih and r['case_id']==cid
        for name,h in r['artifacts'].items():assert sha(out/name)==h
        validate(source,out/'aligned_amino_acids.faa');return r
    assert not (out/'aligned_amino_acids.faa').exists(),'Unreceipted output requires review'
    assert shutil.disk_usage(out).free>=p['minimum_free_disk_gib']*1024**3
    command=[p['mafft'],'--auto','--thread','1','--inputorder',str(source.resolve())];start=time.monotonic();tmp=out/'aligned_amino_acids.faa.partial';log=out/'mafft.log'
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    with tmp.open('w') as f,log.open('w') as err:subprocess.run(command,stdout=f,stderr=err,env=env,check=True,timeout=p['per_case_timeout_seconds'])
    taxa,columns=validate(source,tmp);assert taxa==int(case['taxa']) and sha(source)==ih
    target=out/'aligned_amino_acids.faa';tmp.replace(target)
    r=dict(status='complete_full_group_protein_realignment_pending_readback',case_id=cid,translation_table=case['translation_table'],marker_copy_caveat=case['marker_copy_caveat'],taxa=taxa,columns=columns,input_sha256=ih,plan_sha256=ph,command=command,elapsed_seconds=time.monotonic()-start,artifacts={x.name:sha(x) for x in [target,log]})
    rp.write_text(json.dumps(r,indent=2)+'\n');return r

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();p=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for path,h in p['pins'].items():assert sha(path)==h,path
    verify();root=Path(p['inputs']);r=json.loads((root/'receipt.json').read_text());proof=json.loads(Path(p['readback']).read_text());assert proof['producer_receipt_sha256']==sha(root/'receipt.json') and proof['status']=='passed_full_codon_realignment_source_and_position_readback'
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    with (root/'cases.tsv').open() as f:cases=list(csv.DictReader(f,delimiter='\t'))
    assert len(cases)==len({x['case_id'] for x in cases})==r['cases']==1712
    out=Path(p['output']);out.mkdir(parents=True,exist_ok=True);lock=(out/'run.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);completed=[]
    with ThreadPoolExecutor(max_workers=p['workers']) as pool:
        for i,item in enumerate(pool.map(run_case,[(c,p,ph) for c in cases]),1):
            completed.append(dict(case_id=item['case_id'],receipt_sha256=sha(out/item['case_id']/'receipt.json')))
            if i%10==0:
                state=dict(completed_cases=i,total_cases=len(cases));(out/'state.json').write_text(json.dumps(state)+'\n');print(json.dumps(state),flush=True)
    verify();result=dict(status='complete_full_codon_group_protein_realignments_pending_readback',plan_sha256=ph,cases=len(completed),source_receipt_sha256=sha(root/'receipt.json'),source_readback_sha256=sha(p['readback']),case_receipts=completed,scope='All 1712 full-protein groups with exact source identities and order preserved, one MAFFT auto thread per group. Raw local alignments only; codon projection, occupancy filtering, correspondence sensitivity and independent full readback remain downstream. Auto algorithm may differ with group size. No selection eligibility, orthology or alignment correctness claim.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
