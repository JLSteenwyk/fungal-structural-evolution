#!/usr/bin/env python3
"""Verify every full-group alignment independently of the producer parser."""
import argparse,csv,hashlib,json,re,time
from pathlib import Path
from collections import Counter
import psutil

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def fasta(path):
    records={};key=None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            key=line[1:].split()[0]
            if key in records:raise ValueError('Duplicate FASTA identity')
            records[key]=''
        elif line.strip():
            if key is None:raise ValueError('Sequence before header')
            records[key]+=line.strip().upper()
    return records

def check_sequences(source,aligned):
    assert source and list(source)==list(aligned)
    lengths={len(s) for s in aligned.values()};assert len(lengths)==1 and next(iter(lengths))>0
    assert all(aligned[t].replace('-','')==s for t,s in source.items())
    assert all(set(col)!={'-'} for col in zip(*aligned.values()))
    return len(source),next(iter(lengths)),sum(map(len,source.values()))

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();p=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for path,h in p['pins'].items():assert sha(path)==h,path
    verify();dep=p['producer']
    while True:
        try:
            process=psutil.Process(dep['pid'])
            if abs(process.create_time()-dep['created'])>.01 or process.status()==psutil.STATUS_ZOMBIE:break
            assert process.cmdline()==dep['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    verify();sp=json.loads(Path(p['source_plan']).read_text());root=Path(sp['output']);rp=root/'receipt.json';r=json.loads(rp.read_text());rh=sha(rp);inputs=Path(sp['inputs']);ir=json.loads((inputs/'receipt.json').read_text())
    assert r['status']=='complete_full_codon_group_protein_realignments_pending_readback' and r['plan_sha256']==sha(p['source_plan'])
    assert r['source_receipt_sha256']==sha(inputs/'receipt.json') and r['source_readback_sha256']==sha(sp['readback'])
    for name,h in ir['artifacts'].items():assert sha(inputs/name)==h
    with (inputs/'cases.tsv').open() as f:cases={x['case_id']:x for x in csv.DictReader(f,delimiter='\t')}
    receipts={x['case_id']:x['receipt_sha256'] for x in r['case_receipts']};assert len(receipts)==len(r['case_receipts'])==len(cases)==r['cases']==1712 and set(receipts)==set(cases)
    rows=[];strategies=Counter()
    for cid,case in sorted(cases.items()):
        folder=root/cid;cp=folder/'receipt.json';assert sha(cp)==receipts[cid];cr=json.loads(cp.read_text());source=inputs/cid/'full_amino_acids.faa'
        assert cr['case_id']==cid and cr['status']=='complete_full_group_protein_realignment_pending_readback' and cr['plan_sha256']==sha(p['source_plan'])
        assert cr['input_sha256']==sha(source)
        assert cr['command']==[sp['mafft'],'--auto','--thread','1','--inputorder',str(source.resolve())]
        for key in ['translation_table','marker_copy_caveat']:assert cr[key]==case[key]
        assert set(cr['artifacts'])=={'aligned_amino_acids.faa','mafft.log'}
        for name,h in cr['artifacts'].items():assert sha(folder/name)==h
        taxa,columns,residues=check_sequences(fasta(source),fasta(folder/'aligned_amino_acids.faa'))
        assert taxa==cr['taxa']==int(case['taxa']) and columns==cr['columns'] and residues==int(case['full_amino_acids'])
        strategy=re.findall(r'(?m)^Strategy:\s*\n\s*([^\n]+)',(folder/'mafft.log').read_text());assert len(strategy)==1
        strategies[strategy[0]]+=1
        rows.append(dict(case_id=cid,taxa=taxa,alignment_columns=columns,full_amino_acids=residues,translation_table=case['translation_table'],marker_copy_caveat=case['marker_copy_caveat'],strategy=strategy[0],case_receipt_sha256=receipts[cid]))
    verify();assert sha(rp)==rh
    out=Path(p['output']);out.mkdir(parents=True,exist_ok=False)
    with (out/'cases.tsv').open('w') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    proof=dict(status='passed_full_codon_group_protein_alignment_readback',plan_sha256=ph,producer_receipt_sha256=rh,cases=len(rows),case_taxon_sequences=sum(x['taxa'] for x in rows),full_amino_acids=sum(x['full_amino_acids'] for x in rows),strategy_counts=dict(strategies),artifacts={'cases.tsv':sha(out/'cases.tsv')},scope='Independent manual FASTA parser verifies every source taxon/order and ungapped sequence, rectangular alignment, no all-gap column, exact receipt/command/code/copy fields and all log/output hashes. Strategy text recorded without asserting algorithmic optimality. No codon projection, homology correctness, selection eligibility or biological effect claim.')
    (out/'receipt.json').write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof,indent=2))
if __name__=='__main__':main()
