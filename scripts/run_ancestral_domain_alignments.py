#!/usr/bin/env python3
"""Run both domain boundaries for all candidate clades through two aligners."""
import argparse,csv,fcntl,json,subprocess,time
from pathlib import Path
from Bio import SeqIO
from prepare_case_ancestral_neighborhoods import sha,read


def load(path):
    records=list(SeqIO.parse(path,'fasta'));data={r.id:str(r.seq) for r in records}
    assert len(records)==len(data)
    return data


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);args=parser.parse_args()
    plan=json.loads(args.plan.read_text());digest=sha(args.plan)
    def verify():
        assert sha(args.plan)==digest
        for path,h in plan['pins'].items():assert sha(path)==h,path
    verify();inputs=Path(plan['inputs']);source=json.loads((inputs/'receipt.json').read_text())
    proof=json.loads(Path(plan['input_proof']).read_text());assert proof['source_receipt_sha256']==sha(inputs/'receipt.json')
    assert source['status']=='complete_case_domain_sequence_inputs_with_residue_readback'
    families=read(inputs/'input_summary.tsv');assert len(families)==26
    assert len({(r['family'],r['boundary']) for r in families})==26
    assert {r['boundary'] for r in families}=={'alignment','envelope'}
    assert sha(inputs/'input_summary.tsv')==source['artifacts']['input_summary.tsv']
    output=Path(plan['output']);output.mkdir(parents=True,exist_ok=True)
    lock=(output/'run.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    cp=output/'plan.json'
    if cp.exists():assert json.loads(cp.read_text())==plan
    else:cp.write_text(json.dumps(plan,indent=2)+'\n')
    results=[]
    for family in families:
        fid=family['family']+'-'+family['boundary'];input_path=inputs/family['fasta'];assert sha(input_path)==source['artifacts'][input_path.name]==family['sha256']
        original=load(input_path);assert len(original)==int(family['proteins'])
        for method in ['mafft','famsa']:
            folder=output/(fid+'-'+method);rp=folder/'receipt.json'
            if rp.exists():
                r=json.loads(rp.read_text());assert r['plan_sha256']==digest and r['status']=='complete_exact_sequence_preserving_alignment'
                for p,h in r['artifacts'].items():assert sha(folder/p)==h
                results.append(r);continue
            # An unfinished attempt is retained rather than silently overwritten.
            folder.mkdir(exist_ok=False);alignment=folder/'alignment.faa'
            if method=='mafft':command=[plan['mafft'],'--amino','--localpair','--maxiterate','1000','--thread','4','--threadit','0',str(input_path.resolve())]
            else:command=[plan['famsa'],'-t','4','-keep-duplicates',str(input_path.resolve()),str(alignment.resolve())]
            started=time.monotonic()
            with (folder/'stderr.log').open('w') as err,(folder/('alignment.faa' if method=='mafft' else 'stdout.log')).open('w') as out:
                subprocess.run(command,stdout=out,stderr=err,check=True)
            aligned=load(alignment);assert set(aligned)==set(original)
            assert len({len(s) for s in aligned.values()})==1
            assert {g:s.replace('-','') for g,s in aligned.items()}==original
            width=len(next(iter(aligned.values())));n=len(aligned)
            columns=[]
            for pos in range(width):
                residues=[s[pos] for s in aligned.values()];nongap=sum(x!='-' for x in residues)
                columns.append(dict(column=pos+1,observed_residues=nongap,canonical_residues=sum(x in 'ACDEFGHIKLMNPQRSTVWY' for x in residues),total_sequences=n))
            with (folder/'column_coverage.tsv').open('w') as f:
                w=csv.DictWriter(f,list(columns[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(columns)
            assert all(x['observed_residues']>0 for x in columns)
            verify()
            r=dict(status='complete_exact_sequence_preserving_alignment',family=family['family'],boundary=family['boundary'],input_set=fid,method=method,plan_sha256=digest,input_sha256=sha(input_path),command=command,proteins=n,columns=width,columns_at_least_50pct=sum(2*x['observed_residues']>=n for x in columns),columns_at_least_70pct=sum(10*x['observed_residues']>=7*n for x in columns),columns_at_least_90pct=sum(10*x['observed_residues']>=9*n for x in columns),elapsed_seconds=time.monotonic()-started,artifacts={p.name:sha(p) for p in folder.iterdir()})
            rp.write_text(json.dumps(r,indent=2)+'\n');results.append(r);print(fid,method,'complete',width,'columns',flush=True)
    assert len(results)==52
    result=dict(status='complete_52_ancestral_domain_alignments_pending_independent_audit',plan_sha256=digest,results=results,scope='All 13 candidate domains under both boundary definitions and both alignment methods; every input residue/copy preserved including partial HMM hits. Column occupancy descriptive only, not homology correctness or ancestral inference.')
    (output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
