#!/usr/bin/env python3
"""Run and independently read back SIC indel coding for all candidate alignments."""
import csv,json,re,subprocess
from pathlib import Path
from Bio import SeqIO
from prepare_case_ancestral_neighborhoods import sha


def fasta(p):
    rs=list(SeqIO.parse(p,'fasta'));d={r.id:str(r.seq) for r in rs};assert len(d)==len(rs);return d


def expected(sequences,omit):
    gaps={};unknown={};characters=set()
    for gene,seq in sequences.items():
        gs=[(m.start(),m.end()) for m in re.finditer('-+',seq)]
        terminal=[x for x in gs if x[0]==0 or x[1]==len(seq)] if omit else []
        gaps[gene]=[x for x in gs if x not in terminal];unknown[gene]=terminal+[(m.start(),m.end()) for m in re.finditer('X+',seq)];characters.update(gaps[gene])
    chars=sorted(characters);result={}
    for gene in sequences:
        states=[]
        for a,b in chars:
            state='1' if (a,b) in gaps[gene] else ('?' if any(c<=a and d>=b for c,d in gaps[gene]) else '0')
            if any(c<=a and d>=b for c,d in unknown[gene]) or (state=='1' and any(d==a or c==b for c,d in unknown[gene])):state='?'
            states.append(state)
        result[gene]=''.join(states)
    return chars,result


def encode(binary,input_path,folder,omit):
    folder.mkdir(parents=True,exist_ok=False);folder=folder.resolve();input_path=Path(input_path).resolve()
    settings=dict(_seqFile=str(input_path),_indelOutputInfoFile=str(folder/'characters.txt'),_indelOutputFastaFile=str(folder/'characters.faa'),_nexusFileName=str(folder/'characters.nex'),_logFile=str(folder/'coder.log'),_logValue=4,_codingType='SIC',_isOmitLeadingAndEndingGaps=int(omit))
    params=folder/'parameters.txt';params.write_text(''.join(f'{k} {v}\n' for k,v in settings.items()))
    with (folder/'stdout.log').open('w') as h:subprocess.run([binary,str(params)],cwd=folder,stdout=h,stderr=subprocess.STDOUT,check=True)
    seq=fasta(input_path);chars,reference=expected(seq,omit);observed=fasta(folder/'characters.faa');assert observed==reference,(input_path,omit,'matrix mismatch')
    text=(folder/'characters.txt').read_text();parsed=[tuple(map(int,x)) for x in re.findall(r'character number: (\d+)\s+Start position relative to MSA: (\d+)\s+End position relative to MSA: (\d+)\s+Length: (\d+)',text)]
    assert [(i,a,b,b-a) for i,(a,b) in enumerate(chars)]==parsed,(input_path,omit,'coordinates mismatch')
    return seq,chars,observed


def main():
    pp=Path('metadata/ancestral_indel_coding_plan_20260927.json');plan=json.loads(pp.read_text())
    def verify():
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify();out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    toy=out/'analytic.faa';toy.write_text('>short\nAAAA--AAAAAA\n>long\nAAA----AAAAA\n>partial\nAAAAA--AAAAA\n>terminal\n--AAAAAAAA--\n>unknown\nAAAAXXAAAAAA\n>flanked\nAAAA--XAAAAA\n>intact\nAAAAAAAAAAAA\n')
    for omit in [False,True]:
        seq,chars,states=encode(plan['binary'],toy,out/('analytic-'+str(int(omit))),omit);j=chars.index((4,6));assert states['short'][j]=='1' and states['long'][j]=='?' and states['partial'][j]=='0' and states['unknown'][j]=='?' and states['flanked'][j]=='?' and states['intact'][j]=='0'
        assert ((0,2) in chars)==(not omit) and ((10,12) in chars)==(not omit)
    summaries=[];coordinates=[]
    for job in plan['jobs']:
        assert sha(job['alignment'])==job['alignment_sha256']
        for omit in [False,True]:
            policy='terminal_unknown' if omit else 'terminal_gap';name=job['input_id']+'-'+policy;folder=out/name
            seq,chars,states=encode(plan['binary'],job['alignment'],folder,omit);assert len(seq)==job['proteins']
            flat=''.join(states.values());r=dict(input_id=job['input_id'],family=job['family'],boundary=job['boundary'],method=job['method'],terminal_policy=policy,proteins=len(seq),alignment_columns=job['columns'],characters=len(chars),coded_cells=len(flat),gap_present_cells=flat.count('1'),gap_absent_cells=flat.count('0'),unknown_cells=flat.count('?'),source_alignment_sha256=job['alignment_sha256'])
            for i,(a,b) in enumerate(chars):coordinates.append(dict(input_id=job['input_id'],terminal_policy=policy,character=i+1,start_column=a+1,end_column=b,length=b-a))
            receipt=dict(status='complete_SIC_encoding_full_matrix_readback',job=job,summary=r,plan_sha256=sha(pp),artifacts={f.name:sha(f) for f in folder.iterdir()});(folder/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');summaries.append(r);print(name,len(seq),len(chars),'all_cells_verified',flush=True)
    assert len(summaries)==156
    for name,rows in [('coding_summary.tsv',summaries),('character_coordinates.tsv',coordinates)]:
        with (out/name).open('w') as h:w=csv.DictWriter(h,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    verify();r=dict(status='complete_156_indel_codings_with_full_matrix_and_coordinate_readback',encodings=156,maximum_proteins=max(r['proteins'] for r in summaries),character_records=len(coordinates),coded_cells=sum(r['coded_cells'] for r in summaries),analytic_nested_partial_terminal_unknown_checks='passed',plan_sha256=sha(pp),artifacts={f.name:sha(f) for f in out.iterdir() if f.is_file()},job_receipts={str(f.relative_to(out)):sha(f) for f in out.glob('*/receipt.json')},scope='SIC gap-character matrices, not ancestral indel inference. 1=exact coded gap,0=no exact gap under SIC,?=containing gap/qualifying unknown. Partial overlaps can be0; not equivalent to residue presence. Terminal policies both retained. Multi-position characters can remain dependent; no event count or ancestral probability claim.')
    (out/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k not in ['artifacts','job_receipts']}),flush=True)

if __name__=='__main__':main()
