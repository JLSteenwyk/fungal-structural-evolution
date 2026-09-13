#!/usr/bin/env python3
"""Audit raw symmetry results and independently reconstruct Bowker pair counts."""
import argparse,csv,io,json,math
from itertools import combinations
from pathlib import Path
import numpy as np
from scipy.stats import chi2
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha,read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['screen','inputs','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use immutable output')
    receipt=checked_receipt(a.screen);inputs=checked_receipt(a.inputs);config=json.loads((a.screen/'config.json').read_text())
    if receipt['config_sha256']!=sha(a.screen/'config.json') or config['input_receipt_sha256']!=sha(a.inputs/'receipt.json'):raise ValueError('Source lineage differs')
    rows=read_table(a.screen/'raw_symmetry_summary.tsv');ready={r['marker']:r for r in read_table(a.inputs/'marker_summary.tsv') if r['status']=='ready_for_inference'}
    if {(r['marker'],r['alphabet']) for r in rows}!={(m,label) for m in ready for label in ['aa','3di']} or len(rows)!=2*len(ready):raise ValueError('Screen universe differs')
    codes={x:i for i,x in enumerate('ACDEFGHIKLMNPQRSTVWY')};out=[];totalpairs=0
    for row in rows:
        marker,label=row['marker'],row['alphabet'];folder=a.screen/marker;rp=folder/(label+'.receipt.json');r=json.loads(rp.read_text())
        if sha(rp)!=receipt['entry_receipts'][marker+'/'+label] or r['config_sha256']!=receipt['config_sha256']:raise ValueError('Entry provenance differs')
        for name,h in r['artifacts'].items():
            if sha(folder/name)!=h:raise ValueError('Entry artifact changed')
        source=a.inputs/marker/(label+'.faa')
        if sha(source)!=r['alignment_sha256']:raise ValueError('Alignment changed')
        raw=list(csv.DictReader(io.StringIO('\n'.join(x for x in (folder/(label+'.symtest.csv')).read_text().splitlines() if x.strip() and not x.startswith('#')))))
        if len(raw)!=1 or any(row[k]!=v for k,v in raw[0].items()) or raw[0]!=r['raw_result']:raise ValueError('Raw summary differs')
        seqs=list(SeqIO.parse(source,'fasta'));n=len(seqs);possible=n*(n-1)//2
        if n!=int(row['taxa']) or any(len(x.seq)!=int(row['columns']) for x in seqs):raise ValueError('Dimensions differ')
        encoded=np.array([[codes[c] if c!='?' else -1 for c in str(x.seq)] for x in seqs],dtype=np.int16)
        sig=non=zero=0
        for i,j in combinations(range(n),2):
            x,y=encoded[i],encoded[j];keep=(x>=0)&(y>=0);matrix=np.bincount(x[keep]*20+y[keep],minlength=400).reshape(20,20)
            sums=matrix+matrix.T;valid=np.triu(sums>0,1);df=int(valid.sum())
            if df==0:zero+=1;continue
            differences=matrix-matrix.T;stat=float(np.sum(differences[valid]**2/sums[valid]));pvalue=float(chi2.sf(stat,df))
            if pvalue<config['nominal_pair_threshold']:sig+=1
            else:non+=1
        if sig!=int(row['SymSig']) or non!=int(row['SymNon']) or sig+non+zero!=possible:raise ValueError('Independent Bowker counts differ: '+marker+' '+label+' '+str((sig,non,zero)))
        for test in ['Sym','Mar','Int']:
            ns,nn=int(row[test+'Sig']),int(row[test+'Non']);value=float(row[test+'Pval'])
            if min(ns,nn)<0 or ns+nn>possible or (math.isfinite(value) and not 0<=value<=1):raise ValueError('Invalid diagnostic values')
            out.append({'marker':marker,'alphabet':label,'test':test,'possible_pairs':possible,'reported_tested_pairs':ns+nn,'pairs_not_in_reported_counts':possible-ns-nn,'reported_maximum_test_pvalue':row[test+'Pval'],'maximum_test_estimable':math.isfinite(value),'nominal_maximum_p_below_0_05':math.isfinite(value) and value<.05})
        totalpairs+=possible;print('checked',marker,label,flush=True)
    a.output.mkdir(parents=True);write_table(a.output/'diagnostic_availability.tsv',out)
    r={'status':'passed_raw_symmetry_and_bowker_count_audit','screen_receipt_sha256':sha(a.screen/'receipt.json'),'script_sha256':sha(Path(__file__)),'alignments':len(rows),'taxon_pairs_independently_reconstructed':totalpairs,'diagnostic_rows':len(out),'nonfinite_maximum_test_results':sum(not x['maximum_test_estimable'] for x in out),'scope':'All raw diagnostics and source/output identities checked. Every pairwise Bowker statistic reconstructed and nominal significant/nonsignificant counts matched. Maximum-test p-values and marginal/internal test statistics were not independently recalculated. Sparse/dependent sites can invalidate nominal calibration; no multiple-testing correction or biological rejection decisions performed.','artifacts':{'diagnostic_availability.tsv':sha(a.output/'diagnostic_availability.tsv')}}
    (a.output/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))

if __name__=='__main__':main()
