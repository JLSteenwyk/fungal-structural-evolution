#!/usr/bin/env python3
"""Reconstruct all projected codons and encoded residue-pair intersections."""
import argparse,csv,json,hashlib,itertools,math
from pathlib import Path
from collections import Counter
import numpy as np
from Bio.Data import CodonTable

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def fasta(p):
    result={};key=None
    for line in Path(p).read_text().splitlines():
        if line.startswith('>'):
            key=line[1:].split()[0];assert key not in result;result[key]=''
        elif line.strip():
            assert key is not None;result[key]+=line.strip().upper()
    return result

def encoded(left,right,base):
    a=np.asarray(left,dtype=np.int64);b=np.asarray(right,dtype=np.int64);valid=(a>0)&(b>0)
    return np.unique(a[valid]*base+b[valid])
def metrics(old,raw,new):
    assert np.isin(new,raw).all()
    sr=len(np.intersect1d(old,raw));sn=len(np.intersect1d(old,new));union=len(np.union1d(old,new))
    return dict(original_pairs=len(old),raw_local_pairs=len(raw),retained_local_pairs=len(new),original_preserved_raw=sr,original_preserved_retained=sn,original_absent_from_raw=len(old)-sr,original_preserved_but_filtered=sr-sn,local_retained_not_original=len(new)-sn,original_fraction_retained=sn/len(old) if len(old) else '',retained_pair_jaccard=sn/union if union else '')

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();p=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for path,h in p['pins'].items():assert sha(path)==h,path
    verify();sp=json.loads(Path(p['projection_plan']).read_text());root=Path(sp['output']);inputs=Path(sp['inputs']);alignments=Path(sp['alignments']);rp=root/'receipt.json';r=json.loads(rp.read_text());rh=sha(rp)
    assert r['plan_sha256']==sha(p['projection_plan']) and r['source_alignment_receipt_sha256']==sha(alignments/'receipt.json') and r['source_alignment_readback_sha256']==sha(sp['alignment_readback'])
    for base in [inputs,root]:
        for name,h in json.loads((base/'receipt.json').read_text())['artifacts'].items():assert sha(base/name)==h
    ar=json.loads((alignments/'receipt.json').read_text());casehash={x['case_id']:x['receipt_sha256'] for x in ar['case_receipts']}
    def rows(path):
        with open(path) as f:return list(csv.DictReader(f,delimiter='\t'))
    sources={x['case_id']:x for x in rows(inputs/'cases.tsv')};summaries=rows(root/'cases.tsv');assert len(summaries)==len(sources)==r['cases']==1712 and {x['case_id'] for x in summaries}==set(sources)
    pairrows=rows(root/'pair_correspondence.tsv');pairs={(x['case_id'],x['taxon_a'],x['taxon_b']):x for x in pairrows};assert len(pairs)==len(pairrows)==r['pairs'];checked=observed=seqs=0;dispositions=Counter()
    for summary in summaries:
        cid=summary['case_id'];source=sources[cid];folder=root/cid;cp=alignments/cid/'receipt.json';assert sha(cp)==casehash[cid]
        for name,h in json.loads(cp.read_text())['artifacts'].items():assert sha(alignments/cid/name)==h
        aa=fasta(alignments/cid/'aligned_amino_acids.faa');dna=fasta(inputs/cid/'full_codons.fna');outdna=fasta(folder/'codons.fna');outaa=fasta(folder/'amino_acids.faa');old=json.loads((inputs/cid/'original_residue_positions.json').read_text());stored=json.loads((folder/'raw_called_residue_positions.json').read_text());taxa=list(aa);assert set(taxa)==set(dna)==set(outdna)==set(outaa)==set(old)==set(stored)
        width=len(aa[taxa[0]]);pos=np.zeros((len(taxa),width),dtype=np.int64);triplets={};table=CodonTable.unambiguous_dna_by_id[int(source['translation_table'])]
        for k,t in enumerate(taxa):
            chars=np.array(list(aa[t]));ranks=np.cumsum(chars!='-');codons=[dna[t][i:i+3] for i in range(0,len(dna[t]),3)];assert len(codons)==int(ranks[-1]);projected=[]
            for j,x in enumerate(chars):
                codon=codons[int(ranks[j])-1] if x!='-' else '???';valid=x in 'ACDEFGHIKLMNPQRSTVWY' and set(codon)<=set('ACGT')
                if x!='-' and set(codon)<=set('ACGT'):assert table.forward_table[codon]==x
                pos[k,j]=int(ranks[j]) if valid else 0;projected.append(codon if valid else '???')
            assert stored[t]==pos[k].tolist();triplets[t]=projected
        called=(pos>0).sum(axis=0);keep=np.flatnonzero(called*5>=4*len(taxa));cols=rows(folder/'columns.tsv')
        assert len(cols)==len(keep)
        for j,(col,i) in enumerate(zip(cols,keep),1):assert col==dict(retained_codon_column_1based=str(j),local_alignment_column_1based=str(i+1),called_taxa=str(called[i]))
        for k,t in enumerate(taxa):
            assert outdna[t]==''.join(triplets[t][i] for i in keep)
            assert outaa[t]==''.join(aa[t][i] if pos[k,i] else '?' for i in keep)
        counts=(pos[:,keep]>0).sum(axis=1);status='passes_existing_coverage_gate' if len(keep)>=100 and int(counts.min())>0 and len(taxa)>=4 else 'below_existing_coverage_gate'
        expected=dict(case_id=cid,translation_table=source['translation_table'],taxa=str(len(taxa)),original_columns=source['original_codon_columns'],local_alignment_columns=str(width),retained_codon_columns=str(len(keep)),observed_taxon_codons=str(int(counts.sum())),minimum_observed_codons=str(int(counts.min())),marker_copy_caveat=source['marker_copy_caveat'],coverage_disposition=status,selection_eligibility='not_established');assert summary==expected
        observed+=int(counts.sum());seqs+=len(taxa);dispositions[status]+=1;base=max(len(x)//3 for x in dna.values())+1;index={t:i for i,t in enumerate(taxa)}
        for t,u in itertools.combinations(sorted(taxa),2):
            row=pairs.pop((cid,t,u));i,j=index[t],index[u];m=metrics(encoded(old[t],old[u],base),encoded(pos[i],pos[j],base),encoded(pos[i,keep],pos[j,keep],base))
            for field,value in m.items():
                if value=='':assert row[field]==''
                elif isinstance(value,int):assert int(row[field])==value
                else:assert math.isclose(float(row[field]),value,rel_tol=1e-12,abs_tol=1e-12)
            checked+=1
    assert not pairs and checked==r['pairs'] and observed==r['observed_taxon_codons'] and dict(dispositions)==r['dispositions']
    verify();assert sha(rp)==rh
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    result=dict(status='passed_full_codon_projection_and_correspondence_readback',plan_sha256=ph,producer_receipt_sha256=rh,cases=len(summaries),case_taxon_sequences=seqs,taxon_pairs=checked,observed_taxon_codons=observed,dispositions=dict(dispositions),scope='Independent FASTA parser, cumulative residue ranks, genetic-code lookup and NumPy encoded-pair intersections reconstruct every codon, amino acid, mask, position map, summary and pair statistic. No producer projection/pair helper imported. Checks numerical provenance, not homology correctness or selection eligibility.')
    Path(p['output']).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
