#!/usr/bin/env python3
"""Project whole alignments to exact domain members/intervals and compare columns."""
import csv,gzip,json,hashlib
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
from prepare_case_ancestral_neighborhoods import sha,read


def fasta(path):
    rs=list(SeqIO.parse(path,'fasta'));d={r.id:str(r.seq) for r in rs};assert len(d)==len(rs);return d


def main():
    pp=Path('metadata/whole_domain_ancestral_alignment_comparison_plan_20260927.json');plan=json.loads(pp.read_text())
    def verify():
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify();whole=Path(plan['whole']);domain=Path(plan['domain']);seqroot=Path(plan['sequences']);dr=json.loads((domain/'receipt.json').read_text());wr=json.loads((whole/'receipt.json').read_text());sr=json.loads((seqroot/'receipt.json').read_text())
    wjobs={(r['family'],r['method']):r for r in wr['results']};djobs={(r['family'],r['boundary'],r['method']):r for r in dr['results']}
    dispositions=read(seqroot/'protein_dispositions.tsv');intervals=defaultdict(dict)
    for r in dispositions:
        if r['status']=='extracted_single_target_hit':intervals[r['family'],r['boundary']][r['gene']]=(int(r['start']),int(r['end']))
    maps=defaultdict(lambda:defaultdict(list))
    with gzip.open(plan['coordinates'],'rt') as h:
        for r in csv.DictReader(h,delimiter='\t'):maps[r['family'],r['boundary'],r['method']][int(r['column'])].append((r['gene'],int(r['protein_position'])))
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);summary=[];projections=[];nrows=0
    fields=['family','boundary','whole_method','domain_method','whole_column','domain_column','status','retained_residues','coordinate_signature_sha256']
    with gzip.open(out/'column_correspondence.tsv.gz','wt') as h:
        writer=csv.DictWriter(h,fields,delimiter='\t',lineterminator='\n');writer.writeheader()
        for (family,boundary),bounds in sorted(intervals.items()):
            f=seqroot/(family+'-'+boundary+'.faa');assert sha(f)==sr['artifacts'][f.name];extracted=fasta(f);assert set(extracted)==set(bounds)
            for wm in ['mafft','famsa']:
                path=whole/(family+'-'+wm)/'alignment.faa';assert sha(path)==wjobs[family,wm]['artifacts']['alignment.faa'];aligned=fasta(path);assert set(bounds)<=set(aligned)
                cols=defaultdict(list);residue_count=0
                for gene,(start,end) in bounds.items():
                    pos=0;seen=[]
                    for col,aa in enumerate(aligned[gene],1):
                        if aa=='-':continue
                        pos+=1
                        if start<=pos<=end:cols[col].append((gene,pos));seen.append(aa);residue_count+=1
                    assert ''.join(seen)==extracted[gene] and len(seen)==end-start+1
                projected={tuple(sorted(v)):col for col,v in cols.items()};assert len(projected)==len(cols)
                assert residue_count==sum(map(len,extracted.values()))
                projections.append(dict(family=family,boundary=boundary,method=wm,whole_proteins=len(aligned),domain_proteins=len(bounds),excluded_no_domain_proteins=len(aligned)-len(bounds),whole_columns=len(next(iter(aligned.values()))),projected_columns=len(projected),columns_without_retained_domain_residues=len(next(iter(aligned.values())))-len(projected),retained_residues=residue_count))
                for dm in ['mafft','famsa']:
                    coords=maps[family,boundary,dm];dmap={tuple(sorted(v)):col for col,v in coords.items()};assert len(dmap)==len(coords)==djobs[family,boundary,dm]['columns']
                    assert {item for sig in dmap for item in sig}=={item for sig in projected for item in sig}
                    assert sum(map(len,dmap))==residue_count
                    common=set(projected)&set(dmap);union=set(projected)|set(dmap)
                    for sig in sorted(union):
                        a,b=projected.get(sig),dmap.get(sig);status='matched_projected_column' if a and b else ('only_whole_projection' if a else 'only_domain_alignment')
                        writer.writerow(dict(family=family,boundary=boundary,whole_method=wm,domain_method=dm,whole_column=a or '',domain_column=b or '',status=status,retained_residues=len(sig),coordinate_signature_sha256=hashlib.sha256(json.dumps(sig,separators=(',',':')).encode()).hexdigest()));nrows+=1
                    summary.append(dict(family=family,boundary=boundary,whole_method=wm,domain_method=dm,proteins=len(bounds),whole_projected_columns=len(projected),domain_columns=len(dmap),exact_shared_columns=len(common),only_whole_projection=len(projected)-len(common),only_domain_alignment=len(dmap)-len(common),shared_fraction_whole=len(common)/len(projected),shared_fraction_domain=len(common)/len(dmap)))
    assert len(summary)==104 and len(projections)==52
    for filename,rows in [('comparison_summary.tsv',summary),('projection_summary.tsv',projections)]:
        with (out/filename).open('w') as h:w=csv.DictWriter(h,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    # Complete independent saved-record count/identity and projection checks.
    counts=defaultdict(lambda:defaultdict(int));seen=set()
    with gzip.open(out/'column_correspondence.tsv.gz','rt') as h:
        for r in csv.DictReader(h,delimiter='\t'):
            key=tuple(r[k] for k in ['family','boundary','whole_method','domain_method']);identity=key+(r['coordinate_signature_sha256'],);assert identity not in seen;seen.add(identity);counts[key][r['status']]+=1
            assert bool(r['whole_column'])==(r['status']!='only_domain_alignment') and bool(r['domain_column'])==(r['status']!='only_whole_projection')
    assert len(seen)==nrows
    for r in summary:
        c=counts[tuple(r[k] for k in ['family','boundary','whole_method','domain_method'])]
        assert [c[k] for k in ['matched_projected_column','only_whole_projection','only_domain_alignment']]==[r[k] for k in ['exact_shared_columns','only_whole_projection','only_domain_alignment']]
    verify();r=dict(status='complete_whole_domain_projected_alignment_correspondence',comparisons=104,projections=52,column_union_records=nrows,comparisons_with_identical_projected_columns=sum(x['only_whole_projection']==x['only_domain_alignment']==0 for x in summary),plan_sha256=sha(pp),artifacts={p.name:sha(p) for p in out.iterdir()},scope='Exact coordinate projection to shared annotated proteins and domain intervals; all matched/unmatched column signatures retained. Projection does not realign whole proteins or remove influence of excluded proteins/flanking sequence on the original alignment. Not homology probabilities or isolation of context versus membership effects.')
    (out/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)

if __name__=='__main__':main()
