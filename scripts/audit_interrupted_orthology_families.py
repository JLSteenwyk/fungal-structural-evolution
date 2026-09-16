#!/usr/bin/env python3
"""Inventory interrupted family artifacts without modifying original files."""
import argparse,csv,hashlib,json,math,time
from collections import Counter
from pathlib import Path
from Bio import SeqIO,Phylo


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--working',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    a.output.mkdir(parents=True,exist_ok=False);start=time.time();counts=Counter();total_genes=0
    fields=['family','source_sequences','source_duplicate_ids','source_sha256','alignment_status','alignment_columns','alignment_sha256','tree_status','tree_sha256','issues']
    with (a.output/'families.tsv').open('w') as out:
        writer=csv.DictWriter(out,fields,delimiter='\t',lineterminator='\n');writer.writeheader()
        for i,path in enumerate(sorted((a.working/'Sequences_ids').glob('*.fa')),1):
            seqs={};duplicates=0
            for rec in SeqIO.parse(path,'fasta'):
                duplicates+=rec.id in seqs;seqs[rec.id]=str(rec.seq).upper()
            row={'family':path.stem,'source_sequences':len(seqs),'source_duplicate_ids':duplicates,'source_sha256':sha(path),'alignment_status':'missing','alignment_columns':'','alignment_sha256':'','tree_status':'missing','tree_sha256':'','issues':''};issues=[]
            if duplicates:issues.append('duplicate_source_ids')
            aln=a.working/'Alignments_ids'/path.name
            if aln.exists():
                seen=set();widths=set();invalid=False;order_bad=0
                try:
                    for rec in SeqIO.parse(aln,'fasta'):
                        text=str(rec.seq).upper();widths.add(len(text));invalid|=rec.id in seen or rec.id not in seqs;seen.add(rec.id)
                        if rec.id in seqs:
                            it=iter(seqs[rec.id]);ungapped=text.replace('-','').replace('.','')
                            if not all(letter in it for letter in ungapped):order_bad+=1
                    invalid|=seen!=set(seqs) or len(widths)!=1 or not seen or 0 in widths or bool(order_bad)
                    row['alignment_status']='content_checks_passed' if not invalid else 'requires_review'
                    row['alignment_columns']=next(iter(widths)) if len(widths)==1 else ''
                    if invalid:issues.append('alignment_grid_or_retained_residue_order')
                except (ValueError,UnicodeError) as e:row['alignment_status']='parse_error';issues.append(type(e).__name__)
                row['alignment_sha256']=sha(aln)
            tree=a.working/'Trees_ids'/(path.stem+'.txt')
            if tree.exists():
                row['tree_sha256']=sha(tree)
                if not tree.stat().st_size:row['tree_status']='empty'
                else:
                    try:
                        obj=Phylo.read(tree,'newick');tips=[n.name for n in obj.get_terminals()]
                        valid=len(tips)==len(set(tips)) and set(tips)==set(seqs) and all(n.branch_length is None or (math.isfinite(n.branch_length) and n.branch_length>=0) for n in obj.find_clades())
                        row['tree_status']='tip_branch_checks_passed' if valid else 'requires_review'
                        if not valid:issues.append('tree_tip_grid_or_branch_lengths')
                    except Exception as e:row['tree_status']='parse_or_traversal_error';issues.append(type(e).__name__)
            row['issues']=';'.join(issues);writer.writerow(row);total_genes+=len(seqs)
            counts['alignment_'+row['alignment_status']]+=1;counts['tree_'+row['tree_status']]+=1;counts['source_duplicate_families']+=bool(duplicates)
            if i%1000==0:out.flush();print('Audited',i,flush=True)
    source_ids={p.stem for p in (a.working/'Sequences_ids').glob('*.fa')}
    orphan={sub:sorted(p.stem for p in (a.working/sub).iterdir() if p.is_file() and p.stem not in source_ids) for sub in ['Alignments_ids','Trees_ids']}
    result={'status':'complete_interrupted_family_artifact_inventory','families':i,'source_sequences':total_genes,'counts':dict(counts),'orphan_artifacts':orphan,'elapsed_seconds':time.time()-start,'script_sha256':sha(Path(__file__)),'artifacts':{'families.tsv':sha(a.output/'families.tsv')},'interpretation':'Artifact inventory, not inferred orthology validation. Missing trees may be expected for small families; disposition requires source-workflow review. Retained-residue subsequence checks allow upstream column trimming; they do not validate alignment homology or reconstruct the exact trimming mask. No files replaced.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
