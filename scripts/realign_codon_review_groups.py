#!/usr/bin/env python3
"""Realign flagged groups and compare residue correspondences with global masks."""
import argparse,csv,itertools,json,shutil,statistics,subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from audit_busco_gene_copies import ROOT,sha,read_table
from align_phylogenetic_markers import validate_alignment
AA=set('ACDEFGHIKLMNPQRSTVWY')


def residue_pairs(a,b,columns,valid_a,valid_b):
    if len(a)!=len(b):raise ValueError('Ragged alignment')
    chosen=set(columns);i=j=0;pairs=set();different=0
    for col,(x,y) in enumerate(zip(a,b),1):
        i+=x!='-';j+=y!='-'
        if col in chosen and x in AA and y in AA and i in valid_a and j in valid_b:
            pairs.add((i,j));different+=x!=y
    return pairs,different


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable realignment review')
    selection=ROOT/'metadata/codon_divergence_alignment_review.tsv';review=read_table(selection)
    coverage=ROOT/'results/cds/genus-code-coverage-v1';gr=json.loads((coverage/'receipt.json').read_text())
    for name,h in gr['artifacts'].items():
        if sha(coverage/name)!=h:raise ValueError('Changed coverage source')
    groups={(x['marker'],x['genus_label'],x['translation_table'],x['policy']):x for x in read_table(coverage/'group_codon_coverage.tsv')}
    base=ROOT/'results/cds/marker-codon-alignments-v1';cr=json.loads((base/'receipt.json').read_text())
    if sha(base/'receipt.json')!=gr['source_codon_receipt_sha256']:raise ValueError('Changed codon provenance')
    masks=json.loads((base/'protein_source_columns.json').read_text())
    if sha(base/'protein_source_columns.json')!=cr['artifacts']['protein_source_columns.json']:raise ValueError('Changed source masks')
    index=ROOT/'results/cds/marker-source-index-v1';ix=json.loads((index/'receipt.json').read_text())
    if sha(index/'receipt.json')!=cr['source_index_receipt_sha256']:raise ValueError('Changed source CDS index')
    for name,h in ix['artifacts'].items():
        if sha(index/name)!=h:raise ValueError('Changed indexed source')
    dna={r.id:str(r.seq) for r in SeqIO.parse(index/'marker_source_cds.fna','fasta')}
    arpath=ROOT/'results/phylogeny/mafft-audit-full-v1/receipt.json'
    if sha(arpath)!=cr['protein_alignment_audit_receipt_sha256']:raise ValueError('Changed global protein audit')
    alignment_hashes={r['marker']:r['alignment_sha256'] for r in json.loads(arpath.read_text())['alignment_receipts']}
    mafft=shutil.which('mafft');version=subprocess.run([mafft,'--version'],capture_output=True,text=True,check=True)
    version=(version.stdout+version.stderr).strip();a.output.mkdir(parents=True);summaries=[];pairrows=[];runs=[]
    for n,case in enumerate(review,1):
        marker=case['marker'];genus=case['genus_label'];code=int(case['translation_table']);g=groups[marker,genus,str(code),case['policy']];taxa=sorted(json.loads(g['retained_taxa_json']))
        path=ROOT/'results/phylogeny/alignments-full-v1'/(marker+'.faa')
        if sha(path)!=alignment_hashes[marker]:raise ValueError('Changed global alignment')
        global_aa={r.id:str(r.seq).upper() for r in SeqIO.parse(path,'fasta') if r.id in taxa};valid={}
        source=a.output/(marker+'.'+genus+'.input.faa')
        with source.open('w') as f:
            for t in taxa:
                aa=global_aa[t].replace('-','');coding=dna[marker+'|'+t];translated=str(Seq(coding).translate(table=code))
                if translated.endswith('*'):translated=translated[:-1];coding=coding[:-3]
                if aa!=translated:raise ValueError('Full source translation differs')
                valid[t]={i//3+1 for i in range(0,len(coding),3) if set(coding[i:i+3])<=set('ACGT')}
                f.write('>'+t+'\n'+aa+'\n')
        target=source.with_name(source.name.replace('.input.faa','.aligned.faa'));log=target.with_suffix('.log')
        command=[mafft,'--auto','--thread','1','--inputorder',str(source.resolve())]
        with target.open('w') as f,log.open('w') as err:subprocess.run(command,stdout=f,stderr=err,check=True)
        validate_alignment(source,target);local={r.id:str(r.seq).upper() for r in SeqIO.parse(target,'fasta')}
        columns=[i for i,col in enumerate(zip(*(local[t] for t in taxa)),1) if sum(x in AA for x in col)/len(taxa)>=.8]
        current=[]
        for t,u in itertools.combinations(taxa,2):
            old,dold=residue_pairs(global_aa[t],global_aa[u],masks[marker],valid[t],valid[u])
            new,dnew=residue_pairs(local[t],local[u],columns,valid[t],valid[u])
            shared=old&new;union=old|new
            row={'marker':marker,'genus_label':genus,'translation_table':code,'taxon_a':t,'taxon_b':u,'global_shared_codons':len(old),'local_shared_codons':len(new),'global_aa_difference_fraction':dold/len(old) if old else '', 'local_aa_difference_fraction':dnew/len(new) if new else '', 'retained_original_residue_pairs':len(shared),'residue_pair_jaccard':len(shared)/len(union) if union else '', 'global_pair_fraction_retained':len(shared)/len(old) if old else ''}
            current.append(row);pairrows.append(row)
        baseline=statistics.median(x['global_aa_difference_fraction'] for x in current if x['global_shared_codons']>=100)
        if abs(baseline-float(case['median_amino_acid_difference_fraction']))>1e-12:raise ValueError('Independent global baseline does not reproduce codon pair screen')
        eligible=[x for x in current if x['local_shared_codons']>=100]
        summaries.append({'marker':marker,'genus_label':genus,'translation_table':code,'taxa':len(taxa),'global_mask_columns':len(masks[marker]),'local_80pct_columns':len(columns),'global_median_aa_difference':baseline,'local_median_aa_difference':statistics.median(x['local_aa_difference_fraction'] for x in eligible) if eligible else '', 'local_pairs_at_least_100_codons':len(eligible),'median_original_pair_fraction_retained':statistics.median(x['global_pair_fraction_retained'] for x in current if x['global_shared_codons']),'median_residue_pair_jaccard':statistics.median(x['residue_pair_jaccard'] for x in current if x['residue_pair_jaccard']!='')})
        runs.append({'marker':marker,'genus_label':genus,'command':command,'input_sha256':sha(source),'output_sha256':sha(target),'local_columns':columns});print(n,marker,genus,summaries[-1],flush=True)
    for name,data in [('realignment_summary.tsv',summaries),('pair_correspondence.tsv',pairrows)]:
        with (a.output/name).open('w') as f:
            w=csv.DictWriter(f,list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    result={'status':'complete_targeted_codon_alignment_review','cases':len(summaries),'pairs':len(pairrows),'selection_sha256':sha(selection),'coverage_receipt_sha256':sha(coverage/'receipt.json'),'codon_receipt_sha256':sha(base/'receipt.json'),'mafft_version':version,'mafft_sha256':sha(Path(mafft)),'runs':runs,'script_sha256':sha(Path(__file__)),'interpretation':'Targeted review of the 20 largest observed marker/group median AA differences, not a representative calibration sample or a separate pilot. Exact same full source proteins realigned within the recorded group using MAFFT auto; algorithm choice can change with group size. Local 80-percent occupancy differs from original full-dataset 50-percent mask. Differences can reflect taxon sampling, alignment and site selection. Original residue-pair retention/Jaccard measures correspondence sensitivity; lower divergence does not establish alignment correctness, orthology or absence of biological change. Canonical source codons define comparable residues; original global divergence independently reproduced.','artifacts':{p.name:sha(p) for p in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
