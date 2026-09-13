#!/usr/bin/env python3
"""Split verified domain alignment into repeat-specific matrices and summaries."""
import argparse,csv,json
from collections import Counter
from pathlib import Path
from Bio import SeqIO
from audit_busco_gene_copies import ROOT,sha,read_table
AA=set('ACDEFGHIKLMNPQRSTVWY')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable repeat sensitivity input')
    base=ROOT/'results/phylogeny/tfiib-domain-pairs-v1';r=json.loads((base/'receipt.json').read_text());vp=ROOT/'metadata/tfiib_domain_alignment_readback.json';v=json.loads(vp.read_text())
    if r['status']!='complete_ordered_tfiib_domain_pair_alignment' or v['source_receipt_sha256']!=sha(base/'receipt.json'):raise ValueError('Verified source alignment required')
    for name,h in r['artifacts'].items():
        if sha(base/name)!=h:raise ValueError('Changed domain alignment')
    seqs={x.id:str(x.seq) for x in SeqIO.parse(base/'paired_domains.faa','fasta')};length=r['profile_columns_per_repeat'];mapping={x['entry_id']:x for x in read_table(base/'candidate_pair_audit.tsv')}
    a.output.mkdir(parents=True);summaries=[];rows=[]
    for repeat in [1,2]:
        split={k:s[(repeat-1)*length:repeat*length] for k,s in seqs.items()};target=a.output/('repeat'+str(repeat)+'.faa')
        if {len(s) for s in split.values()}!={length}:raise ValueError('Invalid repeat boundaries')
        with target.open('w') as f:
            for key,seq in sorted(split.items()):f.write('>'+key+'\n'+seq+'\n')
        constants=variable=informative=empty=0
        for col in zip(*split.values()):
            c=Counter(x for x in col if x in AA)
            empty+=not c;constants+=len(c)==1;variable+=len(c)>1;informative+=sum(n>=2 for n in c.values())>=2
        for key,seq in split.items():rows.append({'entry_id':key,'taxon_id':mapping[key]['taxon_id'],'protein_id':mapping[key]['protein_id'],'repeat_order':repeat,'canonical_residues':sum(x in AA for x in seq),'gap_residues':seq.count('-'),'ambiguous_residues':sum(x not in AA and x!='-' for x in seq),'source_BRF1_hits':mapping[key]['BRF1_hits'],'source_gene_status':mapping[key]['gene_mapping_status']})
        summaries.append({'repeat_order':repeat,'proteins':len(split),'columns':length,'distinct_aligned_sequences':len(set(split.values())),'constant_columns':constants,'variable_columns':variable,'parsimony_informative_columns':informative,'all_gap_or_ambiguous_columns':empty,'path':target.name,'sha256':sha(target)})
    with (a.output/'repeat_coverage.tsv').open('w') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    result={'status':'complete_tfiib_repeat_sensitivity_inputs','source_receipt_sha256':sha(base/'receipt.json'),'source_readback_sha256':sha(vp),'repeats':summaries,'script_sha256':sha(Path(__file__)),'interpretation':'Same 1040 gene identities split into the existing 92-column positional-repeat blocks; no additional filtering or new alignment. Descriptive site variation does not establish phylogenetic resolution or independence. Separate-repeat topology/support comparisons remain pending.','artifacts':{p.name:sha(p) for p in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
