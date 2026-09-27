#!/usr/bin/env python3
"""Compare qualified structural states at identical taxon/marker alignment positions."""
import argparse,csv,json
from collections import Counter
from pathlib import Path
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from run_ortholog_pair_guide_comparison import sha


def load(root,marker):
    columns=list(csv.DictReader((root/marker/'columns.tsv').open(),delimiter='\t'))
    if [int(r['paired_column_1based']) for r in columns]!=list(range(1,len(columns)+1)):raise ValueError('Invalid paired column grid')
    positions=[int(r['alignment_column_1based']) for r in columns]
    if positions!=sorted(set(positions)):raise ValueError('Ambiguous original positions')
    data={}
    for label in ['aa','3di']:
        records=list(SeqIO.parse(root/marker/(label+'.faa'),'fasta'));data[label]={r.id:str(r.seq) for r in records}
        if len(data[label])!=len(records) or any(len(s)!=len(columns) for s in data[label].values()):raise ValueError('Invalid alignment dimensions')
    if set(data['aa'])!=set(data['3di']):raise ValueError('Different paired taxa')
    for t in data['aa']:
        if [c=='?' for c in data['aa'][t]]!=[c=='?' for c in data['3di'][t]]:raise ValueError('Different paired masks')
    return positions,data


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['reference','local','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();roots=[a.reference,a.local];receipts=[checked_receipt(x) for x in roots]
    if any(r['status']!='complete_paired_phylogenetic_input_preparation' for r in receipts):raise ValueError('Incomplete paired inputs')
    if receipts[0]['source_receipts']['matrix']!=receipts[1]['source_receipts']['matrix']:raise ValueError('Different source marker alignment')
    pins={str(x/'receipt.json'):sha(x/'receipt.json') for x in roots}
    ready=[]
    for root in roots:
        ready.append({r['marker'] for r in csv.DictReader((root/'marker_summary.tsv').open(),delimiter='\t') if r['status']=='ready_for_inference'})
    a.output.mkdir(parents=True,exist_ok=False);cells=[];confusion=Counter();totals=Counter()
    for marker in sorted(ready[0]&ready[1]):
        datasets=[load(root,marker) for root in roots]
        for root,receipt in zip(roots,receipts):
            for name in ['aa.faa','3di.faa','columns.tsv']:
                path=root/marker/name;h=sha(path)
                if h!=receipt['artifacts'][marker+'/'+name]:raise ValueError('Changed paired artifact')
                pins[str(path)]=h
        (pa,da),(pb,db)=datasets;ia={x:i for i,x in enumerate(pa)};ib={x:i for i,x in enumerate(pb)}
        shared=sorted(set(pa)&set(pb))
        for taxon in sorted(set(da['aa'])&set(db['aa'])):
            counts=Counter()
            for pos in shared:
                x,y=ia[pos],ib[pos];aa,ab=da['aa'][taxon][x],db['aa'][taxon][y]
                if aa=='?' or ab=='?':continue
                counts['jointly_observed_positions']+=1
                if aa!=ab:counts['aa_mismatches']+=1;continue
                sa,sb=da['3di'][taxon][x],db['3di'][taxon][y]
                counts['same_aa_positions']+=1;counts['state_mismatches']+=sa!=sb;confusion[sa,sb]+=1
            record=dict(marker=marker,taxon=taxon,shared_alignment_columns=len(shared),reference_observed=sum(c!='?' for c in da['aa'][taxon]),local_observed=sum(c!='?' for c in db['aa'][taxon]))
            for key in ['jointly_observed_positions','aa_mismatches','same_aa_positions','state_mismatches']:record[key]=counts[key]
            record['state_mismatch_fraction']=counts['state_mismatches']/counts['same_aa_positions'] if counts['same_aa_positions'] else ''
            record['status']='compared' if counts['same_aa_positions'] else 'no_joint_same_aa_positions'
            cells.append(record);totals.update(counts)
    def table(path,rows):
        with path.open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    table(a.output/'taxon_marker_state_comparison.tsv',cells)
    table(a.output/'state_confusion.tsv',[dict(reference_state=x,local_state=y,positions=n) for (x,y),n in sorted(confusion.items())])
    assert sum(confusion.values())==totals['same_aa_positions']
    for path,h in pins.items():
        if sha(path)!=h:raise ValueError('Source changed during comparison')
    result=dict(status='complete_paired_source_state_comparison',source_sha256=pins,script_sha256=sha(__file__),shared_taxon_marker_cells=len(cells),taxa=len({r['taxon'] for r in cells}),markers=len({r['marker'] for r in cells}),counts=dict(totals),pooled_state_mismatch_fraction=totals['state_mismatches']/totals['same_aa_positions'] if totals['same_aa_positions'] else None,artifacts={x.name:sha(x) for x in a.output.glob('*.tsv')},scope='Common qualified taxon/marker/original-alignment positions under each source mask. Only identical observed AA positions enter state confusion. This is source-associated disagreement, not a matched full-protein predictor experiment: exact complete protein identity and model context are not established here. Uneven coverage and dependence prevent extrapolation or significance claims; no physical displacement, branch effect or calibrated error rate.')
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['source_sha256','artifacts']},indent=2))

if __name__=='__main__':main()
