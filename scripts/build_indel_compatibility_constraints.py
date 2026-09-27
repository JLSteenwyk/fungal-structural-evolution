#!/usr/bin/env python3
"""Build mutually exclusive exact-gap sets and check every extant coding."""
import csv,gzip,json
from collections import defaultdict
from pathlib import Path
import numpy as np
from Bio import SeqIO
from prepare_case_ancestral_neighborhoods import sha


def conflict_sets(intervals):
    # Closed gap endpoints; touching runs would merge into a different exact gap.
    seen=set();result=[]
    for position in sorted({a for a,b in intervals}):
        ids=tuple(i for i,(a,b) in enumerate(intervals) if a<=position<=b+1)
        if len(ids)>1 and ids not in seen:
            assert all(max(intervals[i][0],intervals[j][0])<=min(intervals[i][1]+1,intervals[j][1]+1) for i in ids for j in ids)
            result.append((position,ids));seen.add(ids)
    return result


def main():
    toy=[(1,3),(2,4),(4,6),(8,9)]
    sets=conflict_sets(toy);assert any(ids==(0,1,2) for _,ids in sets)
    assert all(3 not in ids for _,ids in sets)
    # Pairwise probability sums can pass even when a three-way exclusion fails.
    assert .4+.4<=1 and sum([.4,.4,.4])>1
    root=Path('results/ancestral/indel-coding-20260927-v1');rp=root/'receipt.json';receipt=json.loads(rp.read_text());cp=root/'character_coordinates.tsv';assert sha(cp)==receipt['artifacts'][cp.name]
    grouped=defaultdict(list)
    with cp.open() as h:
        for row in csv.DictReader(h,delimiter='\t'):grouped[row['input_id']+'-'+row['terminal_policy']].append(row)
    out=Path('results/ancestral/indel-compatibility-constraints-20260927-v1');out.mkdir(parents=True,exist_ok=False)
    summaries=[];pins={str(rp):sha(rp),str(cp):sha(cp)};checks=0
    with gzip.open(out/'mutually_exclusive_gap_sets.jsonl.gz','wt') as handle:
        for rel,digest in sorted(receipt['job_receipts'].items()):
            source=root/rel;assert sha(source)==digest;r=json.loads(source.read_text());name=source.parent.name
            fasta=source.parent/'characters.faa';assert sha(fasta)==r['artifacts']['characters.faa'];pins[str(source)]=sha(source);pins[str(fasta)]=sha(fasta)
            rows=sorted(grouped[name],key=lambda x:int(x['character']));assert [int(x['character']) for x in rows]==list(range(1,r['summary']['characters']+1))
            intervals=[(int(x['start_column']),int(x['end_column'])) for x in rows]
            sequences=[str(x.seq) for x in SeqIO.parse(fasta,'fasta')]
            observed=np.array([list(s) for s in sequences])=='1'
            sets=conflict_sets(intervals)
            for position,ids in sets:
                assert np.all(observed[:,ids].sum(axis=1)<=1),(name,position)
                checks+=len(sequences)
                handle.write(json.dumps(dict(encoding=name,witness_column=position,characters=[i+1 for i in ids]))+'\n')
            summaries.append(dict(encoding=name,characters=len(intervals),mutually_exclusive_sets=len(sets),maximum_set_size=max((len(ids) for _,ids in sets),default=0),proteins=len(sequences)))
    assert len(summaries)==156
    with (out/'encoding_summary.tsv').open('w') as h:w=csv.DictWriter(h,list(summaries[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(summaries)
    result=dict(status='all_156_exact_gap_exclusion_constraints_checked',encodings=156,mutually_exclusive_sets=sum(r['mutually_exclusive_sets'] for r in summaries),extant_sequence_set_checks=checks,maximum_set_size=max(r['maximum_set_size'] for r in summaries),pins=pins,script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in out.iterdir()},scope='Exact maximal gap runs that overlap or touch are mutually exclusive in a single aligned sequence. Every extant observed1 coding checked; unknown/0 states are not claimed as absence/presence of residues. Any joint exact-gap distribution must have probability sum<=1 in each listed set. Sets may be redundant, not independent events. Posterior feasibility and uncertainty-preserving sequence assembly remain pending.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print({k:v for k,v in result.items() if k not in ['pins','artifacts']})

if __name__=='__main__':main()
