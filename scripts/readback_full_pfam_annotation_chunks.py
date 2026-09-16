#!/usr/bin/env python3
"""Independently replay each annotated full-Pfam row against raw HMMER fields."""
import argparse,csv,gzip,itertools,json
from pathlib import Path
from prepare_pfam import digest


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['annotations','plan','output']:ap.add_argument('--'+name,type=Path,required=True)
    a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    r=json.loads((a.annotations/'receipt.json').read_text());p=json.loads(a.plan.read_text())
    if r['plan_sha256']!=digest(a.plan) or r['completed_chunks']!=len(p['chunks']):raise ValueError('Incomplete or mismatched annotation snapshot')
    totals=[]
    for shard in r['shards']:
        key=shard['chunk'];rp=a.annotations/(key+'.annotation.receipt.json');sr=json.loads(rp.read_text());raw=Path(p['search'])/(key+'.domtblout');annotated=a.annotations/(key+'.annotated.tsv.gz')
        if digest(rp)!=shard['receipt_sha256'] or digest(raw)!=sr['raw_table_sha256'] or digest(annotated)!=shard['annotation_sha256']:raise ValueError('Changed shard source/output')
        count=0
        with raw.open() as src,gzip.open(annotated,'rt') as dst:
            source=(line.split(maxsplit=22) for line in src if line.strip() and not line.startswith('#'))
            for f,x in itertools.zip_longest(source,csv.DictReader(dst,delimiter='\t')):
                if f is None or x is None:raise ValueError('Dropped/extra annotated row')
                count+=1
                assert x['hit_id']==f'FULL-{key}-{count:09d}' and x['search_partition']=='additional_full_proteome'
                for name,index in [('sequence_id',0),('pfam_name',3),('pfam_accession',4)]:assert x[name]==f[index]
                for name,index in [('protein_length',2),('hmm_length',5),('hmm_start',15),('hmm_end',16),('alignment_start',17),('alignment_end',18),('envelope_start',19),('envelope_end',20)]:assert int(x[name])==int(f[index])
                for name,index in [('sequence_score',7),('domain_score',13),('conditional_evalue',11),('independent_evalue',12),('posterior_accuracy',21)]:assert float(x[name])==float(f[index])
                assert abs(float(x['hmm_coverage'])-(int(f[16])-int(f[15])+1)/int(f[5]))<1e-12
        assert count==shard['rows']==sr['rows'];totals.append({'chunk':key,'rows_checked':count});print(key,count,flush=True)
    if sum(x['rows_checked'] for x in totals)!=r['annotated_hit_rows']:raise ValueError('Snapshot count differs')
    a.output.mkdir(parents=True)
    result={'status':'passed_full_pfam_annotation_raw_field_readback','chunks_checked':len(totals),'rows_checked':sum(x['rows_checked'] for x in totals),'annotation_receipt_sha256':digest(a.annotations/'receipt.json'),'script_sha256':digest(Path(__file__)),'chunks':totals,'scope':'Every annotated row aligned to raw HMMER order; all source identity, coordinate, score, E-value, posterior-accuracy and HMM-coverage fields independently replayed. Pfam descriptive metadata comes from the pinned parser; not independently parsed here. Does not validate biological assignments or resolve overlaps.','artifacts':{}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
