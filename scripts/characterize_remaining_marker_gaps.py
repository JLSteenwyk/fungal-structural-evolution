#!/usr/bin/env python3
"""Describe marker-model gaps and alternative cached sources without launching predictions."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
from Bio import SeqIO


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8388608),b''):h.update(block)
    return h.hexdigest()


def table(path):
    with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))


def write(path,rows,fields=None):
    with path.open('w') as f:
        w=csv.DictWriter(f,fields or list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)


def length_bin(n):
    for upper in [512,768,1024,1536,2048,3072]:
        if n<=upper:return 'at_most_'+str(upper)
    return 'above_3072'


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True)
    a=p.parse_args();plan=json.loads(a.plan.read_text());pins={str(a.plan):sha(a.plan),**plan['pins']}
    def verify():
        for path,digest in pins.items():
            if sha(path)!=digest:raise ValueError('Changed input: '+path)
    verify();out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    availability=table(plan['availability']);receipt=json.loads(Path(plan['availability_receipt']).read_text())
    if sha(plan['availability'])!=receipt['artifacts']['marker_availability.tsv']:raise ValueError('Unbound availability')
    missing=[r for r in availability if r['availability']=='neither_catalog']
    missing_hashes={r['sequence_sha256'] for r in missing}
    by_marker=defaultdict(dict)
    for r in availability:
        if r['taxon_id'] in by_marker[r['marker']]:raise ValueError('Repeated marker/taxon')
        by_marker[r['marker']][r['taxon_id']]=r
    sequences={};lengths={}
    for marker,taxa in by_marker.items():
        rows=list(SeqIO.parse(Path(plan['unaligned'])/(marker+'.faa'),'fasta'))
        if {r.id for r in rows}!=set(taxa) or len(rows)!=len(taxa):raise ValueError('Sequence/taxon universe differs')
        for row in rows:
            seq=str(row.seq);record=taxa[row.id];digest=hashlib.sha256(seq.encode()).hexdigest()
            if digest!=record['sequence_sha256']:raise ValueError('Sequence hash differs')
            lengths[marker,row.id]=len(seq)
            if digest in missing_hashes:sequences[digest]=seq
    if set(sequences)!=missing_hashes:raise ValueError('Incomplete missing-sequence frame')
    latest={};records=0;loghash=hashlib.sha256()
    with Path(plan['inventory']).open('rb') as f:
        for line in f:
            loghash.update(line);row=json.loads(line);records+=1
            models=[m for m in row.get('models',[]) if m['sequence_sha256'] in missing_hashes] if row['status']=='verified' else []
            if models:latest[row['uniprot_accession']]=models
            else:latest.pop(row['uniprot_accession'],None)
    if loghash.hexdigest()!=pins[plan['inventory']]:raise ValueError('Inventory changed during read')
    candidates={};accessions=defaultdict(set)
    for accession,models in latest.items():
        for m in models:
            if m.get('provider')=='GDM' and m.get('tool')=='AlphaFold Monomer v2.0 pipeline':
                raise ValueError('Missing link has an eligible model in the same frozen catalog')
            key=(m['model_id'],str(m['version']))
            fields=['model_id','version','sequence_sha256','length','provider','tool','path','sha256']
            record={k:m.get(k) for k in fields}
            if key in candidates and candidates[key]!=record:raise ValueError('Conflicting candidate identity')
            candidates[key]=record;accessions[key].add(accession)
    by_sequence=Counter()
    candidate_rows=[]
    for key,m in sorted(candidates.items()):
        if sha(m['path'])!=m['sha256'] or len(sequences[m['sequence_sha256']])!=m['length']:
            raise ValueError('Alternative cached model changed or wrong length')
        by_sequence[m['sequence_sha256']]+=1
        candidate_rows.append(dict(m,uniprot_accessions=','.join(sorted(accessions[key]))))
    gaps=[];unique_counts=Counter();link_counts=Counter()
    for r in missing:
        seq=sequences[r['sequence_sha256']];symbols=''.join(sorted(set(seq)-set('ACDEFGHIKLMNPQRSTVWY')))
        gaps.append(dict(r,length=len(seq),length_bin=length_bin(len(seq)),noncanonical_symbols=symbols,
                         alternative_cached_models=by_sequence[r['sequence_sha256']]))
        link_counts[length_bin(len(seq))]+=1
    seq_rows=[]
    for digest,seq in sorted(sequences.items()):
        unique_counts[length_bin(len(seq))]+=1
        seq_rows.append(dict(sequence_sha256=digest,length=len(seq),length_bin=length_bin(len(seq)),
            noncanonical_symbols=''.join(sorted(set(seq)-set('ACDEFGHIKLMNPQRSTVWY'))),alternative_cached_models=by_sequence[digest]))
    manifest=table(plan['manifest']);lineage=defaultdict(Counter)
    taxon_meta={r['taxon_id']:r for r in manifest}
    for r in manifest:
        c=lineage[r['study_role'],r['lineage'].split(';')[0]];c['taxa']+=1;c['marker_slots']+=len(by_marker)
    for r in availability:
        taxon=taxon_meta[r['taxon_id']];c=lineage[taxon['study_role'],taxon['lineage'].split(';')[0]]
        c['recovered_records']+=1;c[r['availability']]+=1
    lineage_rows=[]
    for (role,group),c in sorted(lineage.items()):
        covered=c['both']+c['esmfold_only']+c['alphafold_only']
        lineage_rows.append(dict(study_role=role,manifest_lineage=group,taxa=c['taxa'],marker_slots=c['marker_slots'],
            recovered_records=c['recovered_records'],unrecovered_sequences=c['marker_slots']-c['recovered_records'],
            model_linked_records=covered,uncovered_recovered_records=c['neither_catalog'],
            model_fraction_of_recovered_records=covered/c['recovered_records'] if c['recovered_records'] else ''))
    write(out/'gap_marker_records.tsv',gaps);write(out/'gap_unique_sequences.tsv',seq_rows)
    write(out/'alternative_cached_models.tsv',candidate_rows,['model_id','version','sequence_sha256','length','provider','tool','path','sha256','uniprot_accessions'])
    write(out/'lineage_marker_availability.tsv',lineage_rows)
    verify()
    result=dict(status='complete_frozen_marker_gap_characterization',inventory_records=records,
        missing_marker_records=len(gaps),missing_unique_sequences=len(sequences),missing_taxa=len({r['taxon_id'] for r in gaps}),
        missing_markers=len({r['marker'] for r in gaps}),link_length_bins=dict(link_counts),unique_length_bins=dict(unique_counts),
        noncanonical_sequences=sum(bool(r['noncanonical_symbols']) for r in seq_rows),
        maximum_missing_sequence_length=max(map(len,sequences.values())),
        alternative_cached_models=len(candidates),missing_sequences_with_alternative_models=len(by_sequence),
        missing_links_with_alternative_models=sum(r['alternative_cached_models']>0 for r in gaps),
        plan_sha256=sha(a.plan),script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in out.iterdir()},
        scope='Full recovered-marker gaps in the specified catalogs, with exact frozen-log checks for alternative cached sources. Length bins are disjoint. Alternative models are not qualified or silently substituted. Not a launched prediction queue, proof of public-database absence, or whole-proteome coverage.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
