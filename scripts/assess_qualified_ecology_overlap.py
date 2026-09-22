#!/usr/bin/env python3
"""Measure common qualified alignment columns among curated ecological taxa."""
import argparse
import csv
import itertools
import json
from pathlib import Path
from Bio import SeqIO
from readback_whole_proteome_family_coverage import sha


def load(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def common_positions(position_sets):
    return set.intersection(*position_sets) if position_sets else set()


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text())
    for path,digest in plan['pins'].items():
        if sha(path)!=digest:raise ValueError('Changed dependency: '+path)
    out=Path(plan['output'])
    if out.exists():raise FileExistsError(out)
    evidence=load(plan['evidence']);taxa=[r['taxon_id'] for r in evidence]
    assert len(taxa)==len(set(taxa))
    evidence_receipt=json.loads(Path(plan['evidence_receipt']).read_text())
    assert evidence_receipt['table_sha256']==sha(plan['evidence'])
    base=Path(plan['inputs']);receipt=json.loads((base/'receipt.json').read_text())
    readback=json.loads(Path(plan['input_readback']).read_text())
    assert readback['status']=='passed_complete_paired_inputs_from_qualified_arrays_readback'
    assert readback['source_receipt_sha256']==sha(base/'receipt.json')
    for name,digest in receipt['artifacts'].items():assert sha(base/name)==digest
    masks={};species_counts={t:dict(eligible_markers=0,qualified_observations=0) for t in taxa}
    for name in sorted(receipt['artifacts']):
        if not name.endswith('/aa.faa'):continue
        marker=name.split('/')[0]
        aa=SeqIO.to_dict(SeqIO.parse(base/name,'fasta'))
        di=SeqIO.to_dict(SeqIO.parse(base/marker/'3di.faa','fasta'))
        assert aa.keys()==di.keys()
        lengths={len(r) for r in aa.values()};assert len(lengths)==1
        marker_masks={}
        for taxon in aa:
            a=str(aa[taxon].seq);d=str(di[taxon].seq)
            assert len(a)==len(d) and set(a)<=set('ACDEFGHIKLMNPQRSTVWY?') and set(d)<=set('ACDEFGHIKLMNPQRSTVWY?')
            apos={i for i,c in enumerate(a) if c!='?'};dpos={i for i,c in enumerate(d) if c!='?'}
            assert apos==dpos
            if taxon in species_counts:
                marker_masks[taxon]=apos
                species_counts[taxon]['eligible_markers']+=1
                species_counts[taxon]['qualified_observations']+=len(apos)
        masks[marker]=marker_masks
    assert len(masks)==readback['markers']
    pairs=[];details=[]
    for a,b in itertools.combinations(evidence,2):
        x,y=a['taxon_id'],b['taxon_id'];shared=[];observations=0;both_eligible=0;nonzero=0
        for marker,positions in masks.items():
            if x not in positions or y not in positions:continue
            both_eligible+=1;n=len(positions[x]&positions[y]);observations+=n;nonzero+=n>0
            if n>=plan['minimum_shared_columns']:shared.append(marker)
            details.append(dict(taxon_a=x,taxon_b=y,marker=marker,shared_qualified_columns=n))
        pairs.append(dict(taxon_a=x,species_a=a['species_name'],state_a=a['state'],taxon_b=y,
                          species_b=b['species_name'],state_b=b['state'],
                          different_published_states=a['state']!=b['state'],
                          both_eligible_markers=both_eligible,markers_with_shared_columns=nonzero,
                          markers_with_at_least_50_shared_columns=len(shared),
                          shared_qualified_columns=observations,marker_ids_at_least_50=';'.join(shared)))
    groups=[]
    for group in sorted({r['provisional_transition_group'] for r in evidence}):
        members=[r['taxon_id'] for r in evidence if r['provisional_transition_group']==group]
        eligible=nonzero=qualified=observations=0
        for positions in masks.values():
            if not all(t in positions for t in members):continue
            eligible+=1;n=len(common_positions([positions[t] for t in members]));observations+=n
            nonzero+=n>0;qualified+=n>=plan['minimum_shared_columns']
        groups.append(dict(provisional_group=group,taxa=len(members),taxon_ids=';'.join(members),
                           all_taxa_eligible_markers=eligible,markers_with_common_columns=nonzero,
                           markers_with_at_least_50_common_columns=qualified,common_qualified_columns=observations))
    per_taxon=[dict(taxon_id=r['taxon_id'],species_name=r['species_name'],state=r['state'],**species_counts[r['taxon_id']]) for r in evidence]
    out.mkdir(parents=True)
    artifacts={}
    for name,rows in [('taxa.tsv',per_taxon),('pairs.tsv',pairs),('pair_marker_columns.tsv',details),('groups.tsv',groups)]:
        with (out/name).open('w') as handle:
            w=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
        artifacts[name]=sha(out/name)
    result=dict(status='complete_qualified_esmfold_ecology_overlap',taxa=len(taxa),
                taxa_with_eligible_markers=sum(v['eligible_markers']>0 for v in species_counts.values()),
                markers=len(masks),pairs=len(pairs),provisional_groups=len(groups),
                pairs_with_at_least_one_50_column_marker=sum(r['markers_with_at_least_50_shared_columns']>0 for r in pairs),
                plan_sha256=sha(args.plan),artifacts=artifacts,
                scope='Shared identical alignment columns under existing AA/3Di confidence and eligibility masks. ESMFold only, minimum 50 common columns is descriptive, not a power criterion. Taxon pairs and provisional groups are not independent transitions; no ecological effect or orthology claim.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
