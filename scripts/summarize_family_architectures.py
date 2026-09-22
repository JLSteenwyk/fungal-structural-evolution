#!/usr/bin/env python3
"""Inventory within-family annotation variation across all competition policies."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import sqlite3
import time

POLICIES = ('alignment_evalue','alignment_bitscore','envelope_evalue','envelope_bitscore')


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8388608),b''):
            h.update(block)
    return h.hexdigest()


class Variation:
    def __init__(self):
        self.counts=Counter()
        self.taxa=set()
        self.signatures=set()
        self.bags=defaultdict(set)
        self.sets=defaultdict(set)
        self.by_taxon=defaultdict(set)

    def add(self,taxon,tokens):
        ordered=tuple(tuple(t) for t in tokens)
        if not ordered:
            raise ValueError('No-hit must not be treated as an architecture')
        bag=tuple(sorted(Counter(ordered).items()))
        kinds=tuple(sorted(set(ordered)))
        self.counts['proteins']+=1
        self.taxa.add(taxon)
        self.signatures.add(ordered)
        self.bags[bag].add(ordered)
        self.sets[kinds].add(bag)
        self.by_taxon[taxon].add(ordered)

    def result(self,prefix):
        return {prefix+k:v for k,v in dict(proteins=self.counts['proteins'],taxa=len(self.taxa),
                ordered_signatures=len(self.signatures),multisets=len(self.bags),model_type_sets=len(self.sets),
                multisets_with_order_variation=sum(len(v)>1 for v in self.bags.values()),
                model_type_sets_with_multiplicity_variation=sum(len(v)>1 for v in self.sets.values()),
                taxa_with_multiple_signatures=sum(len(v)>1 for v in self.by_taxon.values())).items()}


def summarize(rows):
    counts={p:Counter() for p in POLICIES}
    observed={p:Variation() for p in POLICIES}
    conservative={p:Variation() for p in POLICIES}
    taxa=Counter()
    sequences=set()
    for taxon,sequence,raw_hits,agree,text in rows:
        taxa[taxon]+=1
        sequences.add(sequence)
        payload=json.loads(text)
        for policy in POLICIES:
            c=counts[policy]
            c['proteins']+=1
            c['policy_disagreement_proteins']+=not bool(agree)
            if not raw_hits:
                c['no_ga_hit_proteins']+=1
                continue
            entry=payload['policies'][policy]
            alt=payload['alternatives'][entry['alternative_index']]
            tokens=alt['ordered_model_tokens']
            if not tokens:
                raise ValueError('Raw hits but empty candidate architecture')
            observed[policy].add(taxon,tokens)
            flags=dict(rank_tie_proteins=entry['primary_rank_tie_pairs']>0,
                       unresolved_overlap_proteins=entry['unresolved_overlap_pairs']>0,
                       candidate_nested_proteins=entry['candidate_nested_pairs']>0,
                       alignment_overlap_proteins=alt['alignment_overlap_pairs']>0,
                       partial_hmm_proteins=alt['partial_hmm_hits_below_070']>0)
            c.update(flags)
            if agree and not any(flags.values()):
                conservative[policy].add(taxon,tokens)
    result=[]
    for policy in POLICIES:
        c=counts[policy]
        row=dict(policy=policy,proteins=c['proteins'],taxa=len(taxa),unique_sequences=len(sequences),
                 taxa_with_multiple_family_members=sum(n>1 for n in taxa.values()),
                 maximum_members_per_taxon=max(taxa.values(),default=0))
        for key in ['no_ga_hit_proteins','policy_disagreement_proteins','rank_tie_proteins',
                    'unresolved_overlap_proteins','candidate_nested_proteins','alignment_overlap_proteins','partial_hmm_proteins']:
            row[key]=c[key]
        row.update(observed[policy].result('observed_'))
        row.update(conservative[policy].result('conservative_'))
        assert row['observed_proteins']+row['no_ga_hit_proteins']==row['proteins']
        result.append(row)
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan',type=Path,required=True)
    a=ap.parse_args()
    plan=json.loads(a.plan.read_text())
    pins={str(a.plan):sha(a.plan),**plan['pins']}
    def verify():
        for path,digest in pins.items():
            if sha(path)!=digest:
                raise ValueError('Changed source: '+path)
    verify()
    out=Path(plan['output'])
    if out.exists():
        raise FileExistsError(out)
    if shutil.disk_usage(out.parent).free<plan['resources']['minimum_free_disk_gib']*2**30:
        raise ValueError('Insufficient disk')
    out.mkdir()
    started=time.time()
    def state(stage,**counts):
        p=out/'state.tmp'
        p.write_text(json.dumps(dict(stage=stage,elapsed_seconds=time.time()-started,**counts))+'\n')
        p.replace(out/'state.json')
    state('waiting_for_exact_bridge_auditor')
    proc=Path('/proc')/str(plan['auditor_pid'])/'stat'
    while proc.exists():
        try:
            fields=proc.read_text().rsplit(')',1)[1].split()
        except FileNotFoundError:
            break
        if fields[19]!=str(plan['auditor_start_ticks']) or fields[0]=='Z':
            break
        time.sleep(20)
    audited=json.loads(Path(plan['bridge_readback']).read_text())
    source=json.loads(Path(plan['bridge_receipt']).read_text())
    if (audited['status']!='passed_complete_independent_family_domain_bridge_readback'
            or audited['producer_receipt_sha256']!=sha(plan['bridge_receipt'])
            or source['bridge_sha256']!=sha(plan['bridge_database'])
            or source['architecture_database_sha256']!=sha(plan['architecture_database'])):
        raise ValueError('Unverified bridge or architecture database')
    for key in ['bridge_readback','bridge_receipt','bridge_database','architecture_database']:
        pins[plan[key]]=sha(plan[key])
    verify()
    db=sqlite3.connect('file:'+str(Path(plan['bridge_database']).resolve())+'?mode=ro',uri=True)
    db.execute('ATTACH DATABASE ? AS architecture',('file:'+str(Path(plan['architecture_database']).resolve())+'?mode=ro',))
    summaries=[]
    artifacts={}
    for guide in source['guides']:
        name=guide['guide']
        state('family_variation',guide=name)
        sql='''SELECT a.family,p.taxon_id,p.sequence_id,q.raw_hits,q.retained_sets_agree,q.candidate_architectures_json
               FROM assignments a INDEXED BY assignments_family
               JOIN proteins p ON p.native_gene_id=a.native_gene_id
               JOIN architecture.queries q ON q.sequence_id=p.sequence_id
               WHERE a.guide=? ORDER BY a.family'''
        counts=Counter()
        path=out/(name+'_family_architectures.tsv.gz')
        with gzip.open(path,'wt') as f:
            writer=None
            for family,rows in itertools.groupby(db.execute(sql,(name,)),key=lambda r:r[0]):
                metrics=summarize((r[1:] for r in rows))
                for row in metrics:
                    row={'guide':name,'family':family,**row}
                    if writer is None:
                        writer=csv.DictWriter(f,fieldnames=list(row),delimiter='\t')
                        writer.writeheader()
                    writer.writerow(row)
                    counts['rows']+=1
                counts['families']+=1
                counts['proteins']+=metrics[0]['proteins']
                if counts['families']%1000==0:
                    state('family_variation',guide=name,**counts)
            if counts['families']!=guide['families'] or counts['proteins']!=source['proteins']:
                raise ValueError('Family inventory scope differs')
        summaries.append(dict(guide=name,**counts))
        artifacts[path.name]=sha(path)
    db.close()
    verify()
    result=dict(status='complete_family_architecture_variation_inventory',guides=summaries,artifacts=artifacts,
                input_hashes=pins,elapsed_seconds=time.time()-started,
                scope='Descriptive within-family annotation variation across four policies, preserving no-hit and QC flags. Conservative excludes policy disagreement, partial matches, ties, overlap and candidate nesting. No event directions, branch rates, homology validation or statistical significance inferred.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    state(result['status'])


if __name__=='__main__':
    main()
