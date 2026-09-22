#!/usr/bin/env python3
"""Inventory small families exposed to two observed native output omissions."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8388608),b''):h.update(block)
    return h.hexdigest()


def groups(path):
    started=False
    pending=[]
    index=0
    with path.open() as f:
        for line in f:
            if not started:
                if 'begin' in line:started=True
                continue
            if line.strip()==')':
                if pending:raise ValueError('Unterminated group')
                return
            tokens=line.split()
            if not pending:
                if int(tokens.pop(0))!=index:raise ValueError('Nonconsecutive group index')
            pending.extend(tokens)
            if pending[-1]=='$':
                genes=pending[:-1]
                if not genes or len(set(genes))!=len(genes):raise ValueError('Invalid membership')
                yield f'OG{index:07d}',genes
                pending=[];index+=1
    raise ValueError('Missing cluster terminator')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    receipt=json.loads((a.inputs/'receipt.json').read_text())
    summaries=[];pins={str(a.inputs/'receipt.json'):sha(a.inputs/'receipt.json')};artifacts={}
    for summary in receipt['guides']:
        guide=summary['guide'];folder=a.inputs/guide
        manifest=folder/'copied_files.tsv'
        if sha(manifest)!=summary['manifest_sha256']:raise ValueError('Changed input manifest')
        records={row['relative_path']:row for row in csv.DictReader(manifest.open(),delimiter='\t')}
        cluster=folder/'Source/WorkingDirectory/clusters_OrthoFinder.txt_id_pairs.txt'
        if sha(cluster)!=records[str(cluster.relative_to(folder))]['sha256']:raise ValueError('Changed clusters')
        crosswalk=folder/'family_sources.tsv'
        pins.update({str(q):sha(q) for q in [manifest,cluster,crosswalk]})
        counts=Counter();first_singleton=None;path=a.output/(guide+'_exposed_small_families.tsv')
        with path.open('w') as out,crosswalk.open() as cw:
            table=iter(csv.DictReader(cw,delimiter='\t'))
            w=csv.writer(out,delimiter='\t',lineterminator='\n')
            w.writerow(['family','genes','taxa','reason','expected_directed_pairs','membership_sha256'])
            for family,genes in groups(cluster):
                row=next(table)
                digest=hashlib.sha256(('\n'.join(sorted(genes))+'\n').encode()).hexdigest()
                if row['new_family']!=family or int(row['proteins'])!=len(genes) or row['membership_sha256']!=digest:
                    raise ValueError('Crosswalk mismatch')
                counts['families']+=1;counts['proteins']+=len(genes)
                if len(genes)==1 and first_singleton is None:first_singleton=family
                if len(genes) not in [2,3]:continue
                species=[gene.split('_')[0] for gene in genes]
                taxa=len(set(species))
                counts[f'genes{len(genes)}_taxa{taxa}']+=1
                if taxa<2:continue
                reason=('two_gene_pair_not_appended' if len(genes)==2 else
                        'three_gene_family_after_first_singleton' if first_singleton else None)
                if reason:
                    pairs=sum(x!=y for x in species for y in species)
                    counts[reason]+=1;counts['exposed_families']+=1;counts['expected_missing_directed_pairs']+=pairs
                    w.writerow([family,len(genes),taxa,reason,pairs,digest])
            if next(table,None) is not None:raise ValueError('Extra crosswalk row')
        if counts['families']!=summary['families'] or counts['proteins']!=summary['proteins']:
            raise ValueError('Incomplete partition')
        artifacts[path.name]=sha(path)
        summaries.append(dict(guide=guide,first_singleton=first_singleton,**counts))
    for path,digest in pins.items():
        if sha(path)!=digest:raise ValueError('Input changed during assessment')
    result=dict(status='complete_input_exposure_assessment_pending_native_output_readback',guides=summaries,
                input_hashes=pins,artifacts=artifacts,script_sha256=sha(__file__),
                scope='Expected output omissions inferred from installed-code behavior reproduced on synthetic fixtures. Does not assert observed omission counts in still-running production results, repair native outputs, or establish biological orthology.')
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(summaries,indent=2))


if __name__=='__main__':main()
