#!/usr/bin/env python3
"""Merge complete guide-specific discovery with retained multi-protein families."""
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import shutil
from audit_orthology_family_universe import groups


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def combine(expected,original,discovery):
    """Validate first, then return contiguous families with source provenance."""
    original_seen=set();source_labels=set();retained=[];old_singletons={}
    for family,genes in original:
        gene_set=set(genes)
        if family in source_labels or not gene_set or len(gene_set)!=len(genes) or original_seen.intersection(gene_set) or not gene_set<=expected:
            raise ValueError('Invalid original family partition')
        source_labels.add(family);original_seen.update(gene_set)
        if len(gene_set)==1:old_singletons[next(iter(gene_set))]=family
        else:retained.append(('retained', '',family,gene_set))
    retained_genes=set().union(*(r[3] for r in retained)) if retained else set()
    seen=set(retained_genes);new=[];labels=set()
    for clade,family,genes in discovery:
        gene_set=set(genes)
        if (clade,family) in labels or not gene_set or len(gene_set)!=len(genes) or seen.intersection(gene_set) or not gene_set<=expected:
            raise ValueError('Invalid discovery partition or overlap with retained families')
        labels.add((clade,family));seen.update(gene_set);new.append(('discovery',clade,family,gene_set))
    if seen!=expected:raise ValueError('Merged partition omits input proteins')
    return retained+new,old_singletons


def exported_groups(path):
    with path.open() as handle:
        reader=csv.DictReader(handle,delimiter='\t')
        for family,rows in itertools.groupby(reader,key=lambda r:r['local_family']):
            yield family,[r['original_native_id'] for r in rows]


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args()
    plan=json.loads(a.plan.read_text())
    for name,h in plan['pins'].items():
        if sha(name)!=h:raise ValueError('Changed pinned input '+name)
    output=Path(plan['output'])
    if output.exists():raise FileExistsError(output)
    discovery_root=Path(plan['discovery_output'])
    receipt=json.loads((discovery_root/'receipt.json').read_text())
    readback=json.loads((discovery_root/'readback.json').read_text())
    dp=json.loads(Path(plan['discovery_plan']).read_text());dp_hash=sha(plan['discovery_plan'])
    keys={c['clade_id'] for c in dp['cases']}
    if (receipt['status']!='complete_all_planned_clade_discovery' or set(receipt['clade_receipts'])!=keys
        or receipt['plan_sha256']!=dp_hash or readback['status']!='passed_complete_native_clade_discovery_readback'
        or readback['plan_sha256']!=dp_hash or readback['receipt_sha256']!=sha(discovery_root/'receipt.json')
        or {c['clade_id'] for c in readback['clades']}!=keys):raise ValueError('Complete independent discovery readback required')
    mapping=json.loads(Path(plan['clade_mapping']).read_text())
    if set(mapping['clades'])!=keys:raise ValueError('Guide mapping differs from completed cases')
    if shutil.disk_usage(output.parent).free<plan['resources']['minimum_free_disk_gib']*2**30:raise RuntimeError('Disk gate')
    expected=set()
    with Path(plan['sequence_ids']).open() as handle:
        for line in handle:
            gene,sep,label=line.rstrip('\n').partition(': ')
            if not sep or not label or gene in expected:raise ValueError('Invalid full protein mapping')
            expected.add(gene)
    if len(expected)!=plan['expected_proteins']:raise ValueError('Protein universe size differs')
    original=[(f'OG{number:07d}',genes) for number,genes in groups(Path(plan['original_clusters']))]
    if len(original)!=plan['expected_original_families']:raise ValueError('Original family count differs')
    for key,h in receipt['clade_receipts'].items():
        path=discovery_root/key/'receipt.json'
        if sha(path)!=h:raise ValueError('Clade receipt changed')
        cr=json.loads(path.read_text())
        if cr['plan_sha256']!=dp_hash:raise ValueError('Clade plan differs')
        for name,digest in cr['artifacts'].items():
            if sha(path.parent/name)!=digest:raise ValueError('Clade artifact changed')
    output.mkdir(parents=True);summaries=[]
    for guide,entries in mapping['guide_clades'].items():
        selected=[entry['clade_id'] for entry in entries]
        if len(selected)!=len(set(selected)):raise ValueError('Repeated guide clade')
        def discovered():
            for key in sorted(selected):
                for family,genes in exported_groups(discovery_root/key/'family_membership.tsv'):
                    yield key,family,genes
        merged,singletons=combine(expected,original,discovered())
        if len(singletons)!=plan['expected_original_singletons']:raise ValueError('Singleton count differs')
        folder=output/guide;folder.mkdir()
        paths=[folder/'clusters_id_pairs.txt',folder/'family_sources.tsv',folder/'singleton_relocations.tsv']
        reassigned=0;retained_count=0
        with paths[0].open('w') as cluster,paths[1].open('w') as sources,paths[2].open('w') as relocations:
            source_writer=csv.writer(sources,delimiter='\t');source_writer.writerow(['new_family','source_type','source_clade','source_family','proteins','membership_sha256','original_tree_candidate'])
            relocation_writer=csv.writer(relocations,delimiter='\t');relocation_writer.writerow(['original_singleton_family','protein_id','new_family','new_family_size'])
            cluster.write('(mclmatrix\nbegin\n')
            for number,(kind,clade,family,genes) in enumerate(merged):
                label=f'OG{number:07d}';ordered=sorted(genes)
                cluster.write(str(number)+' '+' '.join(ordered)+' $\n')
                membership_hash=hashlib.sha256(('\n'.join(ordered)+'\n').encode()).hexdigest()
                source_writer.writerow([label,kind,clade,family,len(genes),membership_hash,family if kind=='retained' and len(genes)>=3 else ''])
                retained_count+=kind=='retained'
                for gene in ordered:
                    if gene in singletons:relocation_writer.writerow([singletons[gene],gene,label,len(genes)]);reassigned+=1
            cluster.write(')\n')
        if retained_count!=plan['expected_retained_multi_families'] or reassigned!=len(singletons):raise ValueError('Merge accounting differs')
        summaries.append(dict(guide=guide,proteins=len(expected),families=len(merged),retained_multi_families=retained_count,discovery_families=len(merged)-retained_count,old_singletons_relocated=reassigned,artifacts={str(p.relative_to(output)):sha(p) for p in paths}))
        del merged
    result=dict(status='complete_guide_specific_discovery_merge_requires_independent_readback',guides=summaries,plan_sha256=sha(a.plan),script_sha256=sha(Path(__file__)),discovery_receipt_sha256=sha(discovery_root/'receipt.json'),discovery_readback_sha256=sha(discovery_root/'readback.json'),interpretation='Full protein partitions with unique membership, retaining every original multi-protein family and replacing all original singleton slots through complete guide-specific discovery. Family indices are renumbered contiguously; source and singleton crosswalks are mandatory. Original trees are candidates only and are not copied or assumed valid. Independent merged-partition readback, added trees and reconciliation remain required.')
    (output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
