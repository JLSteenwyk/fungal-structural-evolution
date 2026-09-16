#!/usr/bin/env python3
"""Read back all merged memberships and crosswalks using the native MCL parser."""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import resource
from orthofinder.tools import mcl


def sha(p):
    with Path(p).open('rb') as h:return hashlib.file_digest(h,'sha256').hexdigest()


def audit(plan_path,output):
    if output.exists():raise FileExistsError(output)
    plan=json.loads(plan_path.read_text())
    for name,h in plan['pins'].items():
        if sha(name)!=h:raise ValueError('Pinned input differs')
    root=Path(plan['output']);receipt=json.loads((root/'receipt.json').read_text())
    discovery_root=Path(plan['discovery_output'])
    if receipt['plan_sha256']!=sha(plan_path) or receipt['discovery_receipt_sha256']!=sha(discovery_root/'receipt.json') or receipt['discovery_readback_sha256']!=sha(discovery_root/'readback.json'):
        raise ValueError('Receipt provenance differs')
    dc=json.loads((discovery_root/'receipt.json').read_text())
    dr=json.loads((discovery_root/'readback.json').read_text())
    if (receipt['status']!='complete_guide_specific_discovery_merge_requires_independent_readback'
        or dc['status']!='complete_all_planned_clade_discovery' or dr['status']!='passed_complete_native_clade_discovery_readback'
        or dc['plan_sha256']!=sha(plan['discovery_plan']) or dr['receipt_sha256']!=sha(discovery_root/'receipt.json')):raise ValueError('Invalid completion chain')
    mapping=json.loads(Path(plan['clade_mapping']).read_text())
    original=mcl.GetPredictedOGs(plan['original_clusters'])
    old={f'OG{i:07d}':g for i,g in enumerate(original)}
    singletons={next(iter(g)):name for name,g in old.items() if len(g)==1}
    retained={name:g for name,g in old.items() if len(g)>1}
    expected=set()
    with Path(plan['sequence_ids']).open() as h:
        for line in h:
            gene=line.split(': ',1)[0]
            if gene in expected:raise ValueError('Duplicate full protein identity')
            expected.add(gene)
    if len(expected)!=plan['expected_proteins']:raise ValueError('Full protein count differs')
    summaries={r['guide']:r for r in receipt['guides']}
    if set(summaries)!=set(mapping['guide_clades']) or len(summaries)!=len(receipt['guides']):raise ValueError('Guide coverage differs')
    reports=[]
    for guide,entries in mapping['guide_clades'].items():
        summary=summaries[guide];folder=root/guide
        for name,h in summary['artifacts'].items():
            if sha(root/name)!=h:raise ValueError('Merged artifact differs')
        source={('retained','',name):g for name,g in retained.items()}
        for entry in entries:
            key=entry['clade_id'];groups=defaultdict(set)
            crpath=discovery_root/key/'receipt.json'
            if sha(crpath)!=dc['clade_receipts'][key]:raise ValueError('Discovery clade receipt changed')
            cr=json.loads(crpath.read_text())
            export=discovery_root/key/'family_membership.tsv'
            if sha(export)!=cr['artifacts']['family_membership.tsv']:raise ValueError('Discovery export changed')
            with export.open() as h:
                for row in csv.DictReader(h,delimiter='\t'):
                    genes=groups[row['local_family']];gene=row['original_native_id']
                    if gene in genes:raise ValueError('Duplicate source membership')
                    genes.add(gene)
            for family,genes in groups.items():
                source[('discovery',key,family)]=genes
        merged=mcl.GetPredictedOGs(str(folder/'clusters_id_pairs.txt'))
        # MCL's parser stores sets, so independently check serialized multiplicity.
        rows=0
        with (folder/'clusters_id_pairs.txt').open() as h:
            if next(h).strip()!='(mclmatrix' or next(h).strip()!='begin':raise ValueError('Unexpected merged header')
            for line in h:
                if line.strip()==')':
                    if any(x.strip() for x in h):raise ValueError('Extra serialized content')
                    break
                fields=line.split()
                if not fields or fields[-1]!='$' or int(fields[0])!=rows or len(fields)-2!=len(merged[rows]):raise ValueError('Serialized indexing or multiplicity differs')
                rows+=1
            else:raise ValueError('Missing closing marker')
        if rows!=len(merged):raise ValueError('Merged row count differs')
        with (folder/'family_sources.tsv').open() as h:crosswalk=list(csv.DictReader(h,delimiter='\t'))
        if len(crosswalk)!=len(merged):raise ValueError('Crosswalk coverage differs')
        seen=set();used_sources=set();expected_relocations={}
        for i,(genes,row) in enumerate(zip(merged,crosswalk)):
            label=f'OG{i:07d}';key=(row['source_type'],row['source_clade'],row['source_family'])
            if row['new_family']!=label or key in used_sources or key not in source or source[key]!=genes:raise ValueError('Source crosswalk membership differs')
            if not genes or seen.intersection(genes):raise ValueError('Repeated merged protein')
            seen.update(genes);used_sources.add(key)
            if len(genes)!=int(row['proteins']):raise ValueError('Family size differs')
            expected_hash=hashlib.sha256(('\n'.join(sorted(genes))+'\n').encode()).hexdigest()
            if row['membership_sha256']!=expected_hash:raise ValueError('Membership digest differs')
            tree_candidate=key[2] if key[0]=='retained' and len(genes)>=3 else ''
            if row['original_tree_candidate']!=tree_candidate:raise ValueError('Original tree candidate differs')
            for gene in genes.intersection(singletons):expected_relocations[gene]=(singletons[gene],label,len(genes))
        if used_sources!=set(source) or seen!=expected:raise ValueError('Source or protein universe incomplete')
        observed_relocations={}
        with (folder/'singleton_relocations.tsv').open() as h:
            for row in csv.DictReader(h,delimiter='\t'):
                gene=row['protein_id']
                if gene in observed_relocations:raise ValueError('Repeated singleton relocation')
                observed_relocations[gene]=(row['original_singleton_family'],row['new_family'],int(row['new_family_size']))
        if observed_relocations!=expected_relocations or set(expected_relocations)!=set(singletons):raise ValueError('Singleton relocation differs')
        if summary['proteins']!=len(seen) or summary['families']!=len(merged) or summary['retained_multi_families']!=len(retained) or summary['discovery_families']!=len(source)-len(retained) or summary['old_singletons_relocated']!=len(singletons):raise ValueError('Summary totals differ')
        reports.append(dict(guide=guide,proteins=len(seen),families=len(merged),retained_multi_families=len(retained),singletons_relocated=len(singletons),all_memberships_and_crosswalks_checked=True))
        del source,merged,crosswalk,seen
    result=dict(status='passed_complete_guide_discovery_merge_readback',guides=reports,plan_sha256=sha(plan_path),merge_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(Path(__file__)),native_reader_sha256=sha(Path(mcl.__file__)),peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,scope='Every serialized merged family checked against original or discovered source memberships, all full-universe protein identities, contiguous indices, source/tree-candidate crosswalks and all singleton relocations. Does not validate biological homology or copy/install gene trees.')
    output.write_text(json.dumps(result,indent=2)+'\n');return result


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    print(json.dumps(audit(a.plan,a.output),indent=2))


if __name__=='__main__':main()
