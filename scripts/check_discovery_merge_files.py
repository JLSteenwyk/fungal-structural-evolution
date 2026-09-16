#!/usr/bin/env python3
"""Run the complete merge command on completed synthetic discovery artifacts."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
from audit_orthology_family_universe import groups


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    fixture=Path('results/orthology/clade-discovery-runner-fixture-v1').resolve()
    readback=Path('metadata/clade_discovery_runner_fixture_readback.json')
    audit=json.loads(readback.read_text())
    if audit['receipt_sha256']!=sha(fixture/'receipt.json'):raise ValueError('Fixture receipt changed')
    # Add the already completed independent readback at the controller's expected path.
    target=fixture/'readback.json'
    if target.exists() and sha(target)!=sha(readback):raise ValueError('Different fixture readback')
    if not target.exists():shutil.copyfile(readback,target)
    expected_core={'500_0','501_0','502_0'}
    expected_discovery={f'{100+s}_{i}' for s in range(2) for i in range(48)}
    ids=out/'SequenceIDs.txt';ids.write_text(''.join(f'{gene}: synthetic_{gene}\n' for gene in sorted(expected_core|expected_discovery)))
    original=out/'original_clusters.txt';original.write_text('(mclmatrix\nbegin\n0 100_0 $\n1 500_0 501_0 502_0 $\n)\n')
    mapping=out/'clades.json';mapping.write_text(json.dumps(dict(clades={'synthetic_1':{},'synthetic_2':{}},guide_clades={g:[dict(clade_id='synthetic_2')] for g in ['profile','mafft']}),indent=2)+'\n')
    dp=Path('metadata/clade_discovery_runner_fixture_plan.json').resolve()
    pins=[Path('scripts/merge_clade_discovery_partitions.py').resolve(),dp,ids,original,mapping,readback.resolve()]
    plan=dict(output=str(out/'merged'),discovery_output=str(fixture),discovery_plan=str(dp),clade_mapping=str(mapping),sequence_ids=str(ids),original_clusters=str(original),expected_proteins=99,expected_original_families=2,expected_original_singletons=1,expected_retained_multi_families=1,pins={str(p):sha(p) for p in pins},resources=dict(minimum_free_disk_gib=1),scope='Synthetic file-workflow test only; both guide labels deliberately share the same synthetic two-taxon discovery output.')
    planpath=out/'plan.json';planpath.write_text(json.dumps(plan,indent=2)+'\n')
    with (out/'stdout.log').open('w') as log:
        subprocess.run([sys.executable,'scripts/merge_clade_discovery_partitions.py','--plan',str(planpath)],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=60)
    expected_groups=[expected_core]+[{f'{100+s}_{2*f+c}' for s in range(2) for c in range(2)} for f in range(24)]
    reports=[]
    for guide in ['profile','mafft']:
        folder=out/'merged'/guide
        actual=list(groups(folder/'clusters_id_pairs.txt'))
        if [i for i,_ in actual]!=list(range(25)) or {frozenset(g) for _,g in actual}!={frozenset(g) for g in expected_groups}:raise ValueError('Serialized merged groups differ')
        with (folder/'family_sources.tsv').open() as h:crosswalk=list(csv.DictReader(h,delimiter='\t'))
        if crosswalk[0]['new_family']!='OG0000000' or crosswalk[0]['source_family']!='OG0000001' or crosswalk[0]['original_tree_candidate']!='OG0000001':raise ValueError('Retained tree source crosswalk differs')
        with (folder/'singleton_relocations.tsv').open() as h:relocations=list(csv.DictReader(h,delimiter='\t'))
        if len(relocations)!=1 or relocations[0]['protein_id']!='100_0' or relocations[0]['original_singleton_family']!='OG0000000' or relocations[0]['new_family_size']!='4':raise ValueError('Singleton relocation differs')
        relocated=int(relocations[0]['new_family'][2:])
        if '100_0' not in actual[relocated][1]:raise ValueError('Singleton crosswalk points to wrong family')
        reports.append(dict(guide=guide,proteins=99,families=25,retained_tree_source_crosswalk_checked=True,singleton_relocation_checked=True))
    result=dict(status='passed_file_based_discovery_merge_fixture',guides=reports,plan_sha256=sha(planpath),merge_receipt_sha256=sha(out/'merged/receipt.json'),script_sha256=sha(Path(__file__)),scope='Complete file-based merge command on independently audited synthetic native discovery artifacts. Both guide labels intentionally share one discovery output. Checks all serialized membership sets, contiguous indices, retained-tree source mapping and replacement of a singleton by a four-member group. Not a full biological merge or homology validation.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
