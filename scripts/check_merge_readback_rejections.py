#!/usr/bin/env python3
"""Ensure merge readback rejects semantic corruption after checksums are updated."""
import argparse
import csv
import json
from pathlib import Path
import shutil
from readback_guide_discovery_merge import audit,sha


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    base=Path('results/orthology/discovery-merge-fixture-v1');original=json.loads((base/'plan.json').read_text());results=[]
    for case in ['wrong_source_family','wrong_tree_candidate','wrong_singleton_target','duplicate_serialized_gene','missing_serialized_gene']:
        folder=out/case;folder.mkdir();merged=folder/'merged';shutil.copytree(base/'merged',merged)
        if case in ['wrong_source_family','wrong_tree_candidate']:
            path=merged/'profile/family_sources.tsv'
            with path.open() as h:rows=list(csv.DictReader(h,delimiter='\t'))
            rows[0]['source_family' if case=='wrong_source_family' else 'original_tree_candidate']='OG9999999'
            with path.open('w') as h:
                w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
        elif case=='wrong_singleton_target':
            path=merged/'profile/singleton_relocations.tsv'
            with path.open() as h:rows=list(csv.DictReader(h,delimiter='\t'))
            rows[0]['new_family']='OG0000000'
            with path.open('w') as h:
                w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
        else:
            path=merged/'profile/clusters_id_pairs.txt';lines=path.read_text().splitlines();parts=lines[2].split()
            if case=='duplicate_serialized_gene':parts.insert(1,parts[1])
            else:parts.pop(1)
            lines[2]=' '.join(parts);path.write_text('\n'.join(lines)+'\n')
        plan=dict(original,output=str(merged));planpath=folder/'plan.json';planpath.write_text(json.dumps(plan,indent=2)+'\n')
        receipt=json.loads((merged/'receipt.json').read_text());receipt['plan_sha256']=sha(planpath)
        for report in receipt['guides']:
            report['artifacts']={name:sha(merged/name) for name in report['artifacts']}
        (merged/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
        try:audit(planpath,folder/'readback.json')
        except ValueError as exc:
            if (folder/'readback.json').exists():raise AssertionError('Failed case wrote success output')
            results.append(dict(case=case,rejected=True,reason=str(exc)))
        else:raise AssertionError('Corrupt fixture accepted: '+case)
    result=dict(status='passed_rechecksummed_merge_corruption_fixtures',cases=results,script_sha256=sha(Path(__file__)),auditor_sha256=sha(Path('scripts/readback_guide_discovery_merge.py')),scope='Five synthetic output corruptions with updated artifact/plan digests test semantic validation rather than checksum rejection. Original biological and synthetic source artifacts are untouched.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
