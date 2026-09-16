#!/usr/bin/env python3
"""Exercise singleton replacement and complete-identity gates without biological data."""
import argparse
import hashlib
import json
from pathlib import Path
from merge_clade_discovery_partitions import combine


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    expected={'0_0','1_0','0_1','1_1','2_0','2_1','2_2'}
    original=[('OG0000000',['0_0','1_0']),('OG0000001',['0_1']),('OG0000002',['1_1','2_0'])]
    discovery=[('C_A','OG0000000',['0_1','2_1']),('C_B','OG0000000',['2_2'])]
    merged,singletons=combine(expected,original,discovery)
    reference=[{'0_0','1_0'},{'1_1','2_0'},{'0_1','2_1'},{'2_2'}]
    assert [r[3] for r in merged]==reference and singletons=={'0_1':'OG0000001'}
    assert merged[1][2]=='OG0000002' # Source index differs from new contiguous index.
    alternative=[('C_A','OG0000000',['0_1']),('C_B','OG0000000',['2_1','2_2'])]
    other,_=combine(expected,original,alternative)
    assert [r[3] for r in other]==reference[:2]+[{'0_1'},{'2_1','2_2'}]
    assert set.union(*(r[3] for r in merged))==set.union(*(r[3] for r in other))==expected
    rejected=[]
    cases=[
        ('missing_new_protein',original,discovery[:1]),
        ('retained_overlap',original,discovery+[('C_C','OG0',['0_0'])]),
        ('duplicate_across_clades',original,discovery+[('C_C','OG0',['0_1'])]),
        ('unknown_protein',original,discovery+[('C_C','OG0',['9_9'])]),
        ('duplicate_inside_discovery',original,[('C_A','OG0',['0_1','0_1','2_1']),discovery[1]]),
        ('duplicate_discovery_family_label',original,[discovery[0],('C_A','OG0000000',['2_2'])]),
        ('duplicate_original_gene',original+[('OG3',['0_0'])],discovery),
        ('duplicate_original_family_label',original+[('OG0000000',['2_1'])],discovery),
        ('empty_original_family',original+[('OG3',[])],discovery),
        ('empty_discovery_family',original,discovery+[('C_C','OG0',[])]),
    ]
    for label,old,new in cases:
        try:combine(expected,old,new)
        except ValueError:rejected.append(label)
        else:raise AssertionError('Invalid fixture accepted: '+label)
    result=dict(status='passed_discovery_merge_identity_kernel_fixtures',valid_guides=2,proteins_per_guide=len(expected),retained_multi_families=2,replaced_singleton_slots=1,rejected_cases=rejected,
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                merge_script_sha256=hashlib.sha256(Path('scripts/merge_clade_discovery_partitions.py').read_bytes()).hexdigest(),
                scope='Synthetic identity-set fixtures check retained-family preservation, singleton replacement, guide alternatives, source-index crosswalk and ten invalid partitions. Does not execute a full biological merge, file-based provenance gates or independent merged-file readback.')
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
