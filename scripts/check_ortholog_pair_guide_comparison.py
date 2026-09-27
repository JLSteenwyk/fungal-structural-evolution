#!/usr/bin/env python3
"""Independent set fixtures and malformed-stream checks for the native merge."""
import json,random,subprocess,tempfile
from pathlib import Path


def encode(pairs):
    return b''.join(f'{a:06x}{b:06x}0\n{a:06x}{b:06x}1\n'.encode() for a,b in sorted(pairs))


def main():
    rng=random.Random(20260926)
    with tempfile.TemporaryDirectory() as td:
        d=Path(td);exe=d/'merge'
        subprocess.run(['/usr/bin/g++','-O3','-std=c++17','-Wall','-Wextra','-pedantic','scripts/compare_ortholog_pair_streams.cpp','-o',str(exe)],check=True)
        def run(left,right):
            (d/'a').write_bytes(left);(d/'b').write_bytes(right)
            return subprocess.run([str(exe),str(d/'a'),str(d/'b')],capture_output=True,text=True)
        cases=[(set(),set()),({(0,2**24-1)},set()),(set(),{(0,2**24-1)}),({(0,2**24-1)},{(0,2**24-1)})]
        universe=[(a,b) for a in range(20) for b in range(a+1,20)]
        for _ in range(100):cases.append((set(rng.sample(universe,rng.randrange(191))),set(rng.sample(universe,rng.randrange(191)))))
        cases.append(({(0,b) for b in range(1,70001)},{(0,b) for b in range(35000,100001)}))
        for a,b in cases:
            p=run(encode(a),encode(b));assert p.returncode==0,p.stderr;r=json.loads(p.stdout)
            expected=dict(left_pairs=len(a),right_pairs=len(b),shared_pairs=len(a&b),left_only_pairs=len(a-b),right_only_pairs=len(b-a),union_pairs=len(a|b),symmetric_difference_pairs=len(a^b));assert r==expected,(r,expected)
        valid=encode({(0,1)});bad=[valid[:-1],valid+valid,encode({(0,2)})+valid,valid.replace(b'1\n',b'0\n'),b'0000000000000\n0000000000001\n',b'G000000000010\nG000000000011\n',valid[:14]+encode({(0,2)})[14:]]
        for x in bad:assert run(x,b'').returncode!=0
        print(f'Passed {len(cases)} independent-set fixtures including buffer boundary; rejected {len(bad)} malformed streams.')

if __name__=='__main__':main()
