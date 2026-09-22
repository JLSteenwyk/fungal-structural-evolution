#!/usr/bin/env python3
"""Check exact membership controls despite guide-specific family identifiers."""
import csv
import gzip
from pathlib import Path
import sqlite3
import tempfile
from check_family_architecture_readback import fixture
from summarize_family_architectures import summarize
from compare_family_architecture_guides import compare


def main():
    data = fixture()
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp); db = sqlite3.connect(':memory:')
        db.execute('CREATE TABLE families(guide TEXT,family TEXT,membership_sha256 TEXT)')
        names = {'profile':['OG0','OG1'], 'mafft':['OG9','OG8']}
        for guide, ids in names.items():
            for i, name in enumerate(ids):
                db.execute('INSERT INTO families VALUES (?,?,?)',(guide,name,'shared' if i==0 else guide))
        summaries = [dict(guide=g,families=2,proteins=6,rows=8) for g in names]
        original = {}
        for guide, ids in names.items():
            original[guide] = [dict(guide=guide,family=ids[i],**r)
                               for i in range(2) for r in summarize(data[i*3:i*3+3])]
        def write(variant):
            for guide, rr in original.items():
                rows = [dict(r) for r in rr]
                if guide=='mafft':
                    if variant=='changed': rows[0]['observed_multisets'] += 1
                    if variant=='missing': rows.pop()
                    if variant=='duplicate': rows[-1] = rows[0]
                with gzip.open(root/(guide+'_family_architectures.tsv.gz'),'wt') as f:
                    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
        for variant in ['valid','changed','missing','duplicate']:
            write(variant)
            try: rows,shared,checked=compare(root,db,summaries)
            except ValueError: assert variant!='valid'
            else:
                assert variant=='valid'
                assert (len(rows),shared,checked)==(24,1,4)
                for r in rows:
                    assert r['families']==(2 if r['membership_scope']=='all_families' else 1)
                    assert r['proteins']==(6 if r['membership_scope']=='all_families' else 3)
                    if r['membership_scope']=='exact_shared_membership':
                        assert r['observed_families_with_multisets_with_order_variation']==1
                        assert r['observed_families_with_model_type_sets_with_multiplicity_variation']==1
    print('PASS: renumbered exact membership, partition totals, order/multiplicity summaries, changed/missing/duplicate rejection')


if __name__=='__main__':
    main()
