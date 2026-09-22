#!/usr/bin/env python3
"""Check native resolved-tree coverage and reject semantic output corruption."""
import json
from pathlib import Path
import tempfile
from readback_resolved_tree_memberships import audit, tree_membership

ROOT=Path(__file__).resolve().parents[1]


def main():
    fixture=ROOT/'results/orthology/native-restart-traced-six-hour-fixture-v1/prepared_ids_tree'
    source=fixture/'Source/WorkingDirectory';result=fixture/'Results_native_fixture'
    text=(result/'Resolved_Gene_Trees/Resolved_Gene_Trees.txt').read_text()
    cases=dict(valid=text,missing='\n'.join(text.splitlines()[:-1])+'\n',duplicate=text+text.splitlines()[0]+'\n',
        wrong_family_member=text.replace('Taxon0_protein_0_0','Taxon0_protein_0_1',1),
        duplicate_leaf=text.replace('Taxon1_protein_1_0','Taxon0_protein_0_0',1),
        negative_branch=text.replace(':0.1',':-0.1',1),missing_branch=text.replace(':0.1','',1))
    checks=[]
    with tempfile.TemporaryDirectory(prefix='resolved-tree-check-') as tmp:
        root=Path(tmp)
        for name,content in cases.items():
            output=root/name/'Resolved_Gene_Trees';output.mkdir(parents=True)
            (output/'Resolved_Gene_Trees.txt').write_text(content)
            receipt=root/(name+'.json')
            try:r=audit(source,root/name,receipt)
            except ValueError:
                assert name!='valid' and not receipt.exists();checks.append(dict(case=name,passed=True));continue
            assert name=='valid' and r['resolved_trees']==3 and r['resolved_tips']==13,name
            checks.append(dict(case=name,passed=True))
        # An explicitly deep ladder exercises iterative traversal beyond the
        # default Python recursion limit without changing the process limit.
        newick='g0:1'
        for i in range(1,1500):newick='('+newick+',g'+str(i)+':1):1'
        tips,_,_=tree_membership(newick+';',{f'g{i}':f'0_{i}' for i in range(1500)})
        assert tips==1500
    print(json.dumps(dict(status='passed_resolved_tree_readback_fixture_checks',cases=checks,deep_tree_leaves=1500),indent=2))


if __name__=='__main__':main()
