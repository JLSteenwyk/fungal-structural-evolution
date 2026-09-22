#!/usr/bin/env python3
"""Exercise full HOG identity/clade audit on hand-constructed source files."""
import json
from pathlib import Path
import tempfile
from audit_hog_output_identities import audit


def main():
    with tempfile.TemporaryDirectory() as directory:
        root=Path(directory);source=root/'source';result=root/'result';source.mkdir()
        tree=result/'Species_Tree';tree.mkdir(parents=True)
        tables=result/'Phylogenetic_Hierarchical_Orthogroups';tables.mkdir()
        (source/'SpeciesIDs.txt').write_text('0: t0.faa\n1: t1.faa\n2: t2.faa\n')
        (source/'SequenceIDs.txt').write_text('0_0: p0\n1_0: p1\n2_0: p2\n')
        (source/'clusters_OrthoFinder.txt_id_pairs.txt').write_text('(mclmatrix\nbegin\n0 0_0 1_0 2_0 $\n)\n')
        (tree/'SpeciesTree_rooted_node_labels.txt').write_text('((t0:1,t1:1)N1:1,t2:1)N0;\n')
        header='HOG\tOG\tGene Tree Parent Clade\tt0\tt1\tt2\n'
        n0=header+'N0.HOG0000000\tOG0000000\tn0\tp0\tp1\tp2\n'
        n1=header+'N1.HOG0000000\tOG0000000\tn1\tp0\tp1\t\n'
        (tables/'N0.tsv').write_text(n0);(tables/'N1.tsv').write_text(n1)
        audit(source,result,root/'good')
        assert json.loads((root/'good/receipt.json').read_text())['gene_assignments_across_levels']==5
        cases=[('wrong_family',n1.replace('OG0000000','OG0000001'),'HOG changes source family'),
               ('outside_clade',n1.rstrip('\n')+'p2\n','Gene outside HOG species clade'),
               ('duplicate_gene',n1.replace('p0\t','p0, p0\t'),'Gene repeated'),
               ('unknown_protein',n1.replace('p1','absent'),'Unknown HOG protein'),
               ('oversized_cell_semantic_check',n1.replace('p0\t',(', '.join(['p0']*40000))+'\t'),'Gene repeated')]
        rejected=[]
        for name,text,error in cases:
            (tables/'N1.tsv').write_text(text)
            try:audit(source,result,root/name)
            except ValueError as exc:
                assert error in str(exc),str(exc);rejected.append(name)
            else:raise AssertionError(name+' accepted')
        (tables/'N1.tsv').unlink()
        try:audit(source,result,root/'missing')
        except ValueError as exc:
            assert 'Missing or extra' in str(exc);rejected.append('missing_node_table')
        else:raise AssertionError('Missing table accepted')
        print(json.dumps(dict(status='passed',rejected_corruptions=rejected,scope='Hand-constructed identities, family, clade and repeated-assignment fixture; not biological HOG validation.')))


if __name__=='__main__':main()
