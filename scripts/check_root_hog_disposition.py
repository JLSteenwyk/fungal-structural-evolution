#!/usr/bin/env python3
"""Full CLI checks for explained and unexplained root HOG omissions."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    with tempfile.TemporaryDirectory() as folder:
        root=Path(folder);source=root/'source';source.mkdir();result=root/'result';hogs=result/'Phylogenetic_Hierarchical_Orthogroups';flags=result/'Phylogenetically_Misplaced_Genes';hogs.mkdir(parents=True);flags.mkdir()
        (source/'SpeciesIDs.txt').write_text('0: t0.faa\n1: t1.faa\n')
        (source/'SequenceIDs.txt').write_text('0_0: p0\n1_0: p1\n1_1: p2\n0_1: p3\n')
        (source/'clusters_OrthoFinder.txt_id_pairs.txt').write_text('(mclmatrix\nbegin\n0 0_0 1_0 1_1 $\n1 0_1 $\n)\n')
        original='HOG\tOG\tGene Tree Parent Clade\tt0\tt1\nN0.HOG0000000\tOG0000000\tn0\tp0\tp1\n'
        (hogs/'N0.tsv').write_text(original);(flags/'t1.txt').write_text('p2\n')
        def run(name):
            plan={'source':str(source),'result':str(result),'output':str(root/name),'pins':{}}
            path=root/(name+'.json');path.write_text(json.dumps(plan))
            return subprocess.run([sys.executable,'scripts/resolve_root_hog_gene_disposition.py','--plan',str(path)],capture_output=True,text=True)
        r=run('explained');assert r.returncode==0,r.stderr
        d=json.loads((root/'explained/receipt.json').read_text());assert d['missing_exactly_matches_native_flagged'] and d['singleton_genes']==1 and d['root_assigned_genes']==2
        (flags/'t1.txt').write_text('')
        r=run('unexplained');assert r.returncode==0,r.stderr
        d=json.loads((root/'unexplained/receipt.json').read_text());assert not d['missing_exactly_matches_native_flagged'] and d['unexplained_missing']==1
        (flags/'t1.txt').write_text('unknown\n');r=run('badflag');assert r.returncode!=0 and 'unknown flagged genes' in r.stderr
        (flags/'t1.txt').write_text('p2\n');(hogs/'N0.tsv').write_text(original.replace('OG0000000','OG0000001'))
        r=run('badfamily');assert r.returncode!=0 and 'wrong-family' in r.stderr
        print(json.dumps(dict(status='passed_root_hog_disposition_fixtures',checks=['exact_flagged_omission','unexplained_omission_retained','unknown_flag_rejected','wrong_family_rejected'])))


if __name__=='__main__':main()
