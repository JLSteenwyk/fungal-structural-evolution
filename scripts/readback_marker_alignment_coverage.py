#!/usr/bin/env python3
"""Check every coverage cell using independent FASTA parsing and NumPy counts."""
import csv
import hashlib
import json
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    root=Path('results/phylogeny/marker-alignment-coverage-20260927-v1');receipt=json.loads((root/'receipt.json').read_text())
    for p,d in receipt['sources'].items():assert sha(p)==d
    for p,d in receipt['artifacts'].items():assert sha(root/p)==d
    computed={};taxa=None
    for label in ['profile','mafft']:
        folder=Path('results/phylogeny')/(label+'-matrix-50-v1')
        source=json.loads((folder/'receipt.json').read_text())
        for name,d in source['artifacts'].items():assert sha(folder/name)==d
        seqs={};key=None
        for line in (folder/'matrix.faa').read_text().splitlines():
            if line.startswith('>'):key=line[1:].split()[0];assert key not in seqs;seqs[key]=''
            else:seqs[key]+=line.strip()
        names=sorted(seqs)
        if taxa is None:taxa=names
        assert taxa==names
        matrix=np.array([np.frombuffer(seqs[t].encode(),dtype=np.uint8) for t in names])
        observed=np.isin(matrix,np.frombuffer(b'ACDEFGHIKLMNPQRSTVWY',dtype=np.uint8))
        columns=defaultdict(list)
        for r in read(folder/'site_mapping.tsv'):columns[r['marker']].append(int(r['matrix_column_1based'])-1)
        for marker,cols in columns.items():
            counts=observed[:,cols].sum(axis=1);minimum=max(50,(len(cols)*3+9)//10)
            computed[label,marker]={t:(len(cols),minimum,int(n)) for t,n in zip(names,counts)}
    rows=read(root/'marker_taxa.tsv');seen=set();groups=defaultdict(dict);counts=Counter()
    for r in rows:
        marker,taxon=r['marker'],r['taxon_id'];assert (marker,taxon) not in seen;seen.add((marker,taxon))
        passed=[]
        for label in ['profile','mafft']:
            expected=computed[label,marker][taxon]
            assert tuple(int(r[label+'_'+c]) for c in ['columns','minimum','observed'])==expected
            passed.append(expected[2]>=expected[1])
        status={(True,True):'both',(True,False):'profile_only',(False,True):'mafft_only',(False,False):'neither'}[tuple(passed)]
        assert r['disposition']==status;counts[status]+=1;groups[marker][taxon]=status
    assert len(rows)==len(seen)==65750 and len(groups)==125 and all(len(g)==526 for g in groups.values())
    summaries=read(root/'marker_summary.tsv');assert len(summaries)==125
    for r in summaries:
        marker=r['marker'];g=groups[marker];n=Counter(g.values())
        for k in ['both','profile_only','mafft_only','neither']:assert int(r[k])==n[k]
        assert int(r['profile_taxa'])==n['both']+n['profile_only'] and int(r['mafft_taxa'])==n['both']+n['mafft_only']
        for label in ['profile','mafft']:assert int(r[label+'_columns'])==computed[label,marker][taxa[0]][0]
        for column,status in [('common_taxa','both'),('profile_only_taxa','profile_only'),('mafft_only_taxa','mafft_only')]:assert r[column]==';'.join(sorted(t for t,s in g.items() if s==status))
    assert dict(counts)==receipt['dispositions']
    assert receipt['markers_with_different_taxon_sets']==sum(any(s in ['profile_only','mafft_only'] for s in g.values()) for g in groups.values())
    proof=dict(status='passed_full_marker_alignment_coverage_readback',marker_taxon_cells=65750,markers=125,dispositions=dict(counts),source_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(__file__),scope='Independent FASTA parser and NumPy indexed residue counts reconstruct all 65750 coverage cells, membership dispositions, exact common/differential taxon sets and marker totals. Does not validate unfinished MAFFT inference outputs.')
    Path('metadata/marker_alignment_coverage_readback_20260927.json').write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof,indent=2))


if __name__=='__main__':main()
