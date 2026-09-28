#!/usr/bin/env python3
"""Compare all Serendipita pairs per marker on pairwise and five-taxon shared sites."""
import csv
import hashlib
import itertools
import json
import re
from pathlib import Path
from Bio import SeqIO


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, rows):
    with path.open('w') as stream:
        writer = csv.DictWriter(stream, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    oldroot = Path('results/phylogeny/serendipita-identity-context-20260927-v1')
    rp = oldroot / 'receipt.json'
    oldreceipt = json.loads(rp.read_text())
    p = oldroot / 'all_pair_observed_differences.tsv'
    assert sha(p) == oldreceipt['artifacts'][p.name]
    pins = {str(x):sha(x) for x in [rp,p]}
    with p.open() as stream:
        old = {(r['method'],r['taxon_a'],r['taxon_b']):r for r in csv.DictReader(stream,delimiter='\t')}
    taxa = sorted({k[1] for k in old}|{k[2] for k in old})
    assert len(taxa)==5
    canonical = set('ARNDCQEGHILKMFPSTWYV')
    rows=[]
    summaries=[]
    for method in ['profile','mafft']:
        root=Path('results/phylogeny')/(method+'-matrix-50-v1')
        receipt=json.loads((root/'receipt.json').read_text())
        pins[str(root/'receipt.json')]=sha(root/'receipt.json')
        for name in ['matrix.faa','partitions.nex']:
            assert sha(root/name)==receipt['artifacts'][name]
            pins[str(root/name)]=sha(root/name)
        sequences={r.id:str(r.seq) for r in SeqIO.parse(root/'matrix.faa','fasta') if r.id in taxa}
        assert set(sequences)==set(taxa)
        partitions=[(marker,int(start)-1,int(end)) for marker,start,end in
                    re.findall(r'charset marker_(\w+) = (\d+)-(\d+);',(root/'partitions.nex').read_text())]
        assert len(partitions)==receipt['markers']
        assert [i for _,start,end in partitions for i in range(start,end)]==list(range(receipt['columns']))
        shared={i for i in range(receipt['columns']) if all(sequences[t][i] in canonical for t in taxa)}
        for a,b in itertools.combinations(taxa,2):
            valid={i for i in range(receipt['columns']) if sequences[a][i] in canonical and sequences[b][i] in canonical}
            assert shared<=valid
            for policy,positions in [('pairwise_canonical',valid),('all_five_canonical',shared)]:
                marker_rows=[]
                for marker,start,end in partitions:
                    sites=positions&set(range(start,end))
                    differing=sum(sequences[a][i]!=sequences[b][i] for i in sites)
                    marker_rows.append(dict(method=method,taxon_a=a,taxon_b=b,site_policy=policy,
                        marker=marker,marker_columns=end-start,comparable_columns=len(sites),
                        differing_columns=differing,excluded_columns=end-start-len(sites),
                        observed_difference_fraction=differing/len(sites) if sites else 'NA'))
                n=sum(r['comparable_columns'] for r in marker_rows)
                d=sum(r['differing_columns'] for r in marker_rows)
                assert n==len(positions)
                assert d==sum(x!=y for i,(x,y) in enumerate(zip(sequences[a],sequences[b])) if i in positions)
                if policy=='pairwise_canonical':
                    assert (n,d)==(int(old[method,a,b]['paired_canonical_columns']),int(old[method,a,b]['differing_columns']))
                summaries.append(dict(method=method,taxon_a=a,taxon_b=b,site_policy=policy,
                    comparable_columns=n,differing_columns=d,observed_difference_fraction=d/n if n else 'NA',
                    markers=len(marker_rows),markers_without_comparable_sites=sum(r['comparable_columns']==0 for r in marker_rows),
                    markers_with_differences=sum(r['differing_columns']>0 for r in marker_rows)))
                rows.extend(marker_rows)
    out=Path('results/phylogeny/serendipita-shared-marker-sites-20260927-v1')
    out.mkdir(exist_ok=False)
    write(out/'marker_comparisons.tsv',rows)
    write(out/'pair_summary.tsv',summaries)
    result=dict(status='complete_all_serendipita_shared_site_marker_comparisons',pairs=20,
        site_policies=2,summary_rows=len(summaries),marker_rows=len(rows),pins=pins,
        script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in out.iterdir()},
        verification='Partition coverage is exact and disjoint; all per-marker sums match independent full-string counts; pairwise results reproduce the prior frozen table.',
        scope='Within each alignment, all-five policy compares identical positions for every pair. Different alignment methods still use different positions. Missing/zero-marker dispositions retained. Uncorrected descriptive differences do not delimit species or validate orthology, and shared-site restriction changes the sampled residues.')
    rp=out/'receipt.json'
    rp.write_text(json.dumps(result,indent=2)+'\n')
    result.update(completed_receipt_path=str(rp),completed_receipt_sha256=sha(rp))
    Path('metadata/serendipita_shared_marker_sites_completed_20260927.json').write_text(json.dumps(result,indent=2)+'\n')
    for r in summaries:
        if r['site_policy']=='all_five_canonical' and r['taxon_a'].startswith('F260') and r['taxon_b'].startswith('F260'):print(r)


if __name__=='__main__':
    main()
