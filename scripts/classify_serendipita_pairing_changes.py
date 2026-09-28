#!/usr/bin/env python3
"""Separate alternate residue pairing from absent endpoint coverage in all comparisons."""
import collections,csv,gzip,hashlib,itertools,json
from pathlib import Path


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    root=Path('results/phylogeny/serendipita-residue-correspondence-20260927-v1')
    rp=root/'receipt.json';receipt=json.loads(rp.read_text());pins={str(rp):sha(rp)}
    for name,h in receipt['artifacts'].items():assert sha(root/name)==h;pins[str(root/name)]=h
    positions=collections.defaultdict(list)
    with gzip.open(root/'five_taxon_position_tuples.jsonl.gz','rt') as h:
        for r in map(json.loads,h):positions[r['marker'],r['method']].append(tuple(r['protein_positions']));taxa=r['taxa']
    with (root/'pair_correspondence_summary.tsv').open() as h:previous=list(csv.DictReader(h,delimiter='\t'))
    markers=sorted({r['marker'] for r in previous});assert len(markers)==125
    old={(r['marker'],r['taxon_a'],r['taxon_b'],r['correspondence']):int(r['residue_pairs']) for r in previous}
    rows=[]
    for marker in markers:
        for ia,ib in itertools.combinations(range(5),2):
            a,b=taxa[ia],taxa[ib]
            pairs={method:{(s[ia],s[ib]) for s in positions[marker,method]} for method in ['profile','mafft']}
            for method,other in [('profile','mafft'),('mafft','profile')]:
                endpoints=[{p[i] for p in pairs[other]} for i in [0,1]]
                categories=collections.Counter()
                for x,y in pairs[method]-pairs[other]:
                    status=('both_endpoints_retained_different_pairing' if x in endpoints[0] and y in endpoints[1] else
                            'only_a_endpoint_retained' if x in endpoints[0] else
                            'only_b_endpoint_retained' if y in endpoints[1] else 'neither_endpoint_retained')
                    categories[status]+=1
                assert sum(categories.values())==old[marker,a,b,method+'_only']
                rows.append(dict(marker=marker,taxon_a=a,taxon_b=b,method=method,
                    exact_pairs_shared=len(pairs[method]&pairs[other]),method_only_pairs=sum(categories.values()),
                    **{k:categories[k] for k in ['both_endpoints_retained_different_pairing','only_a_endpoint_retained','only_b_endpoint_retained','neither_endpoint_retained']}))
    out=Path('results/phylogeny/serendipita-pairing-changes-20260927-v1');out.mkdir(exist_ok=False)
    p=out/'pairing_dispositions.tsv'
    with p.open('w') as h:w=csv.DictWriter(h,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    result=dict(status='complete_all_serendipita_pairing_change_dispositions',rows=len(rows),markers=125,pins=pins,script_sha256=sha(__file__),artifacts={p.name:sha(p)},scope='Within all-five-canonical site policies, retained endpoints paired to different partners are distinguished from absent endpoints. Endpoint absence conflates masks, profile-state selection and other taxa coverage. No preferred alignment, homology validation, marker exclusion or species delimitation.')
    rp=out/'receipt.json';rp.write_text(json.dumps(result,indent=2)+'\n');result.update(completed_receipt_path=str(rp),completed_receipt_sha256=sha(rp))
    Path('metadata/serendipita_pairing_changes_completed_20260927.json').write_text(json.dumps(result,indent=2)+'\n')
    print([r for r in rows if r['marker']=='5004391at2759' and r['taxon_a']=='F2600215' and r['taxon_b']=='F2600216'])


if __name__=='__main__':main()
