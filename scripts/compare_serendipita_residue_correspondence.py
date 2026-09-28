#!/usr/bin/env python3
"""Compare exact protein-position pairs across both alignments for every marker."""
import collections,csv,gzip,hashlib,itertools,json
from pathlib import Path
from Bio import AlignIO,SeqIO


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    base=Path('results/phylogeny');taxa=['F109899','F2600213','F2600215','F2600216','F65672']
    pins={};sources={};canonical=set('ARNDCQEGHILKMFPSTWYV')
    for method,folder in [('profile','profile-alignments-full-v1'),('mafft','alignments-full-v1')]:
        root=base/(method+'-matrix-50-v1');receipt=json.loads((root/'receipt.json').read_text())
        for name in ['matrix.faa','site_mapping.tsv']:
            assert sha(root/name)==receipt['artifacts'][name];pins[str(root/name)]=sha(root/name)
        pins[str(root/'receipt.json')]=sha(root/'receipt.json')
        matrix={r.id:str(r.seq) for r in SeqIO.parse(root/'matrix.faa','fasta') if r.id in taxa}
        sites=collections.defaultdict(list)
        with (root/'site_mapping.tsv').open() as h:
            for r in csv.DictReader(h,delimiter='\t'):sites[r['marker']].append((int(r['matrix_column_1based'])-1,int(r['alignment_column_1based'])-1))
        alignroot=base/folder;arp=alignroot/'receipt.json';pins[str(arp)]=sha(arp)
        alignments={r['marker']:r for r in json.loads(arp.read_text())['alignments']}
        sources[method]=(matrix,sites,alignroot,alignments)
    assert set(sources['profile'][1])==set(sources['mafft'][1]) and len(sources['profile'][1])==125
    out=base/'serendipita-residue-correspondence-20260927-v1';out.mkdir(exist_ok=False)
    rows=[];tuple_count=0
    with gzip.open(out/'five_taxon_position_tuples.jsonl.gz','wt') as tuples:
        for marker in sorted(sources['profile'][1]):
            rawpath=base/'markers-full-v1/unaligned'/(marker+'.faa');raw={r.id:str(r.seq).upper() for r in SeqIO.parse(rawpath,'fasta') if r.id in taxa};pins[str(rawpath)]=sha(rawpath)
            positions={}
            for method,(matrix,sites,alignroot,alignments) in sources.items():
                ar=alignments[marker];assert sha(rawpath)==ar['input_sha256']
                ap=alignroot/(marker+'.faa');assert sha(ap)==ar['output_sha256'];pins[str(ap)]=sha(ap)
                compact={r.id:str(r.seq).upper() for r in SeqIO.parse(ap,'fasta') if r.id in taxa}
                if method=='profile':
                    sp=alignroot/(marker+'.sto');assert sha(sp)==ar['stockholm_sha256'];pins[str(sp)]=sha(sp)
                    alignment=AlignIO.read(sp,'stockholm')
                    rf=alignment.column_annotations['reference_annotation'];columns=[i for i,s in enumerate(rf) if s not in '.- ']
                    assert [i+1 for i in columns]==ar['retained_stockholm_columns_1based']
                    full={r.id:str(r.seq).upper().replace('.','-') for r in alignment if r.id in taxa}
                    assert all(''.join(full[t][i] for i in columns)==compact[t] for t in compact)
                else:full=compact;columns=list(range(len(next(iter(full.values()))))) if full else []
                maps={}
                for t,s in full.items():
                    assert s.replace('-','')==raw[t]
                    n=0;lookup=[]
                    for aa in s:
                        if aa!='-':n+=1
                        lookup.append(n if aa!='-' else None)
                    assert n==len(raw[t]);maps[t]=lookup
                selected=[]
                for mcol,acol in sites[marker]:
                    letters=[matrix[t][mcol] for t in taxa]
                    for t,aa in zip(taxa,letters):
                        original=compact[t][acol] if t in compact else '-'
                        normalized=original if original in canonical or original=='-' else 'X'
                        assert aa==normalized
                    if not all(aa in canonical for aa in letters):continue
                    signature=tuple(maps[t][columns[acol]] for t in taxa)
                    assert all(raw[t][pos-1]==aa for t,pos,aa in zip(taxa,signature,letters))
                    selected.append(signature)
                    tuples.write(json.dumps(dict(method=method,marker=marker,matrix_column=mcol+1,taxa=taxa,protein_positions=signature))+'\n');tuple_count+=1
                assert len(selected)==len(set(selected));positions[method]=selected
            for ia,ib in itertools.combinations(range(5),2):
                a,b=taxa[ia],taxa[ib]
                pairs={method:{(s[ia],s[ib]) for s in sigs} for method,sigs in positions.items()}
                common=pairs['profile']&pairs['mafft']
                for category,selected in [('both',common),('profile_only',pairs['profile']-common),('mafft_only',pairs['mafft']-common)]:
                    differing=sum(raw[a][i-1]!=raw[b][j-1] for i,j in selected)
                    rows.append(dict(marker=marker,taxon_a=a,taxon_b=b,correspondence=category,residue_pairs=len(selected),differing_pairs=differing))
    with (out/'pair_correspondence_summary.tsv').open('w') as h:
        w=csv.DictWriter(h,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    prior=base/'serendipita-shared-marker-sites-20260927-v1/marker_comparisons.tsv'
    pr=json.loads((prior.parent/'receipt.json').read_text());assert sha(prior)==pr['artifacts'][prior.name];pins[str(prior)]=sha(prior)
    index=collections.defaultdict(dict)
    for r in rows:index[r['marker'],r['taxon_a'],r['taxon_b']][r['correspondence']]=r
    with prior.open() as h:
        for r in csv.DictReader(h,delimiter='\t'):
            if r['site_policy']!='all_five_canonical':continue
            joined=index[r['marker'],r['taxon_a'],r['taxon_b']];subset=[joined['both'],joined[r['method']+'_only']]
            assert sum(x['residue_pairs'] for x in subset)==int(r['comparable_columns'])
            assert sum(x['differing_pairs'] for x in subset)==int(r['differing_columns'])
    result=dict(status='complete_all_serendipita_exact_residue_correspondence',markers=125,pair_category_rows=len(rows),five_taxon_tuples=tuple_count,pins=pins,script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in out.iterdir()},scope='Exact full-Stockholm and MAFFT source residues traced through profile match states and matrix masks. Both-method pairs share identical original protein coordinates. Method-only categories conflate site retention and pairing changes; no alignment correctness, orthology or species boundary is established. All markers and zero categories retained.')
    rp=out/'receipt.json';rp.write_text(json.dumps(result,indent=2)+'\n');result.update(completed_receipt_path=str(rp),completed_receipt_sha256=sha(rp))
    Path('metadata/serendipita_residue_correspondence_completed_20260927.json').write_text(json.dumps(result,indent=2)+'\n')
    print('rows',len(rows),'tuples',tuple_count)
    print([r for r in rows if r['marker']=='5004391at2759' and r['taxon_a']=='F2600215' and r['taxon_b']=='F2600216'])


if __name__=='__main__':main()
