#!/usr/bin/env python3
"""Preserve marker contributions and selected-protein provenance for every comparison."""
import collections
import csv
import hashlib
import json
from pathlib import Path
from Bio import SeqIO


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read(p):
    with Path(p).open() as h:
        return list(csv.DictReader(h, delimiter='\t'))


def write(p, rows):
    with p.open('w') as h:
        w=csv.DictWriter(h,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)


def main():
    root=Path('results/phylogeny/serendipita-shared-marker-sites-20260927-v1')
    rp=root/'receipt.json';r=json.loads(rp.read_text());p=root/'marker_comparisons.tsv'
    assert sha(p)==r['artifacts'][p.name]
    pins={str(x):sha(x) for x in [rp,p]};rows=read(p)
    taxa={r['taxon_a'] for r in rows}|{r['taxon_b'] for r in rows}
    markers={r['marker'] for r in rows};assert len(taxa)==5 and len(markers)==125
    mr=Path('results/phylogeny/markers-full-v1');mp=mr/'protein_mapping.tsv'
    mapping={(r['taxon_id'],r['marker']):r for r in read(mp) if r['taxon_id'] in taxa}
    pins[str(mp)]=sha(mp);pins[str(mr/'receipt.json')]=sha(mr/'receipt.json')
    cr=Path('results/qc/busco-gene-copies-v1');cp=cr/'marker_gene_copies.tsv';receipt=json.loads((cr/'receipt.json').read_text())
    assert sha(cp)==receipt['artifacts'][cp.name]
    assert sha(mr/'receipt.json')==receipt['marker_extraction_receipt_sha256']
    pins.update({str(x):sha(x) for x in [cp,cr/'receipt.json']})
    copies={(r['taxon_id'],r['marker']):r for r in read(cp) if r['taxon_id'] in taxa}
    sequences={}
    for marker in sorted(markers):
        p=mr/'unaligned'/(marker+'.faa');pins[str(p)]=sha(p)
        for record in SeqIO.parse(p,'fasta'):
            if record.id not in taxa:continue
            key=(record.id,marker);sequence=str(record.seq)
            assert hashlib.sha256(sequence.encode()).hexdigest()==mapping[key]['sequence_sha256']
            assert mapping[key]['protein_id'] in json.loads(copies[key]['protein_ids_json'])
            sequences[key]=sequence
    assert set(sequences)==set(mapping)
    for row in rows:
        for side in ['a','b']:
            key=(row['taxon_'+side],row['marker']);m=mapping.get(key)
            row[side+'_selected_protein']=m['protein_id'] if m else ''
            row[side+'_sequence_sha256']=m['sequence_sha256'] if m else ''
            row[side+'_sequence_length']=len(sequences[key]) if key in sequences else 0
            row[side+'_busco_status']=copies[key]['busco_status']
            row[side+'_raw_copy_count']=copies[key]['raw_hit_proteins']
            if int(row['comparable_columns']):assert m is not None
    groups=collections.defaultdict(list)
    for row in rows:groups[tuple(row[k] for k in ['method','taxon_a','taxon_b','site_policy'])].append(row)
    summaries=[]
    for key,group in sorted(groups.items()):
        ordered=sorted(group,key=lambda r:(-int(r['differing_columns']),r['marker']))
        total=sum(int(r['differing_columns']) for r in ordered)
        summary=dict(zip(['method','taxon_a','taxon_b','site_policy'],key))
        summary.update(total_differences=total,markers_with_data=sum(int(r['comparable_columns'])>0 for r in group),
            identical_markers_with_data=sum(int(r['comparable_columns'])>0 and int(r['differing_columns'])==0 for r in group))
        for count in [1,3,5,10]:
            summary['top_'+str(count)+'_difference_fraction']=sum(int(r['differing_columns']) for r in ordered[:count])/total if total else 'NA'
        cumulative=0
        for rank,row in enumerate(ordered,1):
            cumulative+=int(row['differing_columns']);row['difference_rank']=rank
            row['cumulative_difference_fraction']=cumulative/total if total else 'NA'
        assert cumulative==total;summaries.append(summary)
    out=Path('results/phylogeny/serendipita-marker-provenance-20260927-v1');out.mkdir(exist_ok=False)
    write(out/'marker_provenance.tsv',rows);write(out/'concentration_summary.tsv',summaries)
    result=dict(status='complete_all_serendipita_marker_provenance_and_concentration',marker_comparisons=len(rows),summaries=len(summaries),selected_protein_sequences=len(sequences),pins=pins,script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in out.iterdir()},scope='All markers and zero-data rows retained. Rank is descriptive contribution within each pair/method/site policy, not an exclusion rule. Selected sequences match extraction hashes and raw BUSCO protein calls; this does not prove orthology, annotation quality, contamination absence or species boundaries.')
    rp=out/'receipt.json';rp.write_text(json.dumps(result,indent=2)+'\n');result.update(completed_receipt_path=str(rp),completed_receipt_sha256=sha(rp))
    Path('metadata/serendipita_marker_provenance_completed_20260927.json').write_text(json.dumps(result,indent=2)+'\n')
    for r in summaries:
        if r['taxon_a']=='F2600215' and r['taxon_b']=='F2600216' and r['site_policy']=='all_five_canonical':print(r)


if __name__=='__main__':
    main()
