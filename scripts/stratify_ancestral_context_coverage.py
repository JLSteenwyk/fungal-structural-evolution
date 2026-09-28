#!/usr/bin/env python3
"""Count known descendant residues for every exact matched-site comparison."""
import csv,gzip,hashlib,json
from collections import defaultdict
from pathlib import Path
import numpy as np
from Bio import SeqIO


def main():
    pins={}
    def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
    def read(p):p=Path(p);pins[str(p)]=sha(p);return json.loads(p.read_text())
    source=Path('results/ancestral/refined-whole-domain-probability-comparison-20260927-v1')
    sr=read(source/'receipt.json');audit=read('metadata/refined_whole_domain_probability_readback_completed_20260927.json')
    assert audit['source_receipt_sha256']==sha(source/'receipt.json')
    table=source/'matched_site_comparisons.tsv.gz';assert sha(table)==sr['artifacts'][table.name]
    mp=Path('results/ancestral/case-local-trees-20260927-v1/ancestral_node_mapping.tsv');pins[str(mp)]=sha(mp)
    with mp.open() as h:nodes={(r['family'],r['dataset'],int(r['level'])):set(json.loads(r['retained_set_json'])) for r in csv.DictReader(h,delimiter='\t') if r['guide']=='profile'}
    roots={'whole':Path('results/ancestral/case-alignments-20260927-v1'),'domain':Path('results/ancestral/case-domain-alignments-20260927-v1')}
    manifests={k:read(p/'receipt.json') for k,p in roots.items()}
    hashes={}
    for kind,receipt in manifests.items():
        for r in receipt['results']:hashes[kind,r['family'],r.get('boundary','whole'),r['method']]=r['artifacts']['alignment.faa']
    cache={}
    def counts(kind,family,boundary,method,level):
        key=kind,family,boundary,method,level
        if key not in cache:
            name=family+'-'+(boundary+'-' if kind=='domain' else '')+method
            path=roots[kind]/name/'alignment.faa';assert sha(path)==hashes[kind,family,boundary,method]
            records={r.id:str(r.seq) for r in SeqIO.parse(path,'fasta')};desc=nodes[family,kind,level];assert desc<=set(records)
            matrix=np.array([list(records[g]) for g in sorted(desc)])
            known=np.isin(matrix,list('ARNDCQEGHILKMFPSTWYV')).sum(axis=0)
            unknown=(matrix=='X').sum(axis=0);gaps=(matrix=='-').sum(axis=0)
            assert np.all(known+unknown+gaps==len(desc));cache[key]=(known,unknown,gaps,len(desc))
        return cache[key]
    output=Path('results/ancestral/ancestral-context-coverage-20260927-v1');output.mkdir(exist_ok=False)
    totals=defaultdict(lambda:dict(comparisons=0,map_disagreements=0,opposing_at_least_090=0))
    grouped=defaultdict(lambda:dict(comparisons=0,map_disagreements=0,opposing_at_least_090=0))
    n=0
    with gzip.open(table,'rt') as src,gzip.open(output/'site_coverage.tsv.gz','wt') as dst:
        reader=csv.DictReader(src,delimiter='\t')
        fields=reader.fieldnames+['whole_known_descendants','domain_known_descendants','whole_unknown_descendants','domain_unknown_descendants','whole_gap_descendants','domain_gap_descendants','coverage_class']
        writer=csv.DictWriter(dst,fields,delimiter='\t',lineterminator='\n');writer.writeheader()
        for r in reader:
            level=int(r['level']);wc=int(r['whole_column'])-1;dc=int(r['domain_column'])-1
            w=counts('whole',r['family'],'whole',r['whole_method'],level);d=counts('domain',r['family'],r['boundary'],r['domain_method'],level)
            wk=int(w[0][wc]);dk=int(d[0][dc]);assert wk>=dk
            label='both_have_known_descendants' if dk else ('whole_only_has_known_descendants' if wk else 'neither_has_known_descendants')
            for target in [totals[label],grouped[r['family'],r['boundary'],label]]:
                target['comparisons']+=1;target['map_disagreements']+=r['whole_map']!=r['domain_map'];target['opposing_at_least_090']+=r['opposing_at_least_090']=='True'
            writer.writerow(dict(r,whole_known_descendants=wk,domain_known_descendants=dk,whole_unknown_descendants=int(w[1][wc]),domain_unknown_descendants=int(d[1][dc]),whole_gap_descendants=int(w[2][wc]),domain_gap_descendants=int(d[2][dc]),coverage_class=label));n+=1
    assert n==sr['matched_node_sites']
    for name in ['map_disagreements','opposing_at_least_090']:assert sum(t[name] for t in totals.values())==sr[name]
    rows=[dict(family=f,boundary=b,coverage_class=c,**v) for (f,b,c),v in sorted(grouped.items())]
    with (output/'coverage_summary.tsv').open('w') as h:
        w=csv.DictWriter(h,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    result=dict(status='complete_all_matched_site_descendant_coverage',comparisons=n,coverage_classes=dict(totals),pins=pins,script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in output.iterdir()},scope='All existing exact-matched comparisons retained. Known canonical residues counted separately from X and gaps among mapped descendants. This is observed coverage, not ancestral presence probability, phylogenetic effective sample size, or a correctness filter; outside-descendant tips also inform ancestral marginals. Counts repeat models, bounds, methods and nodes.')
    (output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');result.update(completed_receipt_path=str(output/'receipt.json'),completed_receipt_sha256=sha(output/'receipt.json'))
    Path('metadata/ancestral_context_coverage_completed_20260927.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['pins','artifacts']}))


if __name__=='__main__':main()
