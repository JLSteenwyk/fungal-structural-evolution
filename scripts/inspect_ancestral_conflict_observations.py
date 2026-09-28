#!/usr/bin/env python3
"""Partition all observed whole-column characters for localized conflicts."""
import csv,hashlib,json
from collections import Counter
from pathlib import Path
from Bio import SeqIO


def main():
    pins={}
    def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
    def bind(p):p=Path(p);pins[str(p)]=sha(p);return p
    root=Path('results/ancestral/ancestral-context-conflict-localization-20260927-v1')
    receipt=json.loads(bind(root/'receipt.json').read_text());source=bind(root/'conflict_dossiers.json');assert sha(source)==receipt['artifacts'][source.name]
    dossiers=json.loads(source.read_text())
    dp=bind('results/ancestral/case-domain-sequences-20260927-v1/protein_dispositions.tsv')
    with dp.open() as h:dispositions={(r['family'],r['boundary'],r['gene']):r for r in csv.DictReader(h,delimiter='\t')}
    mp=bind('results/ancestral/case-local-trees-20260927-v1/ancestral_node_mapping.tsv')
    with mp.open() as h:nodes={(r['family'],r['source_node']):set(json.loads(r['retained_set_json'])) for r in csv.DictReader(h,delimiter='\t') if r['guide']=='profile' and r['dataset']=='whole'}
    wr=json.loads(bind('results/ancestral/case-alignments-20260927-v1/receipt.json').read_text())
    hashes={(r['family'],r['method']):r['artifacts']['alignment.faa'] for r in wr['results']}
    cache={};results=[]
    for dossier in dossiers:
        fam=dossier['family'];node=dossier['source_node'];expected={(r['gene'],r['protein_position'],r['residue']) for r in dossier['extant_coordinates']}
        contexts=sorted({(r['boundary'],r['whole_method'],int(r['whole_column'])) for r in dossier['all_comparisons']})
        for boundary,method,col in contexts:
            key=fam,method
            if key not in cache:
                p=bind(Path('results/ancestral/case-alignments-20260927-v1')/(fam+'-'+method)/'alignment.faa');assert sha(p)==hashes[key]
                cache[key]={r.id:str(r.seq) for r in SeqIO.parse(p,'fasta')}
            observations=[];retained=set()
            for gene,seq in cache[key].items():
                aa=seq[col-1];position=len(seq[:col].replace('-','')) if aa!='-' else None
                d=dispositions[fam,boundary,gene]
                if aa=='-':category='alignment_gap'
                elif d['status']!='extracted_single_target_hit':category='protein_excluded_from_domain_analysis'
                elif int(d['start'])<=position<=int(d['end']):category='retained_domain_coordinate';retained.add((gene,position,aa))
                else:category='outside_focal_domain_interval'
                observations.append(dict(gene=gene,residue=aa,protein_position=position,category=category,domain_disposition=d['status'],candidate_descendant=gene in nodes[fam,node]))
            assert retained==expected
            counts={}
            for category in sorted({r['category'] for r in observations}):
                selected=[r for r in observations if r['category']==category]
                counts[category]=dict(all_characters=dict(Counter(r['residue'] for r in selected)),descendant_characters=dict(Counter(r['residue'] for r in selected if r['candidate_descendant'])))
            results.append(dict(family=fam,source_node=node,boundary=boundary,whole_method=method,whole_column=col,coordinate_signature_sha256=dossier['coordinate_signature_sha256'],whole_tips=len(observations),category_counts=counts,observations=observations))
    output=Path('results/ancestral/ancestral-conflict-observations-20260927-v1');output.mkdir(exist_ok=False)
    p=output/'observed_contexts.json';p.write_text(json.dumps(results,indent=2)+'\n')
    summary=dict(status='complete_observed_context_partition',contexts=len(results),node_coordinate_groups=len(dossiers),observed_character_records=sum(r['whole_tips'] for r in results),pins=pins,script_sha256=sha(__file__),artifacts={p.name:sha(p)},scope='Observed full-column characters partitioned by domain retention and candidate descendants. Exact retained coordinates reproduce conflict signatures. Counts do not weight phylogeny or identify causes of ancestral probability differences; excluded does not mean erroneous annotation.')
    (output/'receipt.json').write_text(json.dumps(summary,indent=2)+'\n');summary.update(completed_receipt_path=str(output/'receipt.json'),completed_receipt_sha256=sha(output/'receipt.json'))
    Path('metadata/ancestral_conflict_observations_completed_20260927.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary))
    for r in results:
        if r['whole_method']=='mafft':print(r['family'],r['source_node'],r['boundary'],json.dumps(r['category_counts']))


if __name__=='__main__':main()
