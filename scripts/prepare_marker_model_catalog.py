#!/usr/bin/env python3
"""Freeze source-specific marker models for confidence prefetch before alignment mapping ends."""
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from compare_marker_structures import ROOT, sha
from map_marker_structures import DEFAULT_PROVIDER, DEFAULT_TOOL, source_candidates


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inventory',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('Use a new immutable catalog')
    raw=args.inventory.read_bytes()
    latest={r['uniprot_accession']:r for r in map(json.loads,raw.decode().splitlines())}
    by_sequence=defaultdict(list)
    for row in latest.values():
        if row['status']=='verified':
            for model in row['models']:
                by_sequence[model['sequence_sha256']].append(dict(model,uniprot_accession=row['uniprot_accession']))
    source=ROOT/'results/phylogeny/markers-full-v1/protein_mapping.tsv'
    rows=list(csv.DictReader(source.open(),delimiter='\t'))
    selected,links={},[]
    for row in rows:
        candidates=source_candidates(by_sequence.get(row['sequence_sha256'],[]),DEFAULT_PROVIDER,DEFAULT_TOOL)
        if not candidates:continue
        model=max(candidates,key=lambda m:(m['mean_ca_plddt'],m['version'],m['model_id']))
        if model['path'] not in selected:
            path=ROOT/model['path']
            if sha(path)!=model['sha256']:raise ValueError('Changed validated coordinate artifact')
            selected[model['path']]=model
        elif selected[model['path']]['sha256']!=model['sha256']:
            raise ValueError('Conflicting model identities')
        links.append({'marker':row['marker'],'taxon_id':row['taxon_id'],'protein_id':row['protein_id'],
                      'sequence_sha256':row['sequence_sha256'],'model_id':model['model_id'],'version':model['version'],'model_path':model['path']})
    if len({(m['model_id'],m['version']) for m in selected.values()})!=len(selected):
        raise ValueError('Duplicate model/version under distinct paths')
    args.output.mkdir(parents=True)
    (args.output/'model_provenance.json').write_text(json.dumps(list(selected.values()),indent=2)+'\n')
    with (args.output/'marker_model_links.tsv').open('w') as out:
        writer=csv.DictWriter(out,list(links[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(links)
    receipt={'status':'complete_source_specific_model_catalog','inventory_path':str(args.inventory),
        'inventory_sha256':hashlib.sha256(raw).hexdigest(),'marker_mapping_sha256':sha(source),
        'source_policy':{'provider':DEFAULT_PROVIDER,'tool':DEFAULT_TOOL},
        'marker_proteins_screened':len(rows),'marker_proteins_linked':len(links),'distinct_models':len(selected),
        'distinct_taxa_linked':len({r['taxon_id'] for r in links}),
        'sum_squared_model_lengths':sum(m['length']**2 for m in selected.values()),
        'script_sha256':sha(Path(__file__)),'selection_source_sha256':sha(Path(__file__).with_name('map_marker_structures.py')),
        'interpretation':'Prefetch catalog derived from the same frozen inventory and source-selection policy as the running residue mapper. Cached coordinate hashes verified; this does not complete alignment/coordinate mapping or confidence audit. Final mapping model identities and hashes must match before confidence artifacts are bound to that mapping.',
        'artifacts':{p.name:sha(p) for p in args.output.iterdir()}}
    (args.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
