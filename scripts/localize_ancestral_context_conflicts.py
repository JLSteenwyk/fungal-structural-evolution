#!/usr/bin/env python3
"""Retain coordinate and model dependencies of opposing ancestral calls."""
import csv,gzip,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path


def main():
    pins={}
    def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
    root=Path('results/ancestral/refined-whole-domain-probability-comparison-20260927-v1')
    receipt=json.loads((root/'receipt.json').read_text());pins[str(root/'receipt.json')]=sha(root/'receipt.json')
    table=root/'matched_site_comparisons.tsv.gz';assert sha(table)==receipt['artifacts'][table.name]
    readback=Path('metadata/refined_whole_domain_probability_readback_completed_20260927.json')
    audit=json.loads(readback.read_text());assert audit['source_receipt_sha256']==sha(root/'receipt.json');pins[str(readback)]=sha(readback)
    key=lambda r:(r['family'],r['source_node'],r['coordinate_signature_sha256'])
    conflicts=defaultdict(list)
    with gzip.open(table,'rt') as h:
        for r in csv.DictReader(h,delimiter='\t'):
            if r['opposing_at_least_090']=='True':conflicts[key(r)].append(r)
    contexts=defaultdict(list)
    with gzip.open(table,'rt') as h:
        for r in csv.DictReader(h,delimiter='\t'):
            if key(r) in conflicts:contexts[key(r)].append(r)
    cp=Path('results/ancestral/case-domain-alignment-readback-20260927-v1/protein_position_maps.tsv.gz')
    pins[str(cp)]=sha(cp)
    needed={(r['family'],r['boundary'],r['domain_method'],r['domain_column']) for rs in conflicts.values() for r in rs}
    coordinates=defaultdict(list)
    with gzip.open(cp,'rt') as h:
        for r in csv.DictReader(h,delimiter='\t'):
            k=(r['family'],r['boundary'],r['method'],r['column'])
            if k in needed:coordinates[k].append((r['gene'],int(r['protein_position']),r['residue']))
    results=[]
    for k,rs in sorted(conflicts.items()):
        coordinate_sets=set()
        for r in rs:
            coords=tuple(sorted(coordinates[r['family'],r['boundary'],r['domain_method'],r['domain_column']]))
            signature=hashlib.sha256(json.dumps([(g,p) for g,p,a in coords],separators=(',',':')).encode()).hexdigest()
            assert signature==r['coordinate_signature_sha256'];coordinate_sets.add(coords)
        assert len(coordinate_sets)==1
        coords=next(iter(coordinate_sets))
        results.append(dict(family=k[0],source_node=k[1],coordinate_signature_sha256=k[2],opposing_comparisons=len(rs),all_available_comparisons_at_exact_signature=len(contexts[k]),opposing_state_pairs=dict(Counter(r['whole_map']+'->'+r['domain_map'] for r in rs)),models_with_conflict=sorted({r['model'] for r in rs}),boundaries_with_conflict=sorted({r['boundary'] for r in rs}),observed_residue_counts=dict(Counter(a for g,p,a in coords)),extant_coordinates=[dict(gene=g,protein_position=p,residue=a) for g,p,a in coords],all_comparisons=contexts[k]))
    assert sum(r['opposing_comparisons'] for r in results)==receipt['opposing_at_least_090']==104
    output=Path('results/ancestral/ancestral-context-conflict-localization-20260927-v1');output.mkdir(exist_ok=False)
    artifact=output/'conflict_dossiers.json';artifact.write_text(json.dumps(results,indent=2)+'\n')
    result=dict(status='complete_exact_coordinate_conflict_localization',opposing_comparisons=104,node_coordinate_groups=len(results),coordinate_signatures=len({r['coordinate_signature_sha256'] for r in results}),families=sorted({r['family'] for r in results}),pins=pins,script_sha256=sha(__file__),artifacts={artifact.name:sha(artifact)},scope='Grouping exact extant-coordinate signatures and original source nodes removes repeated model/bound/method labels from descriptive counts. Different nodes and distinct signatures can remain dependent; these are not unique biological substitutions, causal effects, or experimentally validated functional sites.')
    (output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    result.update(completed_receipt_path=str(output/'receipt.json'),completed_receipt_sha256=sha(output/'receipt.json'))
    Path('metadata/ancestral_context_conflict_localization_completed_20260927.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
