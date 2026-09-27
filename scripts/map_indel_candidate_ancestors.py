#!/usr/bin/env python3
"""Map all intended ancestor levels onto gap-model vertices by tip partitions."""
import csv
import hashlib
import json
from collections import Counter,defaultdict
from pathlib import Path
from Bio import Phylo
from map_ancestral_fitted_nodes import signatures
from prepare_case_ancestral_neighborhoods import sha


def main():
    pp=Path('metadata/stable_indel_refinement_plan_20260927.json')
    plan=json.loads(pp.read_text());root=Path('results/ancestral/case-local-trees-20260927-v1')
    receipt=json.loads((root/'receipt.json').read_text())
    mp=root/'ancestral_node_mapping.tsv'
    assert sha(mp)==receipt['artifacts'][mp.name]
    pins={str(pp):sha(pp),str(mp):sha(mp),str(root/'receipt.json'):sha(root/'receipt.json')}
    with mp.open() as handle:
        original=list(csv.DictReader(handle,delimiter='\t'))
    source={}
    for row in original:
        key=row['guide'],row['family'],row['dataset']
        if key not in source:
            path=root/('-'.join(key)+'.nwk')
            assert sha(path)==receipt['artifacts'][path.name]
            pins[str(path)]=sha(path)
            tree=Phylo.read(path,'newick');sigs=signatures(tree)
            source[key]=(tree,sigs)
    results=[]
    for job in plan['jobs']:
        tree_path=Path(job['tree']);assert sha(tree_path)==job['pins'][str(tree_path)]
        pins[str(tree_path)]=sha(tree_path)
        target=Phylo.read(tree_path,'newick');target_sigs=signatures(target)
        named={node:f'indel_internal_{i}' for i,node in enumerate(target.find_clades()) if not node.is_terminal()}
        by=defaultdict(list)
        for node,sig in target_sigs.items():
            by[sig].append(node)
        family,boundary=job['job_id'].split('-')[:2]
        dataset='whole' if boundary=='whole' else 'domain'
        for row in original:
            if row['family']!=family or row['dataset']!=dataset:
                continue
            tree,sigs=source[row['guide'],family,dataset]
            assert {n.name for n in tree.get_terminals()}=={n.name for n in target.get_terminals()}
            names=json.loads(row['output_nodes_json']);assert len(names)==1
            nodes=[n for n in sigs if n.name==names[0]];assert len(nodes)==1
            node=nodes[0];sig=sigs[node];matches=by.get(sig,[])
            if len(sig)==2:
                assert node is tree.root and not matches
                status='degree_two_root_position_not_identified'
            else:
                assert len(sig)>=3 and len(matches)==1,(job['job_id'],row['level'])
                status='matched_unrooted_vertex'
            results.append(dict(job_id=job['job_id'],family=family,boundary=boundary,
                guide=row['guide'],level=int(row['level']),source_node=node.name,status=status,
                posterior_node=named[matches[0]] if matches else '',
                incident_tip_partition_sha256=hashlib.sha256(json.dumps(sig,separators=(',',':')).encode()).hexdigest(),
                character_count=job['character_count'],
                data_disposition='nonempty' if job['character_count'] else 'no_coded_characters',
                tree_sha256=sha(tree_path)))
    assert len(results)==len(plan['jobs'])*8
    out=Path('results/ancestral/indel-candidate-node-mapping-20260927-v1');out.mkdir(parents=True,exist_ok=False)
    with (out/'candidate_nodes.tsv').open('w') as handle:
        writer=csv.DictWriter(handle,list(results[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(results)
    result=dict(status='all_312_gap_model_candidate_nodes_mapped',jobs=len(plan['jobs']),
        mapping_rows=len(results),status_counts=dict(Counter(r['status'] for r in results)),
        unique_matched_job_levels=len({(r['job_id'],r['level']) for r in results if r['status']=='matched_unrooted_vertex'}),
        pins=pins,script_sha256=sha(__file__),artifacts={'candidate_nodes.tsv':sha(out/'candidate_nodes.tsv')},
        scope='All four candidate levels and both guide aliases mapped for whole/domain gap trees. Degree-two local roots remain unidentified; empty matrices retained. Maps refer to deterministic production node labels through verified tip partitions, not label equality. No new biological roots, posterior validation or final ancestral sequence claim.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result['mapping_rows'],result['status_counts'],result['unique_matched_job_levels'])


if __name__=='__main__':
    main()
