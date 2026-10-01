#!/usr/bin/env python3
"""Enumerate every fixed graph node against physical comparison identities before launch.

This census checks whole source files, endpoint roles/versions, catalog hashes,
lengths and physical comparison membership. Prior completed proof locators are
bound; raw structures/native optimization/coverage decisions are not re-audited
here. Full production source gates remain required.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def run(plan_path):
    plan_path=Path(plan_path);plan=json.loads(plan_path.read_text());bindings=dict(plan['pins']);bind(bindings,plan_path);verify(bindings)
    input_closure=json.loads(Path(plan['input_completion']).read_text());assert input_closure['status']=='complete_verified_full_expanded_background_native_input_handoff' and input_closure['active_models']==305434 and input_closure['full_pairs']==146172 and input_closure['exact_process_journals_checked']==2
    bind(bindings,input_closure['full_hash_archive'],input_closure['full_hash_archive_sha256']);bind(bindings,input_closure['producer_receipt'],input_closure['producer_receipt_sha256'])
    source=Path(input_closure['producer_receipt']).parent;ir=json.loads(Path(input_closure['producer_receipt']).read_text())
    for name in ['active_models.jsonl.gz','full_background_work_partition.tsv']:bind(bindings,source/name,ir['artifacts'][name])
    model_path=source/'active_models.jsonl.gz';models={}
    with gzip.open(model_path,'rt') as handle:
        for line in handle:
            row=json.loads(line);key=row['model_id'],row['version'];assert key not in models;models[key]=row
    assert len(models)==305434
    target_plan=json.loads(Path(plan['target_plan']).read_text());target_root=Path(target_plan['output']);target_receipt=json.loads((target_root/'receipt.json').read_text());target_closure=json.loads(Path(plan['target_completion']).read_text())
    assert target_closure['status']=='complete_verified_full_expanded_pair_event_and_taxon_coverage_screens' and len(target_closure['services'])==2
    assert sha(target_root/'receipt.json')==target_closure['source_hashes'][str(target_root/'receipt.json')]
    target_table=target_root/'pair_mask_coverage.tsv';bind(bindings,target_table,target_receipt['artifacts'][target_table.name])
    target_pairs={};seen_masks=set()
    with target_table.open() as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            key=row['pair_key'];ends=[(row['model_'+s],int(row['version_'+s])) for s in ['a','b']];values=dict(zip(ends,[int(row['length_'+s]) for s in ['a','b']]))
            assert (key,row['mask']) not in seen_masks;seen_masks.add((key,row['mask']));assert row['mask'] in ['full','plddt70']
            assert key not in target_pairs or target_pairs[key]==values;target_pairs[key]=values
    assert len(target_pairs)==134812 and len(seen_masks)==269624
    background_pairs={}
    with (source/'full_background_work_partition.tsv').open() as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            key=row['pair_key'];ends=[(row['model_'+s],int(row['version_'+s])) for s in ['a','b']];assert key not in background_pairs and all(end in models for end in ends);background_pairs[key]={end:models[end]['length'] for end in ends}
    assert len(background_pairs)==146172
    graph=Path(plan['graph']);gr=json.loads((graph/'receipt.json').read_text());matched=json.loads(Path(plan['matching_completion']).read_text());assert matched['status']=='complete_verified_expanded_background_fixed_matching' and len(matched['services'])==2
    assert sha(graph/'receipt.json')==matched['source_hashes'][str(graph/'receipt.json')];bind(bindings,graph/'receipt.json')
    counts=Counter();used={'target':set(),'background':set()};float_mismatches=Counter();maximum_float_difference=Counter()
    for kind,name,physical,expected in [('target','target_nodes.jsonl',target_pairs,283409),('background','background_nodes.jsonl',background_pairs,318037)]:
        bind(bindings,graph/name,gr['artifacts'][name]);seen=set()
        with (graph/name).open() as handle:
            for index,line in enumerate(handle,1):
                n=json.loads(line);assert n['node_id'] not in seen;seen.add(n['node_id']);ends=[(n['model_id_'+s],n['version_'+s]) for s in ['a','b']]
                assert n['pair_key']==hashlib.sha256(json.dumps(sorted(ends),separators=(',',':')).encode()).hexdigest() and bool(n['same_model'])==(ends[0]==ends[1])
                prefix=kind+':'+n['guide']+':'
                if ends[0]==ends[1]:counts[prefix+'identical_model_no_alignment']+=1;continue
                assert n['pair_key'] in physical and set(physical[n['pair_key']])==set(ends);used[kind].add(n['pair_key']);counts[prefix+'distinct_model_pair']+=1
                counts[prefix+'reversed_model_roles_vs_canonical']+=ends!=sorted(ends)
                for side,end in zip(['a','b'],ends):
                    assert physical[n['pair_key']][end]==n['length_'+side] and end in models;current=models[end]
                    assert current['length']==n['length_'+side] and current['sha256']==n['sha256_'+side] and current['sequence_sha256']==n['sequence_sha256_'+side]
                    for field in ['mean_ca_plddt','fraction_ca_plddt_below50']:
                        a,b=current[field],n[field+'_'+side];assert math.isfinite(a) and math.isfinite(b)
                        if a!=b:float_mismatches[kind+':'+field]+=1;maximum_float_difference[kind+':'+field]=max(maximum_float_difference[kind+':'+field],abs(a-b))
                if index%100000==0:print('Whole matched-node compatibility census',kind,index,'/',expected,flush=True)
        assert len(seen)==expected==gr[kind+'_nodes']
    assert used['target']==set(target_pairs) and used['background']==set(background_pairs)
    verify(bindings)
    result=dict(status='passed_complete_fixed_matching_node_physical_identity_census',plan_sha256=sha(plan_path),target_nodes=283409,background_nodes=318037,target_physical_pairs=134812,background_physical_pairs=146172,source_inventory_models=len(models),counts=dict(counts),used_target_pairs=len(used['target']),used_background_pairs=len(used['background']),confidence_descriptor_mismatch_endpoint_occurrences=dict(float_mismatches),maximum_confidence_descriptor_absolute_differences=dict(maximum_float_difference),source_hashes=bindings,scientific_eligibility=False,scope=__doc__)
    with Path(plan['output']).open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2));return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);run(parser.parse_args().plan)
