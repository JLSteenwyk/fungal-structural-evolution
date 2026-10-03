"""Full source/runtime rooted-clade and candidate mapping replay."""
from collections import Counter
import csv
import json
from pathlib import Path

from independent_native_alignment_replay_sources import load as native_sources
from independent_native_ancestral_alignment import fasta_records
from independent_native_ancestral_topology import match_trees
from independent_native_alignment_replay import digest
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha

SUMMARY_FIELDS=['full_chains','checked_chains','failed_chains','full_groups','complete_groups','unresolved_groups',
    'intact_chains_in_unresolved_groups','effective_input_groups','node_pairs','nonroot_branch_pairs',
    'candidate_node_pairs','candidate_frame_mappings','assumed_root_candidates','negative_branch_pairs_require_review']


def load(plan,path):
    source,bindings=native_sources(plan,path);mp=Path(plan['mapping_table']);bind(bindings,mp,plan['mapping_table_sha256'])
    selected={}
    with mp.open() as handle:
        reader=csv.DictReader(handle,delimiter='\t')
        for row in reader:
            key=row['guide'],row['family'],row['dataset'];selected.setdefault(key,[]).append(row)
    for rows in selected.values():
        assert len(rows)==4 and {int(r['level']) for r in rows}==set(range(4))
    source['mapping']=selected;source['mapping_sha256']=plan['mapping_table_sha256'];verify(bindings)
    return source,bindings


def replay(source,cid):
    entry=source['details'][cid];native=source['native'][cid];chain=native['chain']
    base=dict(chain_id=cid,seed=entry['seed'],model_input_identity=entry['model_input_identity'],
        effective_input_group=chain['effective_input_group'],scientific_eligibility=False)
    if entry['status']=='unresolved_failed_native_chain_retained':
        return dict(**base,status=entry['status'],native_disposition=entry['native_disposition'],node_rows=[],candidate_rows=[])
    source['doc'](entry['sample_audit'],entry['sample_audit_sha256'])
    assert sha(chain['tree'])==chain['tree_sha256']
    matched=match_trees(Path(chain['tree']).read_text(),Path(entry['native_files']['runtime-tree.nwk']['path']).read_text())
    with Path(entry['input_alignment']).open() as handle:observed=fasta_records(handle)
    assert sorted(observed)==matched['tips'] and len(observed)==entry['tips']==chain['proteins']
    assert len(matched['rows'])==entry['runtime_nodes']
    dataset='whole' if '-whole-' in chain['original_configuration_ids'][0] else 'domain'
    mapping=sorted(source['mapping'][('profile',chain['family'],dataset)],key=lambda r:int(r['level']))
    audit=source['doc'](entry['sample_audit'],entry['sample_audit_sha256'])
    assert audit['mapping_sha256']==source['mapping_sha256']
    candidate_rows=[];derived={};tip_bits={tip:1<<i for i,tip in enumerate(matched['tips'])}
    for row in mapping:
        node=row['source_node'];retained=json.loads(row['retained_set_json'])
        assert len(retained)==len(set(retained))==int(row['retained_descendants']) and retained
        assert set(retained)<=set(tip_bits)
        mask=sum(tip_bits[tip] for tip in retained)
        assert matched['source_labels'][node]['mask']==mask
        target=matched['runtime_index'][mask]['label'];derived[node]=target
        level=int(row['level']);is_root=mask==(1<<len(tip_bits))-1
        assert is_root==(level==3)
        candidate_rows.append(dict(level=level,source_node=node,runtime_node=target,descendant_mask_hex=hex(mask),
            retained_descendants=len(retained),retained_tips_sha256=digest(sorted(retained)),
            assumed_root=is_root,scientific_eligibility=False))
    assert len(derived)==len(set(derived.values()))==4
    assert derived==entry['source_to_runtime_candidates']
    seen=set()
    for row in audit['candidate_samples']:
        key=row['iteration'],row['source_node'];assert key not in seen;seen.add(key)
        assert row['iteration'] in range(0,1001,10)
        assert row['runtime_node']==derived[row['source_node']]
        assert int(row['level'])==next(r['level'] for r in candidate_rows if r['source_node']==row['source_node'])
    assert seen=={(i,n) for i in range(0,1001,10) for n in derived}
    return dict(**base,status='all_source_runtime_rooted_clades_and_candidate_mappings_checked_conditional_on_root',
        source_tree=dict(path=chain['tree'],sha256=chain['tree_sha256']),runtime_tree=entry['native_files']['runtime-tree.nwk'],
        input_alignment_sha256=entry['input_alignment_sha256'],sample_audit_sha256=entry['sample_audit_sha256'],
        mapping_table_sha256=source['mapping_sha256'],ordered_tips=matched['tips'],node_rows=matched['rows'],candidate_rows=candidate_rows,
        candidate_frame_mappings=len(seen),maximum_absolute_length_difference=matched['maximum_absolute_length_difference'],
        negative_branch_pairs_require_review=matched['negative_branch_pairs_require_review'],
        source_root=matched['source_root'],runtime_root=matched['runtime_root'],biological_root_accepted=False,
        native_saved_alignments_decoded=False)


def summarize(results,source,plan):
    assert set(results)==set(source['details']);counts=Counter();groups=set();effective=set()
    for cid,row in results.items():
        entry=source['details'][cid];assert row['scientific_eligibility'] is False
        assert row['chain_id']==cid and row['seed']==entry['seed'] and row['model_input_identity']==entry['model_input_identity']
        groups.add(row['model_input_identity']);effective.add(row['effective_input_group'])
        if entry['status']=='unresolved_failed_native_chain_retained':
            assert row['status']==entry['status'] and row['node_rows']==row['candidate_rows']==[]
            counts['failed_chains']+=1;continue
        assert row['status']=='all_source_runtime_rooted_clades_and_candidate_mappings_checked_conditional_on_root'
        counts['checked_chains']+=1;counts['node_pairs']+=len(row['node_rows'])
        counts['nonroot_branch_pairs']+=sum(not r['is_root'] for r in row['node_rows'])
        counts['candidate_node_pairs']+=len(row['candidate_rows']);counts['candidate_frame_mappings']+=row['candidate_frame_mappings']
        counts['assumed_root_candidates']+=sum(r['assumed_root'] for r in row['candidate_rows'])
        counts['negative_branch_pairs_require_review']+=row['negative_branch_pairs_require_review']
        counts['intact_chains_in_unresolved_groups']+=source['groups'][row['model_input_identity']]['status']=='unresolved_failed_native_chain_retained'
    result=dict(full_chains=len(results),full_groups=len(groups),complete_groups=403,unresolved_groups=2,effective_input_groups=len(effective),**dict(counts))
    assert result['checked_chains']==1618 and result['failed_chains']==2 and result['full_groups']==405
    assert result['candidate_node_pairs']==6472 and result['candidate_frame_mappings']==653672
    assert result['assumed_root_candidates']==1618 and result['intact_chains_in_unresolved_groups']==6
    assert result['nonroot_branch_pairs']==result['node_pairs']-1618
    assert all(result[k]==v for k,v in plan['expected'].items())
    return result
