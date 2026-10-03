"""Complete source fixture for independent rooted native clade verification."""
from collections import Counter
import csv
import json
import os
from pathlib import Path

import numpy as np

from check_independent_baliphy_scalar_readback_v2 import setup as scalar_fixture
from independent_native_ancestral_alignment import ALPHABET
from prepare_ancestral_state_traces import project_states as original_project
from prepare_independent_baliphy_scalar_readback_v2 import atomic
from run_ortholog_pair_guide_comparison import sha

def setup(root):
    path=scalar_fixture(root);plan=json.loads(path.read_text());source=root/'source'
    complete_path=Path(plan['diagnostic_completion']);complete=json.loads(complete_path.read_text())
    archive_path=Path(complete['full_hash_archive']);archive=json.loads(archive_path.read_text());bindings=archive['source_hashes']
    def create(p,value):atomic(p,value);bindings[str(p)]=sha(p);return p
    templates=source/'native-templates';templates.mkdir();iterations=list(range(0,1001,10));variants=[]
    observed={f't{i}':'ACDX' for i in range(5)};mapping={f'source-{i}':f'node-{i}' for i in range(4)}
    input_path=templates/'input.faa';input_path.write_text(''.join('>'+t+' description\nACD-X\n' for t in observed));bindings[str(input_path)]=sha(input_path)
    tree_path=templates/'runtime-tree.nwk';tree_path.write_text('((((t0:1,t1:1)node-0:1,t2:1)node-1:1,t3:1)node-2:1,t4:1)node-3;\n');bindings[str(tree_path)]=sha(tree_path)
    source_tree=templates/'source-tree.nwk';source_tree.write_text(tree_path.read_text().replace('node-','source-'));bindings[str(source_tree)]=sha(source_tree)
    mapping_path=templates/'mapping.tsv'
    with mapping_path.open('w') as h:
        writer=csv.DictWriter(h,fieldnames=['guide','family','dataset','level','source_node','retained_descendants','retained_set_json'],delimiter='\t',lineterminator='\n');writer.writeheader()
        for number in range(135):
            for level in range(4):writer.writerow(dict(guide='profile',family=f'family-{number:03}',dataset='whole' if number%2==0 else 'domain',level=level,source_node=f'source-{level}',retained_descendants=level+2,retained_set_json=json.dumps([f't{i}' for i in range(level+2)])))
    bindings[str(mapping_path)]=sha(mapping_path)
    coords=dict(nodes=sorted(mapping),tips=[dict(tip=t,length=4) for t in sorted(observed)],alphabet=ALPHABET,input_alignment=str(input_path),
        input_alignment_sha256=sha(input_path),
        ordering='state axes: saved iteration, source node, concatenated ungapped input residues in listed tip order; positions are one-based within each tip')
    for variant in range(4):
        rng=np.random.default_rng(8100+variant);rows=[];states=[];free=[];candidates=[]
        for index,iteration in enumerate(iterations):
            positions=np.sort(rng.choice(9,size=4,replace=False));tip=['-']*9
            for position,residue in zip(positions,'ACD'+('E' if index%2 else 'X')):tip[position]=residue
            sequences={t:''.join(tip) for t in observed}
            for i in range(4):sequences['node-'+str(i)]=''.join(ALPHABET[n] for n in rng.integers(0,22,size=9))
            # Expected arrays come from the original projection, not the new decoder.
            a,b=original_project(sequences,observed,mapping);states.append(a);free.append(b)
            candidates.extend(dict(iteration=iteration,level=str(int(n.split('-')[-1])),source_node=n,runtime_node=mapping[n],
                ungapped_length=len(sequences[mapping[n]].replace('-',''))) for n in sorted(mapping))
            rows.append('iterations = '+str(iteration)+'\n\n'+''.join('>'+n+' description\n'+s[:4].lower()+'\n'+s[4:].lower()+'\n'
                for n,s in (reversed(list(sequences.items())) if index%2 else sequences.items())))
        native=templates/('native-'+str(variant)+'.fastas');native.write_text('\n'.join(rows));bindings[str(native)]=sha(native)
        values=np.asarray(states,dtype=np.uint8);extra=np.asarray(free,dtype=np.int32)
        arrays=templates/('states-'+str(variant)+'.npz')
        counts={f'counts_after_{cut}':np.stack([np.count_nonzero(values[np.asarray(iterations)>cut]==s,axis=0) for s in range(22)],axis=-1).astype(np.uint16) for cut in [250,500]}
        np.savez_compressed(arrays,states=values,iterations=np.asarray(iterations),unanchored_residue_counts=extra,**counts)
        bindings[str(arrays)]=sha(arrays)
        variants.append(dict(native=native,arrays=arrays,candidates=candidates,unanchored=int(extra.sum())))
    recovery_path=Path(plan['recovery_completion']);recovery=json.loads(recovery_path.read_text())
    overlay_path=Path(recovery['overlay']);overlay=json.loads(overlay_path.read_text())
    producer_path=Path(complete['producer_receipt']);producer=json.loads(producer_path.read_text());states={};details={};totals=Counter()
    for row in overlay['rows']:
        chain=row['chain'];cid=chain['chain_id'];group=row['model_input_identity'];selected=row['selected_disposition']
        number=int(group.split('-')[-1])//3
        chain.update(effective_input_group=f'effective-{number:03}',family=f'family-{number:03}',tree=str(source_tree),tree_sha256=sha(source_tree),original_configuration_ids=[f'family-{number:03}-'+('whole' if number%2==0 else 'domain')+'-synthetic'])
        variant=variants[(chain['seed']-1)%4]
        entry=dict(chain_id=cid,model_input_identity=group,seed=chain['seed'],native_disposition=selected['status'],scientific_eligibility=False)
        if selected['status']!='all_saved_alignments_and_candidate_nodes_checked':
            entry['status']='unresolved_failed_native_chain_retained';details[cid]=entry;totals['failed_chains']+=1
            states[cid]={};continue
        chain.update(alignment=str(input_path),alignment_sha256=sha(input_path),proteins=5)
        folder=source/group/(cid+'-native');folder.mkdir()
        native=folder/'C1.P1.fastas';os.link(variant['native'],native);bindings[str(native)]=sha(native)
        tree=folder/'runtime-tree.nwk';os.link(tree_path,tree);bindings[str(tree)]=sha(tree)
        receipt=create(folder/'receipt.json',dict(exit_code=0,artifacts={'C1.P1.fastas':sha(native),'runtime-tree.nwk':sha(tree)}))
        audit_path=Path(selected['sample_audit']);audit=json.loads(audit_path.read_text())
        audit.update(iterations=1000,mapping_sha256=sha(mapping_path),attempt_receipt=str(receipt),attempt_receipt_sha256=sha(receipt),candidate_samples=variant['candidates'])
        create(audit_path,audit);selected['sample_audit_sha256']=sha(audit_path)
        ap=folder/'states.npz';os.link(variant['arrays'],ap);bindings[str(ap)]=sha(ap)
        cp=create(folder/'coordinates.json',coords)
        sr=create(folder/'states-receipt.json',dict(summaries=[dict(samples=101,source_audit_sha256=sha(audit_path),
            candidate_anchor_coordinates=80,state_observations=8080,unanchored_residue_observations=variant['unanchored'],
            artifacts={str(ap):sha(ap),str(cp):sha(cp)})]))
        states[cid]=dict(receipt=str(sr),receipt_sha256=sha(sr))
        entry.update(status='full_native_source_metadata_and_runtime_labels_checked_sample_decode_pending',
            sample_audit=str(audit_path),sample_audit_sha256=sha(audit_path),selected_attempt_receipt=str(receipt),selected_attempt_receipt_sha256=sha(receipt),
            input_alignment=str(input_path),input_alignment_sha256=sha(input_path),
            native_files={name:dict(path=str(p),sha256=sha(p),bytes_at_observation=p.stat().st_size) for name,p in [('C1.P1.fastas',native),('runtime-tree.nwk',tree)]},
            source_to_runtime_candidates=mapping,original_state_array=dict(path=str(ap),sha256=sha(ap)),coordinates=dict(path=str(cp),sha256=sha(cp)),
            tips=5,runtime_nodes=9,expected_saved_alignments=101,candidate_anchor_coordinates=80,state_observations=8080,
            inherited_unanchored_residue_observations=variant['unanchored'])
        details[cid]=entry
        totals.update(dict(checked_chains=1,saved_alignments_to_decode=101,candidate_frames_to_decode=404,
            state_observations=8080,candidate_anchor_coordinates=80,inherited_unanchored_residue_observations=variant['unanchored'],
            native_alignment_bytes=native.stat().st_size))
    create(overlay_path,overlay);recovery['overlay_sha256']=sha(overlay_path);create(recovery_path,recovery)
    producer['states']=states;create(producer_path,producer)
    archive['source_hashes']=bindings;atomic(archive_path,archive)
    complete.update(bound_source_hashes=len(bindings),full_hash_archive_sha256=sha(archive_path),producer_receipt_sha256=sha(producer_path));atomic(complete_path,complete)
    dp=source/'native-inventory-details.json';atomic(dp,details)
    hp=source/'native-inventory-hashes.json';atomic(hp,dict(status='full_native_alignment_inventory_source_hash_archive',source_hashes=bindings))
    ip=source/'native-inventory.json'
    atomic(ip,dict(status='complete_full_independent_native_alignment_inventory_sample_replay_pending',scientific_eligibility=False,
        full_chains=1620,checked_chains=1618,failed_chains=2,full_groups=405,complete_groups=403,totals=dict(totals),
        full_chain_details=str(dp),full_chain_details_sha256=sha(dp),full_source_hash_archive=str(hp),full_source_hash_archive_sha256=sha(hp),source_hashes_rechecked=len(bindings)))
    ic=source/'native-inventory-completed.json'
    atomic(ic,dict(status='verified_original_native_alignment_inventory_completion',inventory=str(ip),inventory_sha256=sha(ip),
        details_sha256=sha(dp),source_hash_archive_sha256=sha(hp),bound_sources=len(bindings),native_alignment_samples_replayed=False,scientific_eligibility=False,
        original_terminal_handle=dict(status='verified_original_terminal_success_with_bound_completed_artifacts',scope='Synthetic fixture, not actual process evidence.')))
    plan.update(native_inventory=str(ip),inventory_completion=str(ic),mapping_table=str(mapping_path),mapping_table_sha256=sha(mapping_path),expected=dict(full_chains=1620,checked_chains=1618,failed_chains=2,
        saved_alignments=163418,candidate_frames=653672,state_observations=totals['state_observations'],candidate_anchor_coordinates=totals['candidate_anchor_coordinates'],
        unanchored_residue_observations=totals['inherited_unanchored_residue_observations'],cutoff_count_cells=44*totals['candidate_anchor_coordinates'],intact_chains_in_unresolved_groups=6),
        scope='Complete synthetic1620chain/405group native workflow, four genuinely varying projection templates. Original projection creates expected arrays. Synthetic source journals are not production provenance.')
    plan['expected']=dict(full_chains=1620,checked_chains=1618,failed_chains=2,full_groups=405,complete_groups=403,unresolved_groups=2,effective_input_groups=135,node_pairs=14562,nonroot_branch_pairs=12944,candidate_node_pairs=6472,candidate_frame_mappings=653672,assumed_root_candidates=1618,negative_branch_pairs_require_review=0,intact_chains_in_unresolved_groups=6)
    plan['scope']='Full1620chain/405group/135input synthetic rooted-clade software workflow, both failures/six intact unresolved retained; no actual native journals or biological root qualification.'
    atomic(path,plan);return path
