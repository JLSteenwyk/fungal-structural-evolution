#!/usr/bin/env python3
"""Full 1620-chain native decoding/checkpoint/source-corruption contracts."""
from collections import Counter
import argparse
import copy
import json
import os
from pathlib import Path
import tempfile

import numpy as np

from check_independent_baliphy_scalar_readback_v2 import setup as scalar_fixture
from independent_native_ancestral_alignment import ALPHABET
from independent_native_alignment_replay_sources import load
from independent_native_alignment_replay import replay_chain, count_states
from prepare_ancestral_state_traces import project_states as original_project
from prepare_independent_native_alignment_replay import run
from readback_independent_native_alignment_replay import run as readback
from prepare_independent_baliphy_scalar_readback_v2 import atomic
from run_ortholog_pair_guide_comparison import sha


def setup(root):
    path=scalar_fixture(root);plan=json.loads(path.read_text());source=root/'source'
    complete_path=Path(plan['diagnostic_completion']);complete=json.loads(complete_path.read_text())
    archive_path=Path(complete['full_hash_archive']);archive=json.loads(archive_path.read_text());bindings=archive['source_hashes']
    def create(p,value):atomic(p,value);bindings[str(p)]=sha(p);return p
    templates=source/'native-templates';templates.mkdir();iterations=list(range(0,1001,10));variants=[]
    observed={'tip':'ACDX'};mapping={f'source-{i}':f'node-{i}' for i in range(4)}
    input_path=templates/'input.faa';input_path.write_text('>tip description\nACD-X\n');bindings[str(input_path)]=sha(input_path)
    tree_path=templates/'runtime-tree.nwk';tree_path.write_text('((((tip:1)node-0:1)node-1:1)node-2:1)node-3;\n');bindings[str(tree_path)]=sha(tree_path)
    coords=dict(nodes=sorted(mapping),tips=[dict(tip='tip',length=4)],alphabet=ALPHABET,input_alignment=str(input_path),
        input_alignment_sha256=sha(input_path),
        ordering='state axes: saved iteration, source node, concatenated ungapped input residues in listed tip order; positions are one-based within each tip')
    for variant in range(4):
        rng=np.random.default_rng(8100+variant);rows=[];states=[];free=[];candidates=[]
        for index,iteration in enumerate(iterations):
            positions=np.sort(rng.choice(9,size=4,replace=False));tip=['-']*9
            for position,residue in zip(positions,'ACD'+('E' if index%2 else 'X')):tip[position]=residue
            sequences={'tip':''.join(tip)}
            for i in range(4):sequences['node-'+str(i)]=''.join(ALPHABET[n] for n in rng.integers(0,22,size=9))
            # Expected arrays come from the original projection, not the new decoder.
            a,b=original_project(sequences,observed,mapping);states.append(a);free.append(b)
            candidates.extend(dict(iteration=iteration,source_node=n,runtime_node=mapping[n],
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
        variant=variants[(chain['seed']-1)%4]
        entry=dict(chain_id=cid,model_input_identity=group,seed=chain['seed'],native_disposition=selected['status'],scientific_eligibility=False)
        if selected['status']!='all_saved_alignments_and_candidate_nodes_checked':
            entry['status']='unresolved_failed_native_chain_retained';details[cid]=entry;totals['failed_chains']+=1
            states[cid]={};continue
        chain.update(alignment=str(input_path),alignment_sha256=sha(input_path),proteins=1)
        folder=source/group/(cid+'-native');folder.mkdir()
        native=folder/'C1.P1.fastas';os.link(variant['native'],native);bindings[str(native)]=sha(native)
        tree=folder/'runtime-tree.nwk';os.link(tree_path,tree);bindings[str(tree)]=sha(tree)
        receipt=create(folder/'receipt.json',dict(exit_code=0,artifacts={'C1.P1.fastas':sha(native),'runtime-tree.nwk':sha(tree)}))
        audit_path=Path(selected['sample_audit']);audit=json.loads(audit_path.read_text())
        audit.update(iterations=1000,attempt_receipt=str(receipt),attempt_receipt_sha256=sha(receipt),candidate_samples=variant['candidates'])
        create(audit_path,audit);selected['sample_audit_sha256']=sha(audit_path)
        ap=folder/'states.npz';os.link(variant['arrays'],ap);bindings[str(ap)]=sha(ap)
        cp=create(folder/'coordinates.json',coords)
        sr=create(folder/'states-receipt.json',dict(summaries=[dict(samples=101,source_audit_sha256=sha(audit_path),
            candidate_anchor_coordinates=16,state_observations=1616,unanchored_residue_observations=variant['unanchored'],
            artifacts={str(ap):sha(ap),str(cp):sha(cp)})]))
        states[cid]=dict(receipt=str(sr),receipt_sha256=sha(sr))
        entry.update(status='full_native_source_metadata_and_runtime_labels_checked_sample_decode_pending',
            sample_audit=str(audit_path),sample_audit_sha256=sha(audit_path),selected_attempt_receipt=str(receipt),selected_attempt_receipt_sha256=sha(receipt),
            input_alignment=str(input_path),input_alignment_sha256=sha(input_path),
            native_files={name:dict(path=str(p),sha256=sha(p),bytes_at_observation=p.stat().st_size) for name,p in [('C1.P1.fastas',native),('runtime-tree.nwk',tree)]},
            source_to_runtime_candidates=mapping,original_state_array=dict(path=str(ap),sha256=sha(ap)),coordinates=dict(path=str(cp),sha256=sha(cp)),
            tips=1,runtime_nodes=5,expected_saved_alignments=101,candidate_anchor_coordinates=16,state_observations=1616,
            inherited_unanchored_residue_observations=variant['unanchored'])
        details[cid]=entry
        totals.update(dict(checked_chains=1,saved_alignments_to_decode=101,candidate_frames_to_decode=404,
            state_observations=1616,candidate_anchor_coordinates=16,inherited_unanchored_residue_observations=variant['unanchored'],
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
    plan.update(native_inventory=str(ip),inventory_completion=str(ic),expected=dict(full_chains=1620,checked_chains=1618,failed_chains=2,
        saved_alignments=163418,candidate_frames=653672,state_observations=totals['state_observations'],candidate_anchor_coordinates=totals['candidate_anchor_coordinates'],
        unanchored_residue_observations=totals['inherited_unanchored_residue_observations'],cutoff_count_cells=44*totals['candidate_anchor_coordinates'],intact_chains_in_unresolved_groups=6),
        scope='Complete synthetic1620chain/405group native workflow, four genuinely varying projection templates. Original projection creates expected arrays. Synthetic source journals are not production provenance.')
    atomic(path,plan);return path


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    rng=np.random.default_rng(250500);histograms=0
    for draws in [50,75]:
        for length in [1,7,4097]:
            values=rng.integers(0,22,(draws,4,length),dtype=np.uint8)
            expected=np.stack([np.count_nonzero(values==s,axis=0) for s in range(22)],axis=-1)
            assert np.array_equal(count_states(values),expected);histograms+=1
    rejected=[];source_rejected=[]
    with tempfile.TemporaryDirectory(prefix='independent-native-replay-contract-') as directory:
        root=Path(directory);path=setup(root)
        try:run(path,stop_after_chains=3)
        except InterruptedError:pass
        else:raise AssertionError('Interruption not exercised')
        before={str(p):sha(p) for p in (root/'output/chains').glob('*')}
        result=run(path);assert all(sha(p)==v for p,v in before.items())
        assert result['saved_alignments']==163418 and result['intact_chains_in_unresolved_groups']==6
        readback(path,root/'output/readback.json')
        for action in ['completed_producer','alternate_reader']:
            try:run(path) if action=='completed_producer' else readback(path,root/'alternate.json')
            except AssertionError:pass
            else:raise AssertionError('Completed restart accepted')
        (root/'output/readback.json').unlink()
        rootout=root/'output';manifest=rootout/'chain_manifest.json';receipt=rootout/'receipt.json'
        mb,rb=manifest.read_bytes(),receipt.read_bytes();entries=json.loads(mb);first=entries[0]
        cp=rootout/first['path'];cb=cp.read_bytes();record=json.loads(cb);frames=rootout/'chains'/record['result']['frames_path'];fb=frames.read_bytes()
        rows=[json.loads(line) for line in fb.decode().splitlines()]
        changes=['missing_frame','duplicate_frame','alignment_digest','projection_digest','unanchored_count','node_length',
            'promote_science','cutoff_digest','full_state_digest','foreign_seed','false_failure','hide_cutoff','omit_failed_chain','omit_intact_unresolved_chain']
        for change in changes:
            altered=copy.deepcopy(rows);document=json.loads(cb);parts=json.loads(mb)
            if change=='missing_frame':altered.pop()
            elif change=='duplicate_frame':altered[-1]=copy.deepcopy(altered[0])
            elif change=='alignment_digest':altered[0]['alignment_sha256']='0'*64
            elif change=='projection_digest':altered[0]['projected_states_sha256']='0'*64
            elif change=='unanchored_count':altered[0]['unanchored_residue_counts'][0]+=1
            elif change=='node_length':altered[0]['node_lengths']['source-0']+=1
            elif change=='promote_science':altered[0]['scientific_eligibility']=True
            elif change=='cutoff_digest':document['result']['cutoffs']['250']['counts_sha256']='0'*64
            elif change=='full_state_digest':document['result']['full_projected_states_sha256']='0'*64
            elif change=='foreign_seed':document['result']['seed']+=1
            elif change=='false_failure':document['result']['status']='unresolved_failed_native_chain_retained';document['result']['cutoffs']={}
            elif change=='hide_cutoff':del document['result']['cutoffs']['500']
            elif change=='omit_failed_chain':parts=[r for r in parts if r['chain_id']!='group-0404-chain0']
            elif change=='omit_intact_unresolved_chain':parts=[r for r in parts if r['chain_id']!='group-0404-chain1']
            frames.write_text(''.join(json.dumps(r)+'\n' for r in altered));document['result']['frames_sha256']=sha(frames);atomic(cp,document)
            parts[0]['sha256']=sha(cp);atomic(manifest,parts);report=json.loads(rb)
            for artifact in [cp,frames,manifest]:report['artifacts'][str(artifact.relative_to(rootout))]=sha(artifact)
            atomic(receipt,report)
            try:readback(path,root/'bad-reader.json')
            except (AssertionError,ValueError):rejected.append(change)
            else:raise AssertionError('Rehashed false export accepted: '+change)
            frames.write_bytes(fb);cp.write_bytes(cb);manifest.write_bytes(mb);receipt.write_bytes(rb)
        source,_=load(json.loads(path.read_text()),path);cid=first['chain_id'];entry=source['details'][cid]
        original_path=Path(entry['original_state_array']['path']);array_bytes=original_path.read_bytes()
        with np.load(original_path,allow_pickle=False) as saved:arrays={k:saved[k] for k in saved.files}
        for change in ['projected_state','unanchored_array','cutoff_250_counts','cutoff_500_counts','saved_iterations']:
            altered={k:v.copy() for k,v in arrays.items()}
            if change=='projected_state':altered['states'][0,0,0]=(int(altered['states'][0,0,0])+1)%22
            elif change=='unanchored_array':altered['unanchored_residue_counts'][0,0]+=1
            elif change=='cutoff_250_counts':altered['counts_after_250'][0,0,0]+=1
            elif change=='cutoff_500_counts':altered['counts_after_500'][0,0,21]+=1
            else:altered['iterations'][0]=10
            original_path.unlink();np.savez_compressed(original_path,**altered)
            try:replay_chain(source,cid,frames,readback=True)
            except (AssertionError,ValueError):source_rejected.append(change)
            else:raise AssertionError('Changed source array accepted: '+change)
            original_path.write_bytes(array_bytes)
        native=Path(entry['native_files']['C1.P1.fastas']['path']);native_bytes=native.read_bytes()
        for change in ['duplicate_native_iteration','changed_observed_residue','missing_runtime_node']:
            text=native_bytes.decode()
            if change=='duplicate_native_iteration':text=text.replace('iterations = 10\n','iterations = 0\n',1)
            elif change=='changed_observed_residue':
                # Change the first concrete A residue in the tip record only.
                label=text.index('>tip description\n');end=text.index('>',label+1)
                part=text[label:end];text=text[:label]+part.replace('a','c',1)+text[end:]
            else:
                label=text.index('>node-0 description\n');end=text.index('>',label+1)
                text=text[:label]+text[end:]
            native.unlink();native.write_text(text)
            try:replay_chain(source,cid,frames,readback=True)
            except (AssertionError,ValueError):source_rejected.append(change)
            else:raise AssertionError('Changed native sample accepted: '+change)
            native.write_bytes(native_bytes)
    paths=[__file__,*[f'scripts/{n}.py' for n in ['independent_native_ancestral_alignment','independent_native_alignment_replay_sources',
        'independent_native_alignment_replay','prepare_independent_native_alignment_replay','readback_independent_native_alignment_replay',
        'check_independent_baliphy_scalar_readback_v2','prepare_ancestral_state_traces','ancestral_residue_anchors']]]
    result=dict(status='passed_full_independent_native_alignment_replay_contracts',histogram_fixtures=histograms,
        full_chains=1620,checked_chains=1618,failed_chains=2,full_groups=405,intact_chains_in_unresolved_groups=6,
        saved_alignments=result['saved_alignments'],candidate_frames=result['candidate_frames'],state_observations=result['state_observations'],
        cutoff_count_cells=result['cutoff_count_cells'],unchanged_interrupted_checkpoints_reused=True,completed_producer_and_alternate_reader_restart_refused=True,
        rehashed_false_exports_rejected=rejected,source_array_and_native_alterations_rejected=source_rejected,
        source_hashes={str(p):sha(p) for p in paths},scientific_eligibility=False,production_native_samples_replayed=False,
        scope='Complete1620chain/405group software source/producer/serialized-reader/checkpoint grid with four varying native templates. '
              'Expected states/counts use original projection and separate per-label histograms, not new decoder. Source journals are synthetic; no production/posterior proof.')
    with args.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
