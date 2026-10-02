#!/usr/bin/env python3
"""Observe all prior live identities and new handoffs without accepting partial output."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
from background_measurement_union_sources import closed_source
from reference_measurement_union_sources import verify
from record_project_runtime_checkpoint_v4 import fingerprint,live_record,journal_terminal
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    basis=Path('metadata/project_live_launch_inventory_20261002_v3.json');inventory=json.loads(basis.read_text())
    stages={};bindings={};completion_paths=[('metadata/full_triad_context_joint_directions_completed_20261002.json','complete_verified_full_triad_context_joint_directions'),('metadata/full_matching_case_index_completed_20261002.json','complete_verified_full_matching_logical_case_index')]
    for path,status in completion_paths:
        if not Path(path).exists():continue
        c=closed_source(path,status,status+'_archive',2,bindings)
        assert c['scientific_eligibility'] is False
        stages[path]=dict(status=status,completion_sha256=sha(path),full_hash_archive=c['full_hash_archive'],full_hash_archive_sha256=c['full_hash_archive_sha256'],bound_source_hashes=c['bound_source_hashes'],exact_original_journals=2)
    verify(bindings)
    handles={}
    for r in inventory['original_live_handles']:
        for ref in r['launch_records']:assert sha(ref['path'])==ref['sha256'],ref['path']
        launch=r['launch_records'][0]['path'];handles[r['pid']]=dict(pid=r['pid'],created=r['created'],cmdline=r['cmdline'],launch=launch,origin='prior_original_live_inventory')
    new_inventories=['metadata/full_matching_case_index_launches_20261002.json','metadata/full_expanded_measurement_catalog_launches_20261002.json','metadata/full_expanded_measurement_catalog_v2_launches_20261002.json','metadata/full_expanded_case_measurements_launches_20261002.json','metadata/full_expanded_measurement_catalog_v3_launches_20261002.json','metadata/full_expanded_case_measurements_v2_launches_20261002.json']
    for path in new_inventories:
        for launch in json.loads(Path(path).read_text())['launches']:
            r=json.loads(Path(launch).read_text());assert r['pid'] not in handles
            handles[r['pid']]=dict(pid=r['pid'],created=r['created'],cmdline=r['cmdline'],launch=launch,origin=path)
    failed={}
    for path in ['metadata/full_expanded_measurement_catalog_attempt_20261002.json','metadata/full_expanded_measurement_catalog_attempt_20261002_v2.json']:
        for r in json.loads(Path(path).read_text())['original_failed_handles']:failed[r['launch']]=r
    live=[];terminal=[];failures=[]
    for r in handles.values():
        launch=json.loads(Path(r['launch']).read_text());assert all(r[k]==launch[k] for k in ['pid','created','cmdline'])
        if 'plan' in launch:assert sha(launch['plan'])==launch['plan_sha256']
        proc=fingerprint(r)
        if proc is not None:
            assert r['launch'] not in failed;live.append(live_record(r,proc))
        elif r['launch'] in failed:failures.append(journal_terminal(r,expected_failure=True))
        else:
            record=journal_terminal(r)
            if r['launch'].startswith('metadata/full_matching_case_index'):
                complete='metadata/full_matching_case_index_completed_20261002.json' in stages
                record['status']='original_success_with_full_case_provenance_verified' if complete else 'original_case_stage_success_pending_full_provenance_closure'
            elif 'full_triad_context_joint_directions' in r['launch']:
                assert 'metadata/full_triad_context_joint_directions_completed_20261002.json' in stages
                record['status']='original_success_with_full_context_joint_provenance_verified'
            else:raise AssertionError(('New terminal handle requires a completed-stage audit',r['launch']))
            terminal.append(record)
    result=dict(status='verified_expanded_measurement_execution_and_original_live_inventory',checked_utc=datetime.now(timezone.utc).isoformat(),basis_inventory=str(basis),basis_inventory_sha256=sha(basis),source_completions_fully_reverified=stages,distinct_closed_source_artifact_bindings_checked=len(bindings),handles_checked=len(handles),exact_live_handles=len(live),new_terminal_successes=len(terminal),preserved_new_failed_invocations=len(failures),live_handles=live,terminal_success_handles=terminal,preserved_failed_handles=failures,source_pins={p:sha(p) for p in [__file__,str(basis),*new_inventories,'scripts/record_project_runtime_checkpoint_v4.py']},gpu_inference_paused=True,all_eight_aims_incomplete=True,scope='Recheck every one of31prior original live identities plus18new index/catalog/join handles. Completed joint-context and available completed case proofs require complete archive hashes/source bindings, not only locator hashes. A passed case producer/reader with closure still live remains provisional. All nine new failed catalog/dependency invocations are preserved with actual original journals; they are not native prediction failures. Older70terminal successes/ten historical failures remain at the separately linked prior full runtime checkpoint, not counted as newly reverified here. No job mutation, restarts, GPU activation or new charges.')
    with a.output.open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['live_handles','terminal_success_handles','preserved_failed_handles','scope','source_pins']},indent=2),flush=True)


if __name__=='__main__':main()
