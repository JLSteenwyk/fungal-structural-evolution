#!/usr/bin/env python3
"""Read-only audit of original operator/input/design/recovery handles and receipts."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from record_project_runtime_checkpoint_v4 import fingerprint, live_record, journal_terminal
from run_ortholog_pair_guide_comparison import sha


def closed(path, expected):
    path = Path(path)
    if not path.exists():
        return dict(accepted=False, expected_completion=str(path), reason='Final compact completion is absent')
    receipt = json.loads(path.read_text())
    assert receipt['status'] == expected and receipt['scientific_eligibility'] is False
    archive = Path(receipt['full_hash_archive'])
    assert sha(archive) == receipt['full_hash_archive_sha256']
    proof = json.loads(archive.read_text())
    assert proof['status'] == expected+'_archive' and proof['scientific_eligibility'] is False
    assert len(proof['source_hashes']) == receipt['bound_source_hashes']
    assert len(proof['services']) == receipt['exact_process_journals_checked'] == 2
    assert all(receipt[k] == v for k,v in proof['summary'].items())
    for key in ['producer_receipt', 'independent_readback']:
        assert sha(receipt[key]) == receipt[key+'_sha256'] == proof['source_hashes'][receipt[key]]
    journals=[]
    for service in proof['services']:
        launch_path=service['launch']; launch=json.loads(Path(launch_path).read_text())
        assert sha(launch_path)==proof['source_hashes'][launch_path]
        launch['launch']=launch_path
        observed=journal_terminal(launch)
        # journalctl JSON object key order can vary between reads. The sealed
        # historical byte digest remains evidence of its original observation;
        # compare the original semantic completion evidence on this new read.
        observed['archived_journal_sha256_at_original_observation']=service['journal_sha256']
        assert observed['captured_process_messages']==service['captured_process_messages']
        assert observed['completion_resource_records']==service['completion_resource_records']
        journals.append(observed)
    return dict(accepted=True, path=str(path), sha256=sha(path), status=receipt['status'],
        full_hash_archive=str(archive), full_hash_archive_sha256=receipt['full_hash_archive_sha256'],
        original_closure_full_bindings_checked=receipt['bound_source_hashes'],
        original_journals_revalidated=journals, summary=proof['summary'],
        scope='Sealed full archive checksum, compact/summary/receipt linkage and actual original '
              'completion/resource journals revalidated. Full source and artifact hashing was '
              'performed by the original producer, independent reader and closure; it is not '
              'repeated here. No scientific model or effect is accepted by this handoff.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--previous',type=Path,default=Path('metadata/full_entity_operator_queue_checkpoint_20261002.json'))
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(); previous=json.loads(args.previous.read_text())
    checked_pins={}
    assert sha(previous['source_plan'])==previous['source_plan_sha256']
    plan=json.loads(Path(previous['source_plan']).read_text())
    for path, digest in plan['pins'].items():
        assert sha(path)==digest,path
        checked_pins[path]=digest
    handles=[]
    for record in previous['handles']:
        launch=json.loads(Path(record['launch']).read_text())
        assert all(record[k]==launch[k] for k in ['pid','created','cmdline'])
        assert sha(launch['plan'])==launch['plan_sha256']
        proc=fingerprint(record)
        if proc is not None:handles.append(live_record(record,proc))
        else:
            terminal=journal_terminal(record)
            terminal['status']='verified_original_terminal_success_resource_journal'
            handles.append(terminal)
    software=[]
    for path in ['metadata/covariance_basis_audit_validation_20261002_v2.json',
                 'metadata/full_entity_operator_target_node_redundancy_20261002_v2.json']:
        report=json.loads(Path(path).read_text())
        for source,digest in report['source_hashes'].items():assert sha(source)==digest,source
        software.append(dict(path=path,sha256=sha(path),status=report['status']))
    stages={key:closed(path,status) for key,path,status in [
        ('model_inputs','metadata/full_expanded_model_inputs_completed_20261002.json','complete_verified_full_expanded_model_inputs'),
        ('operator_bank','metadata/full_entity_operator_bank_completed_20261002.json','complete_verified_full_entity_operator_bank'),
        ('design_inventory','metadata/full_expanded_model_designs_v2_completed_20261002.json','complete_verified_full_expanded_model_designs')]}
    result=dict(status='verified_original_handles_and_full_stage_closure_checkpoint',
        checked_utc=datetime.now(timezone.utc).isoformat(),previous_checkpoint=str(args.previous),
        previous_checkpoint_sha256=sha(args.previous),source_plan=previous['source_plan'],
        source_plan_sha256=previous['source_plan_sha256'],pins_verified=checked_pins,
        source_hashes={str(Path(__file__)):sha(Path(__file__)),
            'scripts/record_project_runtime_checkpoint_v4.py':sha('scripts/record_project_runtime_checkpoint_v4.py')},
        handles=handles,new_software_and_selection_proofs=software,stages=stages,
        all_eight_scientific_aims_incomplete=True,gpu_inference_paused=True,
        scope='Exact original identities and immutable plans retained. Completed handles require '
              'actual invocation-linked completion/resource journals; artifact acceptance is '
              'separately recorded by full stage closures. No settings changes or restarts.')
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(status=result['status'],handles=len(handles),
        stages={k:v['accepted'] for k,v in stages.items()})))


if __name__=='__main__':main()
