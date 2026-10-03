#!/usr/bin/env python3
"""Close coverage computation using full hashes and original stage execution evidence."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import subprocess

import psutil

from audit_selected_taxon_identity_snapshot_v2 import sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--transports',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--archive',type=Path,required=True)
    a=p.parse_args();assert not a.output.exists() and not a.archive.exists()
    waits=json.loads(a.transports.read_text());assert set(waits['stages'])=={'software','producer','readback'}
    bindings={str(a.transports):sha(a.transports),str(Path(__file__)):sha(__file__)}
    services=[];receipts={}
    def bind(path,digest=None):
        actual=sha(path);assert digest is None or actual==digest,path
        assert str(path) not in bindings or bindings[str(path)]==actual,path
        bindings[str(path)]=actual
    for stage,transport in waits['stages'].items():
        assert transport['original_wait_tool_exit_code']==0 and transport['original_wait_session_id']>0
        epath=Path(transport['execution']);bind(epath,transport['execution_sha256'])
        e=json.loads(epath.read_text());assert e['stage']==stage and e['exit_code']==0 and not e['timed_out']
        assert e['status']=='exited_zero_with_receipt'
        assert e['actual_cgroup_limits']=={'cpu.max':'200000 100000','memory.max':str(8*2**30),'memory.swap.max':'0'}
        assert e['invocation_id']==transport['invocation_id']
        for person in ['wrapper','child']:
            identity=e[person]
            try:
                proc=psutil.Process(identity['pid'])
                assert proc.create_time()!=identity['created'] or proc.status()==psutil.STATUS_ZOMBIE
            except psutil.NoSuchProcess:pass
        receipt=Path(e['receipt']);bind(receipt,e['receipt_sha256']);r=json.loads(receipt.read_text())
        receipts[stage]=r
        for record in [e,r]:
            for field in ['source_hashes','artifacts']:
                for path,digest in record.get(field,{}).items():bind(path,digest)
        records=[json.loads(line) for line in subprocess.check_output(
            ['journalctl','--user','--unit',transport['unit'],'--no-pager','-o','json'],text=True,timeout=30).splitlines()]
        matched=[row for row in records if row.get('USER_INVOCATION_ID',row.get('_SYSTEMD_INVOCATION_ID'))==e['invocation_id']]
        started=[row for row in matched if row.get('JOB_RESULT')=='done' and row.get('MESSAGE','').startswith('Started ')]
        resource=[row for row in matched if 'CPU_USAGE_NSEC' in row and 'MEMORY_SWAP_PEAK' in row]
        assert started and resource and all(int(row['CPU_USAGE_NSEC'])>0 for row in resource)
        assert all(int(row['MEMORY_SWAP_PEAK'])==0 for row in resource)
        services.append(dict(stage=stage,unit=transport['unit'],invocation_id=e['invocation_id'],
            original_wrapper=e['wrapper'],original_child=e['child'],original_wait_tool_exit_code=0,
            captured_invocation_journal_records=matched,actual_cgroup_limits=e['actual_cgroup_limits'],
            child_cpu_seconds=e['child_cpu_seconds'],child_peak_rss_bytes=e['child_peak_rss_bytes'],wall_seconds=e['wall_seconds'],
            memory_scope='Raw manager MEMORY_PEAK preserved as reported; no native/group peak requirement inferred. Child RSS is separately measured.'))
    assert receipts['software']['synthetic_retained_sets']==512
    assert receipts['producer']['projection_cases']==receipts['readback']['projection_cases']==17500
    assert receipts['producer']['branch_projection_cells']==receipts['readback']['branch_projection_cells']==9047500
    assert receipts['readback']['status']=='passed_full_structural_marker_raw_tree_projection_readback'
    assert receipts['producer']['status_counts']==receipts['readback']['status_counts']
    archive=dict(source_hashes=bindings,original_stage_services=services,
        scope='Three actual original wait-tool completions, saved original child exits/identities, actual cgroup limits, full source/artifact hashes and matching original manager invocation startup/resource records. Collected systemd default-success properties are not used as completion proof.')
    with a.archive.open('x') as f:f.write(json.dumps(archive,indent=2)+'\n')
    for path,digest in bindings.items():assert sha(path)==digest,path
    result=dict(status='complete_verified_full_structural_marker_tree_projection',checked_utc=datetime.now(timezone.utc).isoformat(),
        tree_views=70,cohorts=5,sources=2,marker_slots=125,full_panel_entries=526,full_panel_outgroups=25,
        original_view_internal_branches=36190,projection_cases=17500,branch_projection_cells=9047500,
        status_counts=receipts['readback']['status_counts'],source_bindings=len(bindings),
        original_stage_execution_proofs=3,full_hash_archive=str(a.archive),full_hash_archive_sha256=sha(a.archive),
        scientific_eligibility=False,
        scope='Complete full observed-marker coverage projection and independent actual raw-tree/pruning/path-length readback, with all original internal branches, all125slots, both sources and all70candidateviews retained. Coverage geometry only: no fitted structural/sequence rates, accepted topology/root/model, statistical estimability, posterior qualification or biological aim completion. No native/GPU restart, prediction or charges.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
