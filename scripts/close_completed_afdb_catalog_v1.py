#!/usr/bin/env python3
"""Close the original full catalog/reader and publish fixed-universe coverage."""
import argparse
import ast
from collections import Counter
import csv
from datetime import datetime,timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind,verify


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt',type=Path,required=True)
    args=parser.parse_args();assert not args.receipt.exists()
    plan_path=Path('metadata/completed_afdb_catalog_plan_20261005_v1.json')
    plan=json.loads(plan_path.read_text());root=Path(plan['output'])
    raw_path=root/'receipt.json';producer_path=Path('metadata/completed_afdb_catalog_20261005_v1.json')
    reader_path=Path('metadata/completed_afdb_catalog_readback_20261005_v1.json')
    raw,producer,reader=[json.loads(path.read_text()) for path in [raw_path,producer_path,reader_path]]
    assert producer['status']=='completed_full_afdb_catalog_pending_independent_readback'
    assert reader['status']=='passed_full_proteome_sequence_and_model_selection_readback'
    assert producer['raw_catalog_receipt_sha256']==reader['producer_receipt_sha256']==sha(raw_path)
    assert producer['taxa']==reader['taxa']==raw['taxa']==526
    assert producer['proteins_screened']==reader['proteins_screened']==raw['proteins_screened']==5815847
    assert producer['proteins_linked']==reader['protein_links']==raw['proteins_linked']
    assert producer['unique_models']==reader['models']==raw['unique_models']
    pins={};transports=[]
    for prefix,session,unit,validation,queued in [
        ('completed_afdb_catalog',20607,'fungal-completed-afdb-catalog-20261005-v1.service',producer_path,False),
        ('completed_afdb_catalog_readback',32067,'fungal-completed-afdb-catalog-readback-20261005-v1.service',reader_path,True)]:
        ep=Path(f'metadata/{prefix}_execution_20261005_v1.json')
        tp=Path(f'metadata/{prefix}_transport_20261005_v1.json')
        pp=Path(f'metadata/{prefix}_original_tool_payloads_20261005_v1.json')
        lp=Path(f'metadata/{prefix}_launch_20261005_v1.json')
        e,t,tool,launch=[json.loads(path.read_text()) for path in [ep,tp,pp,lp]]
        assert e['status']=='exited_zero_with_receipt' and e['exit_code']==0 and not e['timed_out']
        assert e['receipt_sha256']==t['validation_sha256']==sha(validation)
        assert t['original_tool_session_id']==tool['original_tool_session_id']==launch['original_tool_session_id']==session
        assert tool['initial']['session_id']==session and tool['terminal']['exit_code']==0
        assert e['invocation_id']==launch['invocation_id']
        assert t['whole_wrapper_initial_and_terminal_payloads_matched']
        assert t['manager_start_records']==t['manager_completion_records']==1
        for mapping in [e['source_hashes'],e['artifacts'],t['source_hashes']]:
            for path,digest in mapping.items():bind(pins,path,digest)
        for path in [ep,tp,pp,lp]:bind(pins,path)
        text=subprocess.check_output(['journalctl','--user','-u',unit,'--all','-o','json','--no-pager'],text=True)
        rows=[json.loads(line) for line in text.splitlines()]
        inv=e['invocation_id'];rows=[r for r in rows if inv in [r.get('_SYSTEMD_INVOCATION_ID'),r.get('USER_INVOCATION_ID')]]
        exact=[r for r in rows if r.get('_PID')==str(e['wrapper']['pid']) and r.get('_CMDLINE')==' '.join(e['wrapper']['cmdline'])]
        assert len(exact)==2
        assert json.loads(exact[0]['MESSAGE'])==dict(original_wrapper=e['wrapper'],invocation_id=inv)
        assert json.loads(exact[1]['MESSAGE'])=={k:v for k,v in e.items() if k not in ['source_hashes','artifacts','command','wrapper','child','scope']}
        assert not any('Main process exited' in r.get('MESSAGE','') or 'Failed with result' in r.get('MESSAGE','') for r in rows)
        if queued:
            outer=[r for r in rows if r.get('_PID')==str(launch['pid']) and r.get('_CMDLINE')==' '.join(launch['cmdline'])]
            starts=[r for r in outer if r.get('MESSAGE','').startswith('running_original_journal_verified_command ')]
            assert len(starts)==1
            wait=json.loads(Path(launch['plan']).read_text())
            command=ast.literal_eval(starts[0]['MESSAGE'].split(' ',1)[1])
            assert command==wait['command']==e['wrapper']['cmdline']
            assert any(r.get('MESSAGE','').startswith('waiting_for_original_dependency_completion ') for r in outer)
        else:
            assert e['wrapper']==dict(pid=launch['pid'],created=launch['created'],cmdline=launch['cmdline'])
        jp=Path(f'metadata/{prefix}_original_whole_journal_20261005_v1.jsonl')
        with jp.open('x') as handle:
            for row in rows:handle.write(json.dumps(row,sort_keys=True)+'\n')
        bind(pins,jp)
        transports.append(dict(unit=unit,invocation_id=inv,original_tool_session_id=session,
            native_exit_code=0,wrapper_exit_code=0,queued_outer_controller_verified=queued,
            full_original_payloads_verified=True))
    for mapping in [producer['source_hashes'],reader['source_hashes'],plan['pins']]:
        for path,digest in mapping.items():bind(pins,path,digest)
    for path in [raw_path,producer_path,reader_path,plan_path,Path(__file__)]:bind(pins,path)
    with (root/'taxon_coverage.tsv').open() as handle:coverage=list(csv.DictReader(handle,delimiter='\t'))
    with Path(plan['manifest']).open() as handle:manifest={r['taxon_id']:r for r in csv.DictReader(handle,delimiter='\t')}
    assert len(coverage)==len({r['taxon_id'] for r in coverage})==len(manifest)==526
    summary={};unlinked=0
    for row in coverage:
        taxon=manifest[row['taxon_id']];assert row['study_role']==taxon['study_role']
        assert row['species_name']==taxon['species_name']
        n,k,m=[int(row[key]) for key in ['representative_proteins','proteins_with_model','proteins_without_catalog_model']]
        assert n>0 and k>=0 and m>=0 and k+m==n
        role=row['study_role'];summary.setdefault(role,Counter())
        summary[role].update(taxa=1,taxa_with_models=int(k>0),representative_proteins=n,
                             proteins_with_model=k,proteins_without_catalog_model=m)
        unlinked+=m
    assert sum(s['taxa'] for s in summary.values())==526
    assert sum(s['proteins_with_model'] for s in summary.values())==raw['proteins_linked']
    assert sum(s['representative_proteins'] for s in summary.values())==5815847
    assert sum(s['taxa_with_models'] for s in summary.values())==raw['taxa_with_models']
    verify(pins)
    result=dict(status='complete_verified_full_completed_retrieval_catalog_refresh',
        checked_utc=datetime.now(timezone.utc).isoformat(),taxa=526,proteins_screened=5815847,
        proteins_linked=raw['proteins_linked'],unique_models=raw['unique_models'],
        proteins_without_catalog_model=unlinked,taxa_with_models=raw['taxa_with_models'],
        coordinate_bytes_verified=raw['coordinate_bytes_verified'],inventory_records=raw['inventory_records'],
        study_role_coverage={role:dict(counts) for role,counts in summary.items()},
        original_transports=transports,complete_bound_files=len(pins),source_hashes=pins,
        scientific_eligibility=False,confidence_qualified_atlas_complete=False,
        all_eight_aims_incomplete=True,new_predictions=0,gpu=False,
        scope='Full526taxon/5,815,847representative-protein catalog and independent complete source/link/model '
              'readback closed against both actual original waits, whole wrapper/native/outer-queue/manager '
              'journals and all bindings. Every selected coordinate byte hash was checked by producer; reader '
              'independently reconstructs latest status, source-specific ranking, all links/models/taxon counts '
              'without rehashing every coordinate. Coverage is snapshot/source availability before confidence/PAE, '
              'not public-database absence, experimental validation, ESMFold union, clustering, homology or '
              'completed structural/evolutionary atlas. Missing models remain explicit for every original taxon.')
    with args.receipt.open('x') as handle:json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
