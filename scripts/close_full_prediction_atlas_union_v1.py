#!/usr/bin/env python3
"""Close original full atlas union/readback and release every taxon row."""
import argparse
from collections import Counter,defaultdict
import csv
from datetime import datetime,timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind,verify


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--receipt',type=Path,required=True);p.add_argument('--taxon-table',type=Path,required=True)
    a=p.parse_args();assert not a.receipt.exists() and not a.taxon_table.exists()
    pp=Path('metadata/full_prediction_atlas_union_plan_20261005_v1.json');plan=json.loads(pp.read_text())
    producer_path=Path('metadata/full_prediction_atlas_union_20261005_v1.json')
    reader_path=Path('metadata/full_prediction_atlas_union_readback_20261005_v1.json')
    producer,reader=[json.loads(q.read_text()) for q in [producer_path,reader_path]]
    assert producer['status']=='completed_full_prediction_atlas_union_pending_independent_readback'
    assert reader['status']=='passed_full_prediction_atlas_union_independent_source_replay'
    assert reader['producer_receipt_sha256']==sha(producer_path)
    for key in ['taxa','representative_proteins','model_counts','availability','esmfold_models_with_representative_links','esmfold_models_without_representative_links']:
        assert producer[key]==reader[key],key
    assert reader['taxa']==526 and reader['representative_proteins']==5815847
    pins={};transports=[]
    for prefix,validation,session in [('full_prediction_atlas_union',producer_path,58078),('full_prediction_atlas_union_readback',reader_path,4670)]:
        ep=Path(f'metadata/{prefix}_execution_20261005_v1.json')
        tp=Path(f'metadata/{prefix}_transport_20261005_v1.json')
        tool_path=Path(f'metadata/{prefix}_original_tool_payloads_20261005_v1.json')
        lp=Path(f'metadata/{prefix}_launch_20261005_v1.json')
        e,t,tool,launch=[json.loads(q.read_text()) for q in [ep,tp,tool_path,lp]]
        assert e['exit_code']==0 and not e['timed_out']
        assert e['receipt_sha256']==t['validation_sha256']==sha(validation)
        assert tool['initial']['session_id']==tool['original_tool_session_id']==t['original_tool_session_id']==launch['original_tool_session_id']==session
        assert tool['terminal']['exit_code']==t['original_tool_terminal_exit_code']==0
        assert e['wrapper']==dict(pid=launch['pid'],created=launch['created'],cmdline=launch['cmdline'])
        assert e['invocation_id']==launch['invocation_id']==t['invocation_id']
        assert t['whole_wrapper_initial_and_terminal_payloads_matched']
        assert t['manager_start_records']==t['manager_completion_records']==1
        # Transport binds full original wrapper/native header+terminal and manager records.
        for mapping in [e['source_hashes'],e['artifacts'],t['source_hashes']]:
            for path,digest in mapping.items():bind(pins,path,digest)
        for path in [ep,tp,tool_path,lp]:bind(pins,path)
        transports.append(dict(unit=t['unit'],invocation_id=t['invocation_id'],
            original_tool_session_id=session,actual_terminal_exit_code=0,whole_original_payloads_verified=True))
    for mapping in [producer['source_hashes'],reader['source_hashes'],plan['pins']]:
        for path,digest in mapping.items():bind(pins,path,digest)
    root=Path(plan['output']);source=root/'taxon_source_coverage.tsv'
    with source.open() as h:rows=list(csv.DictReader(h,delimiter='\t'))
    with Path(plan['manifest']).open() as h:manifest={r['taxon_id']:r for r in csv.DictReader(h,delimiter='\t')}
    assert len(rows)==len({r['taxon_id'] for r in rows})==len(manifest)==526
    summary=defaultdict(Counter);availability=Counter()
    for r in rows:
        m=manifest[r['taxon_id']];assert (r['species_name'],r['study_role'])==(m['species_name'],m['study_role'])
        counts={k:int(r[k]) for k in ['afdb_only','esmfold_only','both','neither']};n=int(r['representative_proteins'])
        assert all(v>=0 for v in counts.values()) and sum(counts.values())==n>0
        modeled=n-counts['neither'];availability.update(counts)
        summary[r['study_role']].update(taxa=1,taxa_with_model=int(modeled>0),representative_proteins=n,proteins_with_any_model=modeled,**counts)
    assert dict(availability)==reader['availability']
    assert summary['ingroup']['taxa']==501 and summary['outgroup']['taxa']==25
    assert sum(s['representative_proteins'] for s in summary.values())==5815847
    for path in [pp,producer_path,reader_path,Path(__file__)]:bind(pins,path)
    verify(pins)
    with a.taxon_table.open('xb') as h:h.write(source.read_bytes())
    assert sha(a.taxon_table)==sha(source);bind(pins,a.taxon_table)
    modeled=5815847-reader['availability']['neither']
    result=dict(status='complete_verified_full_prediction_atlas_availability_union',
        checked_utc=datetime.now(timezone.utc).isoformat(),taxa=526,representative_proteins=5815847,
        model_counts=reader['model_counts'],availability=reader['availability'],proteins_with_any_model=modeled,
        any_source_coverage_fraction=modeled/5815847,study_role_coverage={k:dict(v) for k,v in summary.items()},
        esmfold_models_with_representative_links=reader['esmfold_models_with_representative_links'],
        esmfold_models_without_representative_links=reader['esmfold_models_without_representative_links'],
        esmfold_coordinate_bytes_rehashed=producer['esmfold_coordinate_bytes_rehashed'],
        esmfold_pae_bytes_rehashed=producer['esmfold_pae_bytes_rehashed'],
        public_taxon_table=str(a.taxon_table),original_transports=transports,
        complete_bound_files=len(pins),source_hashes=pins,scientific_eligibility=False,
        confidence_qualified_atlas_complete=False,all_eight_aims_incomplete=True,new_predictions=0,gpu=False,
        scope='Complete AFDB/ESMFold source model/metadata and every5,815,847representative-protein '
              'record independently reconstructed with actual original producer/reader waits/full '
              'wrapper/manager transports and all source/output bindings. Public526row coverage table '
              'is a byte-identical source copy. Predictor overlaps and source-unlinked alternative '
              'models retained without cross-source choice. Availability before residue confidence/PAE '
              'qualification, domain/structural comparisons, homology and all evolutionary aims.')
    with a.receipt.open('x') as h:json.dump(result,h,indent=2,allow_nan=False);h.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
