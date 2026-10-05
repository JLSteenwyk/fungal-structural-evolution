#!/usr/bin/env python3
"""Independently reconstruct complete source-model PAE availability and the entire retrieval queue."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True)
    a=p.parse_args();assert not a.receipt.exists()
    plan=json.loads(a.plan.read_text());pins=dict(plan['pins']);verify(pins)
    producer_path,tp=Path(plan['producer_receipt']),Path(plan['producer_transport'])
    producer,t=[json.loads(fp.read_text()) for fp in (producer_path,tp)]
    assert producer['status']=='complete_full_atlas_pae_cache_and_exact_model_availability_snapshot'
    assert t['status']=='verified_original_software_wait_and_whole_wrapper_payloads'
    assert t['validation_sha256']==sha(producer_path) and t['original_tool_terminal_exit_code']==0
    assert t['whole_wrapper_initial_and_terminal_payloads_matched']
    assert t['manager_start_records']==t['manager_completion_records']==1
    for mapping in (producer['source_hashes'],t['source_hashes']):
        for path,digest in mapping.items():bind(pins,path,digest)
    verify(pins)
    root=Path(plan['inventory']);cache={}
    with (root/'cache_metadata.jsonl').open() as handle:
        for line in handle:
            entry=json.loads(line);r=entry['receipt'];identity=(r['model_id'],int(r['version']))
            assert identity not in cache and sha(entry['receipt_path'])==entry['receipt_sha256']
            assert json.loads(Path(entry['receipt_path']).read_text())==r
            assert sha(r['path'])==r['gzip_sha256'];cache[identity]=entry
    assert len(cache)==producer['cache_entries']
    counts,sources,cells=Counter(),Counter(),Counter();missing_cells=0;missing_count=0;used=set()
    with Path(plan['models']).open() as original,gzip.open(root/'model_pae_dispositions.jsonl.gz','rt') as output,(root/'missing_afdb_models.jsonl').open() as queue:
        emitted=(json.loads(line) for line in output);requested=(json.loads(line) for line in queue)
        for line in original:
            raw=json.loads(line);source,m=raw['source'],raw['model'];row=next(emitted)
            assert row['source']==source and row['model']==m
            if source=='ESMFold':
                expected='local_prediction_candidate_pending_matrix_validation';assert 'cache' not in row
                cells[source]+=m['length']*m['length']
            else:
                assert source=='AFDB';identity=(m['model_id'],m['version']);cached=cache.get(identity)
                endpoint='https://alphafold.ebi.ac.uk/files/'+m['model_id']+'-predicted_aligned_error_v'+str(m['version'])+'.json'
                if cached is None:
                    assert 'cache' not in row
                    if m.get('pae_url')==endpoint:
                        expected='missing_cached_matrix_retrievable';assert next(requested)==m
                        missing_count+=1;missing_cells+=m['length']*m['length']
                    else:expected='missing_cached_matrix_without_matching_advertised_url'
                else:
                    used.add(identity);assert row['cache']==cached;r=cached['receipt']
                    identity_matches=(r['sequence_sha256'],r['length'],r['url'])==(m['sequence_sha256'],m['length'],endpoint)
                    if identity_matches and m.get('pae_url')==endpoint:
                        expected='exact_cached_candidate_pending_matrix_validation';cells[source]+=m['length']*m['length']
                    else:expected='cache_provenance_mismatch'
            assert row['status']==expected;counts[source+':'+expected]+=1;sources[source]+=1
        assert next(emitted,None) is None and next(requested,None) is None
    assert dict(sources)==producer['sources'] and dict(counts)==producer['dispositions']
    assert dict(cells)==producer['candidate_directional_matrix_cells'] and missing_cells==producer['missing_retrievable_matrix_cells']
    assert len(used)==producer['cache_entries_matching_any_selected_model']
    for path in (a.plan,producer_path,tp,Path(__file__)):bind(pins,path)
    verify(pins)
    result=dict(status='passed_full_atlas_pae_model_and_missing_queue_reconstruction',
                checked_utc=datetime.now(timezone.utc).isoformat(),source_models=dict(sources),dispositions=dict(counts),
                candidate_directional_matrix_cells=dict(cells),missing_retrievable_afdb_models=missing_count,
                missing_retrievable_matrix_cells=missing_cells,cache_entries=len(cache),cache_entries_outside_current_model_inventory=len(cache)-len(used),
                source_hashes=pins,scientific_eligibility=False,gpu=False,new_predictions=0,
                scope='Complete independent every-source-model/cache identity/URL and ordered missing-queue reconstruction, '
                      'all availability and matrix-size totals. No producer classifier imported. Matrix content, '
                      'native context confidence, accuracy calibration and biological inference remain unqualified.')
    with a.receipt.open('x') as handle:json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
