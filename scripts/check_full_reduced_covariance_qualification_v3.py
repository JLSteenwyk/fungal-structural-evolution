#!/usr/bin/env python3
"""Complete alternate-source contracts with byte-exact fixture restoration."""
import argparse
from copy import deepcopy
from datetime import datetime, timezone
import itertools
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from check_full_reduced_covariance_qualification import fixture, full_fixture, closed, write, rejected
from full_reduced_covariance_sources_v2 import load
from reference_measurement_union_sources import verify
from run_full_reduced_covariance_qualification_v2 import run, SUMMARY_FIELDS
from run_parallel_covariance_readback import STATUS as PARALLEL_STATUS


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--receipt',type=Path,required=True);a=p.parse_args();root=a.output.resolve();root.mkdir(exist_ok=False)
    assert not a.receipt.exists()
    fixtures={key:fixture(key[0],key[1],100+int(key[0])) for key in itertools.product([False,True],['signed','unsigned'])}
    trees=['mafft_guide','pmsf_mafft_profile','pmsf_profile_mafft','pmsf_profile_profile','profile_guide']
    fixture_root=root/'fixture';fixture_root.mkdir();baseline_path,_=full_fixture(fixture_root,fixtures,trees)
    baseline=json.loads(baseline_path.read_text());oldplan=Path(baseline['qualification_plan']);oldroot=Path(json.loads(oldplan.read_text())['output'])
    expected=baseline['expected'];old_summary=dict(**expected,trees=trees,loading_modes=['signed','unsigned'])
    producer=dict(status='complete_full_uniform_covariance_qualification_pending_independent_readback',plan_sha256=sha(oldplan),
        source_contract='explicit-synthetic-original-contract',scientific_eligibility=False,**old_summary)
    (oldroot/'receipt.json').write_text(json.dumps(producer)+'\n')
    alternate=root/'parallel';alternate.mkdir();parallel_plan=root/'parallel-plan.json'
    write(parallel_plan,dict(source_plan=str(oldplan),output=str(alternate),resources=dict(workers=8,maximum_pending_cohorts=16)))
    independent=alternate/'readback.json'
    proof=dict(status=PARALLEL_STATUS,plan_sha256=sha(oldplan),parallel_plan_sha256=sha(parallel_plan),
        producer_receipt_sha256=sha(oldroot/'receipt.json'),source_contract=producer['source_contract'],parallel_workers=8,
        maximum_observed_pending_cohorts=4,actual_cgroup_limits={'cpu.max':'800000 100000','memory.max':str(64*2**30),'memory.swap.max':'0'},
        scientific_eligibility=False,**old_summary)
    write(independent,proof)
    completion=closed(root,'alternate','complete_verified_full_uniform_covariance_qualification',
        [oldplan,parallel_plan,independent,*[oldroot/n for n in ['receipt.json','stage_plan.json','design_covariance_audits.jsonl.gz','setting_audit_links.tsv.gz']]],old_summary)
    compact=json.loads(completion.read_text());compact.update(independent_readback=str(independent),independent_readback_sha256=sha(independent),
        producer_receipt=str(oldroot/'receipt.json'),producer_receipt_sha256=sha(oldroot/'receipt.json'))
    completion.write_text(json.dumps(compact,indent=2)+'\n')
    # The original reader path is absent, exactly as it can be while the real
    # serial controller is live. Preserve its fixture bytes in another file.
    original_reader=oldroot/'readback.json';preserved=root/'preserved-original-reader-fixture.json';original_reader.rename(preserved)
    plan=dict(baseline,qualification_completion=str(completion),parallel_readback_plan=str(parallel_plan),output=str(root/'qualified'))
    pp=root/'new-plan.json';write(pp,plan)
    produced=run(pp);read=run(pp,reader=True)
    assert all(produced[k]==read[k] for k in SUMMARY_FIELDS)
    assert produced['audit_rows']==600 and produced['setting_audit_links']==1200
    assert produced['retained_basis_audit_counts']=={'4':300,'5':300}
    assert produced['link_status_counts']['constant_response']==600
    assert produced['audit_status_counts']=={'qualified_exact_retained_uniform_covariance_basis':600}
    assert not original_reader.exists() and preserved.exists()
    positive_snapshot={p:p.read_bytes() for p in [independent,completion,Path(compact['full_hash_archive'])]}
    negatives=[]
    rejected(lambda:run(pp));negatives.append('completed_producer_restart')
    rejected(lambda:run(pp,reader=True));negatives.append('completed_reader_restart')
    def replace_reader(change):
        altered=deepcopy(proof);change(altered);independent.write_text(json.dumps(altered,indent=2)+'\n')
        # Rehash the synthetic closed archive so malformed identities reach
        # the adapter's semantic checks rather than only its hash check.
        archive_path=Path(compact['full_hash_archive']);archive=json.loads(archive_path.read_text())
        archive['source_hashes'][str(independent)]=sha(independent);archive_path.write_text(json.dumps(archive,indent=2)+'\n')
        altered_compact=dict(compact,full_hash_archive_sha256=sha(archive_path),independent_readback_sha256=sha(independent))
        completion.write_text(json.dumps(altered_compact,indent=2)+'\n')
        rejected(lambda:load(plan,pp))
        for p,raw in positive_snapshot.items():p.write_bytes(raw)
    for name,change in [
        ('unclosed_parallel_reader',lambda r:r.__setitem__('status','pending')),
        ('wrong_original_plan',lambda r:r.__setitem__('plan_sha256','0'*64)),
        ('wrong_parallel_plan',lambda r:r.__setitem__('parallel_plan_sha256','0'*64)),
        ('wrong_producer_receipt',lambda r:r.__setitem__('producer_receipt_sha256','0'*64)),
        ('changed_source_contract',lambda r:r.__setitem__('source_contract','foreign')),
        ('missing_original_settings',lambda r:r.__setitem__('model_setting_rows',119)),
        ('missing_tree',lambda r:r.__setitem__('trees',r['trees'][:-1])),
        ('missing_loading_mode',lambda r:r.__setitem__('loading_modes',['signed'])),
        ('wrong_worker_count',lambda r:r.__setitem__('parallel_workers',1)),
        ('unbounded_pending_cohorts',lambda r:r.__setitem__('maximum_observed_pending_cohorts',17)),
        ('changed_live_resource_caps',lambda r:r['actual_cgroup_limits'].__setitem__('memory.max',str(32*2**30))),
        ('promoted_scientific_claim',lambda r:r.__setitem__('scientific_eligibility',True))]:
        replace_reader(change);negatives.append(name)
    completion.write_text(json.dumps({**compact,'status':'pending_original_journals'}))
    rejected(lambda:load(plan,pp));negatives.append('unclosed_alternate_source');completion.write_bytes(positive_snapshot[completion])
    altered=dict(plan,parallel_readback_plan=str(oldplan));rejected(lambda:load(altered,pp));negatives.append('wrong_reader_namespace')
    assert all(p.read_bytes()==raw for p,raw in positive_snapshot.items())
    for name in ['receipt.json','readback.json']:
        verify(json.loads((root/'qualified'/name).read_text())['source_hashes'])
    bindings={str(p):sha(p) for p in [Path(__file__),Path('scripts/full_reduced_covariance_sources_v2.py'),
        Path('scripts/run_full_reduced_covariance_qualification_v2.py'),Path('scripts/reduced_covariance_basis.py')]}
    for p in root.rglob('*'):
        if p.is_file():bindings[str(p)]=sha(p)
    verify(bindings)
    result=dict(status='passed_full_exact_retained_covariance_parallel_source_v3_software_contracts',checked_utc=datetime.now(timezone.utc).isoformat(),
        synthetic_pipeline_audits=600,synthetic_pipeline_setting_links=1200,retained_basis_audit_counts={'4':300,'5':300},
        all_five_trees_and_both_modes=True,inherited_numerical_envelopes_unchanged=True,constant_response_links_retained=600,
        alternate_source_cases_rejected=negatives,original_reader_path_not_required_or_written=True,byte_exact_post_tampering_fixture_restoration=True,complete_positive_source_bindings_rehashed=True,
        source_and_journal_fixtures_synthetic=True,scientific_eligibility=False,source_hashes=bindings,
        scope='Complete600audit/1200link retained-basis producer/readback after explicitly synthetic alternate original-source closure. Both four/five kernels, both modes/allfive trees/constant states retained. Original reader path absent and never recreated;16 rehashed/malformed alternate identities, source states and restarts rejected. Frozen reducer/gesvd/error envelopes unchanged. No production source closure, fungal pilot, fit or accepted biological effect.')
    write(a.receipt,result);print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
