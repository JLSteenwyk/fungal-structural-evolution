#!/usr/bin/env python3
"""Validate the actual entire failed-role job matrix and reject scope/resource changes."""
import argparse
import copy
from datetime import datetime,timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from baliphy_joint_fasta_v7_failure_jobs_v1 import failed_sources,build_jobs,validate_jobs
from reference_measurement_union_sources import bind,verify


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True)
    a=p.parse_args();assert not a.receipt.exists();a.output.mkdir(exist_ok=False)
    originals,pins,mapping=failed_sources();jobs=build_jobs(originals,a.output/'models')
    positive=validate_jobs(jobs,originals)
    serialized=a.output/'jobs.json'
    with serialized.open('x') as f:json.dump(jobs,f,indent=2);f.write('\n')
    assert json.loads(serialized.read_text())==jobs
    rejected=[]
    mutations=['drop_role','duplicate_role','source_role','chain_id','seed','prior','effective_input',
               'aliases','taxa','alignment_hash','tree_hash','iterations','cpu','wall','as','file','stack',
               'memory_reservation','program_hash','extra_native_argument']
    for case in mutations:
        altered=copy.deepcopy(jobs);j=altered[0];c=j['config']['command']
        if case=='drop_role':altered.pop()
        elif case=='duplicate_role':altered[-1]=copy.deepcopy(altered[0])
        elif case=='source_role':j['source_v6_chain_id']='foreign'
        elif case=='chain_id':j['chain']['chain_id']+='-foreign'
        elif case=='seed':j['chain']['seed']+=1;c[c.index('--seed')+1]=str(j['chain']['seed'])
        elif case=='prior':j['chain']['prior_label']='narrow'
        elif case=='effective_input':j['chain']['effective_input_group']='foreign'
        elif case=='aliases':j['chain']['original_configuration_ids'].pop()
        elif case=='taxa':j['chain']['proteins']-=1
        elif case=='alignment_hash':j['chain']['alignment_sha256']='0'*64
        elif case=='tree_hash':j['chain']['tree_sha256']='0'*64
        elif case=='iterations':c[c.index('--iterations')+1]='19'
        elif case in ['cpu','as','file','stack']:
            flag={'cpu':'--cpu=','as':'--as=','file':'--fsize=','stack':'--stack='}[case]
            index=next(i for i,v in enumerate(c) if v.startswith(flag));c[index]=flag+'1'
        elif case=='wall':j['config']['timeout_seconds']+=1
        elif case=='memory_reservation':j['memory_reservation_bytes']-=1
        elif case=='program_hash':j['chain']['program_sha256']='0'*64
        elif case=='extra_native_argument':c+=['--foreign']
        try:validate_jobs(altered,originals)
        except (AssertionError,KeyError,ValueError):rejected.append(case)
        else:raise AssertionError('Altered job accepted: '+case)
    assert len(rejected)==20
    for path in [Path(__file__),Path('scripts/baliphy_joint_fasta_v7_failure_jobs_v1.py'),
                 Path('scripts/baliphy_joint_fasta_lines_v7.py'),serialized,*sorted((a.output/'models').glob('*.hs'))]:bind(pins,path)
    verify(pins)
    result=dict(status='passed_exact_all24_failed_role_candidate_job_contracts_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(),scope_check=positive,
        rejected_mutations=rejected,actual_original_failed_receipts_reverified=24,
        jobs=str(serialized),mapping=mapping,source_hashes=pins,new_native_runs=0,
        original_jobs_restarted=False,scientific_eligibility=False,posterior_qualified=False,gpu=False,
        scope='Actual24failed-role source records, input/model/configuration/seed/resource binding '
              'and serialized job construction are checked, with20scope/capacity negatives. '
              'No native run or mock biological outcome; this qualifies construction only. '
              'Full first-input output closure remains a separate prerequisite for native admission.')
    with a.receipt.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
