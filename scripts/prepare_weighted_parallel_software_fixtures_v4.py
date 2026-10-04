#!/usr/bin/env python3
"""Prepare complete private synthetic source fixtures under their original caps."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from check_full_weighted_covariance_qualification_v2 import fixture
from prepare_baliphy_scalar_v6_preflight_v1 import project_sources
from reference_measurement_union_sources import bind, verify

AUDIT_CASES = ['raw_gram', 'projected_gram', 'raw_envelope', 'projected_envelope',
    'diagnostic_rank', 'retained_names', 'uniform_class', 'control_label',
    'diagonal_hash', 'source_design_status', 'promote_science', 'omit_last_audit',
    'duplicate_audit', 'foreign_record_count']
LINK_CASES = ['foreign_link', 'missing_link', 'duplicate_link',
    'changed_original_setting', 'changed_setting_ordinal']


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    a=p.parse_args(); assert not a.receipt.exists(); a.output.mkdir(exist_ok=False)
    pins={}; modules=project_sources(pins,[Path(__file__),
        Path('scripts/check_weighted_parallel_software_v2.py'),
        Path('scripts/full_weighted_covariance_qualification_parallel_v1.py')])
    grids=[]; negative=[]
    for name,qualified,exception in [('unqualified',False,False),
            ('qualified-common-pair',True,False),('qualified-pair-exception',True,True)]:
        root=a.output/name
        original=fixture(root,qualified,exception)
        spec=json.loads(original.read_text())
        base=dict(spec, pins={**spec['pins'],**pins}, resources=dict(cpus=2,memory_gib=32,
            workers=2,worker_address_space_gib=12,reservation_capacity_gib=24,
            parent_headroom_gib=8,minimum_free_disk_gib=128))
        paths={}
        for version in ['serial','parallel']:
            plan=dict(base,output=str(root/(version+'-numerical')))
            path=root/(version+'.plan.json')
            with path.open('x') as f: json.dump(plan,f,indent=2);f.write('\n')
            paths[version]=str(path);bind(pins,path)
        grids.append(dict(name=name,qualified=qualified,pair_exception=exception,**paths))
        if name=='unqualified':
            for case in AUDIT_CASES+LINK_CASES:
                plan=dict(base,output=str(root/('negative-'+case)))
                path=root/('negative-'+case+'.plan.json')
                with path.open('x') as f:json.dump(plan,f,indent=2);f.write('\n')
                negative.append(dict(case=case,plan=str(path)));bind(pins,path)
    verify(pins)
    result=dict(status='prepared_three_complete_private_parallel_software_grids_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(),grids=grids,
        negatives=negative,project_modules=sorted(map(str,modules)),source_hashes=pins,
        synthetic_source_closures=True,scientific_eligibility=False,
        scope='Three complete synthetic source grids constructed and census/readback executed '
              'under their original two-CPU/16-GiB caps. Numerical plans are prepared before '
              'execution for two-CPU/32-GiB serial/parallel qualification. Separate output '
              'namespaces retain all private corruption cases. No biological pilot or production fit.')
    with a.receipt.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','project_modules']},indent=2))


if __name__=='__main__':main()
