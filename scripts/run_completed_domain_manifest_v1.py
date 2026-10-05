#!/usr/bin/env python3
"""Run unchanged complete domain-boundary manifest from the closed registry."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True)
    a=p.parse_args();assert not a.receipt.exists()
    plan=json.loads(a.plan.read_text());pins=dict(plan['pins']);verify(pins)
    completed=json.loads(Path(plan['registry_closure']).read_text())
    assert completed['status']=='complete_verified_full_refreshed_afdb_domain_registry'
    assert completed['counts']['models']==2935733 and completed['counts']['protein_links']==2994868
    subprocess.run([sys.executable,'scripts/prepare_domain_extraction_manifest_v2.py',
        '--registry',plan['registry'],'--audit',plan['registry_readback'],'--output',plan['output']],check=True)
    root=Path(plan['output']);raw_path=root/'receipt.json';raw=json.loads(raw_path.read_text())
    assert raw['status']=='complete_all_candidate_domain_extraction_manifest'
    assert raw['registry_readback_sha256']==sha(plan['registry_readback'])
    for name,digest in raw['artifacts'].items():bind(pins,root/name,digest)
    for path in [a.plan,raw_path,Path(__file__)]:bind(pins,path)
    verify(pins)
    result=dict(status='completed_full_refreshed_domain_manifest_pending_independent_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(),source_models=raw['source_models'],
        unique_intervals=raw['unique_intervals'],boundary_links=raw['boundary_links'],
        candidate_model_hit_pairs=raw['candidate_model_hit_pairs'],unique_interval_residues=raw['unique_interval_residues'],
        raw_receipt_sha256=sha(raw_path),source_hashes=pins,scientific_eligibility=False,gpu=False,new_predictions=0,
        scope='Unchanged union of all four qualified Domain policy hits and both alignment/envelope '
              'boundaries over the complete closed current AFDB registry. Every association retained '
              'with reversible sequence/span deduplication. Whole independent SQL union and actual '
              'original completion required before extraction. No validated boundaries, residue '
              'confidence/PAE or evolutionary event acceptance.')
    with a.receipt.open('x') as h:json.dump(result,h,indent=2,allow_nan=False);h.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
