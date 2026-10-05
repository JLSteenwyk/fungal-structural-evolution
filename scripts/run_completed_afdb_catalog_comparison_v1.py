#!/usr/bin/env python3
"""Run the unchanged full catalog comparison only after verified atlas closure."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import subprocess
import sys

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind,verify


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--receipt',type=Path,required=True)
    args=parser.parse_args();assert not args.receipt.exists()
    plan=json.loads(args.plan.read_text());pins=dict(plan['pins']);verify(pins)
    closure=json.loads(Path(plan['new_catalog_closure']).read_text())
    assert closure['status']=='complete_verified_full_completed_retrieval_catalog_refresh'
    assert closure['taxa']==526 and closure['proteins_screened']==5815847
    assert all(row['full_original_payloads_verified'] for row in closure['original_transports'])
    subprocess.run([sys.executable,'scripts/compare_whole_proteome_catalogs.py','--plan',str(args.plan)],check=True)
    root=Path(plan['output']);raw=root/'receipt.json';result=json.loads(raw.read_text())
    assert result['status']=='complete_validated_catalog_comparison' and result['plan_sha256']==sha(args.plan)
    for path,digest in result['sources'].items():bind(pins,path,digest)
    for name,digest in result['artifacts'].items():bind(pins,root/name,digest)
    for path in [args.plan,raw,Path(__file__)]:bind(pins,path)
    verify(pins)
    proof=dict(status='completed_full_afdb_catalog_comparison_pending_readback',checked_utc=datetime.now(timezone.utc).isoformat(),
        raw_comparison_receipt=str(raw),raw_comparison_receipt_sha256=sha(raw),source_plan_sha256=sha(args.plan),
        old_links=result['old_links'],new_links=result['new_links'],net_link_change=result['net_link_change'],
        dispositions=result['dispositions'],unlinked_in_both=result['unlinked_in_both'],source_hashes=pins,
        scientific_eligibility=False,new_predictions=0,gpu=False,
        scope='Unchanged qualified full catalog comparison joins every old/new taxon-protein link, '
              'checks fixed526taxon/5,815,847protein universe and retains gains/losses/replacements/unchanged '
              'models. New catalog actual original producer/reader closure mandatory. Independent full '
              'comparison replay and original stage transport closure remain required; no confidence/PAE '
              'qualification, new inference, clustering, homology or biological result.')
    with args.receipt.open('x') as handle:json.dump(proof,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps({k:v for k,v in proof.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
