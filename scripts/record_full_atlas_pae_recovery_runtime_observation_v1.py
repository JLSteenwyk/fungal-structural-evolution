#!/usr/bin/env python3
"""Record a read-only throughput observation for a live full-atlas PAE recovery."""
import argparse, hashlib, json
from datetime import datetime, timezone
from pathlib import Path


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--state',type=Path,required=True);p.add_argument('--plan',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True);a=p.parse_args()
    if a.receipt.exists():raise FileExistsError('Fresh receipt required')
    state=json.loads(a.state.read_text());plan=json.loads(a.plan.read_text());counts=state['counts'];attempted=counts['attempted'];elapsed=state['elapsed_seconds'];expected=state['expected_models']
    if state['stage']!='retrieving_recovery_queue' or attempted<=0 or elapsed<=0 or expected<attempted or attempted!=counts['verified']+counts['failed']:raise ValueError('Invalid live recovery state')
    rate=attempted/elapsed; remaining=expected-attempted
    result=dict(status='recorded_live_full_atlas_pae_recovery_throughput_observation',checked_utc=datetime.now(timezone.utc).isoformat(),source_state=str(a.state),source_state_sha256=sha(a.state),recovery_plan=str(a.plan),recovery_plan_sha256=sha(a.plan),attempted=attempted,verified=counts['verified'],failed=counts['failed'],elapsed_seconds=elapsed,observed_attempts_per_second=rate,remaining_models=remaining,linear_remaining_seconds=remaining/rate,linear_remaining_days=remaining/rate/86400,http_workers=state['http_workers'],future_window=state['future_window'],scope='A read-only single-time throughput observation. The linear extrapolation assumes the observed average rate persists; it is not a completion guarantee, cost estimate, coverage result, or biological inference.')
    a.receipt.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
