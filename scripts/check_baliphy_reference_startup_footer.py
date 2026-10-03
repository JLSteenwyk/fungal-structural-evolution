#!/usr/bin/env python3
"""Check real timing-footer startup records and reject unrecognized stdout."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from baliphy_reference_startup_readback_v2 import inspect_v2, parse_initial_model


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); np = Path('metadata/baliphy_reference_preflight_plan_20261003.json')
    plan = json.loads(np.read_text()); jobs = json.loads(Path(plan['jobs']).read_text())
    root = Path(plan['output']); actual = []; bindings = {str(np): sha(np)}; first_footer = None
    for job in jobs:
        cp = root / 'chains' / (job['chain']['chain_id'] + '.json')
        if not cp.exists(): continue
        saved = json.loads(cp.read_text()); receipt = Path(saved['attempt_receipt'])
        row = inspect_v2(job, receipt, sha(np), 'software-footer-fixture')
        assert row['status'] == 'reference_startup_homology_density_and_representation_checked', row
        has_footer = row['stdout_timing_summary'] is not None
        if has_footer and first_footer is None: first_footer = receipt.parent / 'stdout.log'
        actual.append(dict(chain_id=row['chain_id'], original_disposition=saved['status'],
            revised_status=row['status'], has_footer=has_footer, receipt=str(receipt), receipt_sha256=sha(receipt)))
        for name, h in json.loads(receipt.read_text())['artifacts'].items(): bindings[str(receipt.parent / name)] = h
        bindings[str(receipt)] = sha(receipt); bindings[str(cp)] = sha(cp)
        if len(actual) >= 24: break
    assert len(actual) >= 16 and first_footer is not None
    text = first_footer.read_text(); value, end = json.JSONDecoder().raw_decode(text.lstrip())
    initial = text.lstrip()[:end]; footer = text.lstrip()[end:].strip()
    rejected = []
    for name, changed in [
        ('second_json', initial + '\n{}'), ('unknown_suffix', initial + '\nerror: unexpected'),
        ('footer_then_message', text + '\nextra output'), ('incomplete_footer', initial + '\nWork:'),
        ('negative_cpu', text.replace('(CPU) time:', '(CPU) time: -1s (-1s) #')),
        ('invalid_date', initial + '\n' + footer.replace('start:', 'start: invalid ')),
        ('changed_footer_label', text.replace('Work:', 'Unknown:')),
    ]:
        try: parse_initial_model(changed)
        except (AssertionError, KeyError, ValueError): rejected.append(name)
        else: raise AssertionError('Unexpected stdout was accepted: ' + name)
    assert parse_initial_model(initial + '\n\n')[1] is None
    for path in ['scripts/baliphy_reference_startup_readback_v2.py', __file__,
                 'scripts/run_baliphy_reference_preflight.py', 'scripts/baliphy_reference_initialization.py']:
        bindings[str(path)] = sha(path)
    result = dict(status='passed_actual_native_reference_startup_footer_contracts',
        checked_utc=datetime.now(timezone.utc).isoformat(), actual_full_grid_startups_checked=len(actual),
        actual_timing_footer_records_checked=sum(x['has_footer'] for x in actual),
        malformed_stdout_rejected=rejected, actual_startups=actual, source_hashes=bindings,
        scientific_eligibility=False,
        scope='Real preserved native preflight outputs and seven malformed stdout contracts. Only the exact native Work footer accepted; extant homology, fixed-tip/free-ancestor lengths, same original model source and degree-aware density rechecked. Full1620footer replay pending; no native reruns or acceptance of original allocation failures.')
    with a.output.open('x') as f: f.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['actual_startups', 'source_hashes']}), flush=True)


if __name__ == '__main__': main()
