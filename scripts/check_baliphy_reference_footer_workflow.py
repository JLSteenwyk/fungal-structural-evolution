#!/usr/bin/env python3
"""Full-grid serialization contracts for footer replay, with explicit fake rows."""
import argparse
import copy
from datetime import datetime, timezone
import json
from pathlib import Path
import tempfile
from unittest.mock import patch

from ancestral_chain_attempt import sha
import baliphy_reference_startup_readback_v2 as replay
from run_baliphy_reference_preflight import summarize


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); native_path = Path('metadata/baliphy_reference_preflight_plan_20261003.json')
    native = json.loads(native_path.read_text()); jobs = json.loads(Path(native['jobs']).read_text()); rows = []
    for i, job in enumerate(jobs):
        c = job['chain']; failed = i in [0, 4]; footer = not failed and i < 14
        row = dict(chain_id=c['chain_id'], effective_input_group=c['effective_input_group'],
            model_input_identity=c['effective_input_group'] + '-' + c['prior_label'],
            chain_role=c['chain'], original_configuration_ids=c['original_configuration_ids'],
            fresh_seed=job['fresh_seed'], scientific_eligibility=False,
            original_startup_disposition=('unsuccessful_native_startup_retained' if failed else
                'invalid_native_startup_retained' if footer else 'reference_startup_homology_density_and_representation_checked'),
            status=('unsuccessful_native_startup_retained' if failed else
                    'reference_startup_homology_density_and_representation_checked'))
        if footer: row['stdout_timing_summary'] = dict(reported_cpu_seconds=1)
        rows.append(row)
    summary = summarize(rows)
    summary.update(native_timing_footers_checked=12, v1_parser_dispositions_reclassified=12)
    assert (summary['validated_startups'], summary['unsuccessful_startups'],
            summary['complete_startup_quartets'], summary['unresolved_startup_quartets']) == (1618, 2, 403, 2)
    rejected = []
    with tempfile.TemporaryDirectory(prefix='reference-footer-workflow-') as folder:
        root = Path(folder); plan_path = root / 'plan.json'; plan = dict(output=str(root / 'stage'), scope='Software fixture only')
        plan_path.write_text(json.dumps(plan))
        bindings = {str(native_path): sha(native_path)}
        with patch.object(replay, 'reconstruct', return_value=(plan, rows, summary, bindings)):
            producer = replay.run(plan_path); assert all(producer[k] == v for k, v in summary.items())
            try: replay.run(plan_path)
            except AssertionError: rejected.append('completed_producer_restart')
            else: raise AssertionError('Completed replay restarted')
            exported = Path(plan['output']) / 'dispositions.json'; original = exported.read_bytes()
            for name in ['missing_row', 'duplicated_row', 'changed_seed', 'removed_failed_chain',
                         'changed_alias', 'removed_timing', 'changed_source_disposition', 'scientific_acceptance']:
                bad = copy.deepcopy(rows)
                if name == 'missing_row': bad.pop()
                elif name == 'duplicated_row': bad[1] = copy.deepcopy(bad[0])
                elif name == 'changed_seed': bad[1]['fresh_seed'] += 1
                elif name == 'removed_failed_chain': bad[0]['status'] = 'reference_startup_homology_density_and_representation_checked'
                elif name == 'changed_alias': bad[1]['original_configuration_ids'] = ['invented-alias']
                elif name == 'removed_timing': bad[1].pop('stdout_timing_summary')
                elif name == 'changed_source_disposition': bad[1]['original_startup_disposition'] = 'invented-status'
                else: bad[1]['scientific_eligibility'] = True
                exported.write_text(json.dumps(bad))
                try: replay.run(plan_path, reader=True)
                except AssertionError: rejected.append(name)
                else: raise AssertionError('Altered serialized replay accepted: ' + name)
                exported.write_bytes(original)
            reader = replay.run(plan_path, reader=True)
            assert all(reader[k] == v for k, v in summary.items())
            try: replay.run(plan_path, reader=True)
            except AssertionError: rejected.append('completed_reader_restart')
            else: raise AssertionError('Completed reader restarted')
    sources = [Path(__file__), Path('scripts/baliphy_reference_startup_readback_v2.py'),
               Path('scripts/run_baliphy_reference_preflight.py'), native_path]
    result = dict(status='passed_full_reference_footer_serialization_software_contracts',
        checked_utc=datetime.now(timezone.utc).isoformat(), full_chains=1620, full_quartets=405,
        simulated_failed_chains_retained=2, simulated_footer_rows=12,
        full_producer_and_reader_serialization_exercised=True,
        altered_exports_and_completed_restarts_rejected=rejected,
        source_hashes={str(x): sha(x) for x in sources}, scientific_eligibility=False,
        scope='Full1620real metadata identities with artificial status/timing rows. Native reconstruction intentionally mocked to test complete serialization, false exports and completed restart refusal. Not a native-output, journal, likelihood or scientific proof. Actual native records and full production replay remain separate gates.')
    with a.output.open('x') as f: f.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__': main()
