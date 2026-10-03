#!/usr/bin/env python3
"""Exercise the complete 1,620-role generator and reject altered full-grid designs."""
import argparse
import copy
from datetime import datetime, timezone
import json
from pathlib import Path
import tempfile

from ancestral_chain_attempt import sha
from baliphy_reference_initialization import transform
from prepare_baliphy_reference_preflight import build
from run_baliphy_reference_preflight import summarize


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    cp = Path('results/ancestral/baliphy-independent-chain-inputs-20260927-v1/chain_inputs.json')
    qp = Path('results/ancestral/full-baliphy-horizon-resources-20261003-v1/proposed.jsonl')
    chains = json.loads(cp.read_text()); proposed = [json.loads(x) for x in qp.read_text().splitlines()]
    binary = Path('data/software_audits/baliphy-4.3-20260927/install/bali-phy-4.3/bin/bali-phy').resolve()
    rejected = []
    with tempfile.TemporaryDirectory(prefix='reference-startup-grid-') as folder:
        root = Path(folder); first = root / 'all-inputs'; first.mkdir()
        jobs, models, bindings = build(chains, proposed, first, binary, Path('/usr/bin/prlimit'), [])
        assert len(jobs) == 1620 and len(models) == 405
        assert all(j['config']['command'][-3:] == ['--test', '--log-format', 'json'] for j in jobs)
        assert all('--iterations' not in j['config']['command'] for j in jobs)
        assert all(Path(j['program']).read_text() == transform(Path(j['chain']['program']).read_text()) for j in jobs)
        assert all(sha(name) == h for name, h in bindings.items())
        for name in ['missing_chain', 'duplicate_chain', 'duplicate_fresh_seed', 'source_seed_reused',
                     'invalid_seed', 'alias_changed', 'role_changed', 'prior_changed',
                     'matrix_changed_within_quartet', 'source_program_hash_changed']:
            cs, qs = copy.deepcopy(chains), copy.deepcopy(proposed)
            if name == 'missing_chain': cs.pop()
            elif name == 'duplicate_chain': cs[1] = copy.deepcopy(cs[0])
            elif name == 'duplicate_fresh_seed': qs[1]['fresh_seed'] = qs[0]['fresh_seed']
            elif name == 'source_seed_reused': qs[0]['fresh_seed'] = cs[0]['seed']
            elif name == 'invalid_seed': qs[0]['fresh_seed'] = 0
            elif name == 'alias_changed': qs[0]['original_configuration_ids'] = ['invented-alias']
            elif name == 'role_changed': cs[0]['chain'] = 4
            elif name == 'prior_changed': cs[0]['prior_label'] = 'invented-prior'
            elif name == 'matrix_changed_within_quartet': cs[1]['alignment'] = 'invented-alignment.faa'
            else: cs[0]['program_sha256'] = '0' * 64
            out = root / name; out.mkdir()
            try: build(cs, qs, out, binary, Path('/usr/bin/prlimit'), [])
            except AssertionError: rejected.append(name)
            else: raise AssertionError('Altered complete design was accepted: ' + name)
        rows = []
        failed_ids = {chains[0]['chain_id'], chains[4]['chain_id']}
        for c in chains:
            rows.append(dict(chain_id=c['chain_id'], effective_input_group=c['effective_input_group'],
                model_input_identity=c['effective_input_group'] + '-' + c['prior_label'],
                chain_role=c['chain'], original_configuration_ids=c['original_configuration_ids'],
                status=('unsuccessful_native_startup_retained' if c['chain_id'] in failed_ids
                        else 'reference_startup_homology_density_and_representation_checked')))
        summary = summarize(rows)
        assert summary['full_chains'] == 1620 and summary['full_quartets'] == 405
        assert summary['validated_startups'] == 1618 and summary['unsuccessful_startups'] == 2
        assert summary['complete_startup_quartets'] == 403 and summary['unresolved_startup_quartets'] == 2
        assert summary['effective_inputs'] == 135 and summary['original_configuration_aliases'] == 324
        try: summarize(rows[:-1])
        except AssertionError: rejected.append('missing_serialized_disposition')
        else: raise AssertionError('Dropped startup disposition accepted')
    sources = [Path(__file__), Path('scripts/baliphy_reference_initialization.py'),
               Path('scripts/prepare_baliphy_reference_preflight.py'),
               Path('scripts/run_baliphy_reference_preflight.py'),
               Path('scripts/readback_baliphy_reference_preflight.py'), cp, qp]
    result = dict(status='passed_full_reference_startup_grid_software_contracts',
        checked_utc=datetime.now(timezone.utc).isoformat(), full_chains=1620, full_quartets=405,
        effective_inputs=135, original_configuration_aliases=324,
        complete_generator_and_all_programs_exercised=True,
        simulated_failed_chains_retained=2, simulated_unresolved_quartets=2,
        altered_designs_and_missing_dispositions_rejected=rejected,
        source_hashes={str(x): sha(x) for x in sources}, scientific_eligibility=False,
        scope='Complete real input generator, exact reversible transformations and full design mutations. Artificial status rows test failure accounting only; no native startup, density, readback or production journal proof claimed here. Native software fixtures and actual full startup/readback remain separate gates.')
    with a.output.open('x') as f: f.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__': main()
