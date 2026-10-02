#!/usr/bin/env python3
"""Check the corrected coordinate contract against every saved quartet."""
from copy import deepcopy
from datetime import datetime, timezone
import json
from pathlib import Path
from close_baliphy_postprocessing_20261002_v2 import coordinate_count
from run_ortholog_pair_guide_comparison import sha


def main():
    categorical = Path('results/ancestral/full-categorical-diagnostics-20260928-v1/receipt.json')
    states_path = Path('results/ancestral/full-anchored-state-traces-20260928-v2/receipt.json')
    states = json.loads(states_path.read_text())['chains']
    r = json.loads(categorical.read_text())
    bindings = {str(categorical): sha(categorical), str(states_path): sha(states_path)}
    checked = original_failures = chain_coordinates = summaries = 0
    first = None
    for x in r['groups'].values():
        if x['status'] != 'categorical_reports_complete_not_posterior_qualification':
            assert x['status'] == 'unresolved_no_verified_extracted_quartet'
            continue
        child_path = Path(x['receipt']); assert sha(child_path) == x['receipt_sha256']
        child = json.loads(child_path.read_text()); mp = Path(child['manifest'])
        assert sha(mp) == child['manifest_sha256']
        m = json.loads(mp.read_text()); count = coordinate_count(m['coordinates'])
        if first is None: first = deepcopy(m['coordinates'])
        bindings[str(child_path)] = x['receipt_sha256']; bindings[str(mp)] = child['manifest_sha256']
        for c in m['chains']:
            sx = states[c['chain_id']]; sp = Path(sx['receipt'])
            assert sha(sp) == sx['receipt_sha256']
            ss = json.loads(sp.read_text())['summaries'][0]
            cp = [q for q in ss['artifacts'] if Path(q).name == 'coordinates.json']
            assert len(cp) == 1 and sha(cp[0]) == ss['artifacts'][cp[0]]
            assert json.loads(Path(cp[0]).read_text()) == m['coordinates']
            assert ss['candidate_anchor_coordinates'] == count
            bindings[str(sp)] = sx['receipt_sha256']; bindings[cp[0]] = ss['artifacts'][cp[0]]
            chain_coordinates += 1
        rp = Path(child['report']); assert sha(rp) == child['report_sha256']
        report = json.loads(rp.read_text()); bindings[str(rp)] = child['report_sha256']
        for cutoff, info in report['outputs'].items():
            sp = rp.parent / ('discard-' + cutoff) / 'summary.json'
            summary = json.loads(sp.read_text()); bindings[str(sp)] = sha(sp)
            assert summary['coordinates'] == info['coordinates'] == count
            assert sum(summary['coordinate_status_counts'].values()) == count
            assert summary['retained_samples_per_chain'] == (75 if cutoff == '250' else 50)
            original_failures += summary['coordinates'] != len(m['coordinates'])
            summaries += 1
        checked += 1
    assert checked == 402 and chain_coordinates == 1608 and summaries == original_failures == 804
    rejected = []
    for name in ['bad_alphabet', 'duplicate_nodes', 'missing_node', 'duplicate_tip', 'zero_length', 'negative_length', 'boolean_length', 'float_length', 'unsorted_tips']:
        bad = deepcopy(first)
        if name == 'bad_alphabet': bad['alphabet'] = 'ACD'
        elif name == 'duplicate_nodes': bad['nodes'][1] = bad['nodes'][0]
        elif name == 'missing_node': bad['nodes'].pop()
        elif name == 'duplicate_tip': bad['tips'][1]['tip'] = bad['tips'][0]['tip']
        elif name == 'zero_length': bad['tips'][0]['length'] = 0
        elif name == 'negative_length': bad['tips'][0]['length'] = -1
        elif name == 'boolean_length': bad['tips'][0]['length'] = True
        elif name == 'float_length': bad['tips'][0]['length'] = 1.5
        elif name == 'unsorted_tips': bad['tips'].reverse()
        try: coordinate_count(bad)
        except AssertionError: rejected.append(name)
        else: raise AssertionError('Accepted invalid metadata: ' + name)
    result = dict(status='passed_full_saved_quartet_coordinate_schema_contract', checked_utc=datetime.now(timezone.utc).isoformat(),
        complete_quartets=checked, chain_coordinate_metadata_checked=chain_coordinates, cutoff_summaries_checked=summaries,
        original_bad_count_reproduced=original_failures, malformed_metadata_rejected=rejected,
        source_hashes=bindings, script_pins={str(p): sha(p) for p in [Path(__file__), 'scripts/close_baliphy_postprocessing_20261002_v2.py']},
        scope='Every 402 saved quartet coordinate schema/count and all 1608 contributing chain coordinate files checked. Full raw array/report provenance hashes and original process journals remain required by the production closer. This contract does not qualify MCMC convergence.')
    with Path('metadata/baliphy_full_postprocessing_coordinate_validation_20261002_v2.json').open('x') as handle:
        handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'script_pins']}, indent=2))


if __name__ == '__main__': main()
