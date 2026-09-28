#!/usr/bin/env python3
"""Report serialized replay discrepancies without conflating precision with fit validation."""
import argparse
import csv
from decimal import Decimal
import json
from pathlib import Path
import re
import subprocess
from ancestral_chain_attempt import sha, write_json


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    assert not args.output.exists()
    pp = Path('metadata/fastml_complete_replay_plan_20260928.json')
    plan = json.loads(pp.read_text())
    for p, h in plan['pins'].items():
        assert sha(p) == h
    state = dict(x.split('=', 1) for x in subprocess.check_output(
        ['systemctl', '--user', 'show', 'fungal-fastml-complete-replay-20260928.service',
         '-p', 'ActiveState', '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
    assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0')
    root = Path(plan['output'])
    receipt = json.loads((root / 'receipt.json').read_text())
    assert receipt['plan_sha256'] == sha(pp)
    for name, h in receipt['artifacts'].items():
        assert sha(root / name) == h
    source = json.loads(Path(plan['source_plan']).read_text())
    rows = []
    pins = {}
    for original in receipt['jobs']:
        assert json.loads((root / (original['job_id'] + '.json')).read_text()) == original
        folder = Path(source['output']) / original['job_id']
        rp = folder / 'receipt.json'
        assert sha(rp) == original['source_receipt_sha256']
        raw = folder / 'RESULTS/EstimatedParameters.txt'
        assert sha(raw) == json.loads(rp.read_text())['artifacts']['RESULTS/EstimatedParameters.txt']
        literal = re.search(r'Log-likelihood=\s*([\deE.+-]+)', raw.read_text()).group(1)
        assert float(literal) == original['reported_log_likelihood']
        half_unit = float(Decimal(5).scaleb(Decimal(literal).as_tuple().exponent - 1))
        row = {k: original[k] for k in ['job_id', 'probability_rows', 'maximum_probability_difference',
               'reported_log_likelihood', 'corrected_log_likelihood', 'likelihood_difference']}
        row.update(reported_likelihood_literal=literal, reported_rounding_half_unit=half_unit,
                   exceeds_reported_rounding_only=abs(row['likelihood_difference']) > half_unit,
                   initial_tree_denominator_difference=original.get('initial_tree_denominator_difference', ''))
        rows.append(row)
        pins[str(raw)] = sha(raw)
    assert len(rows) == 153 and sum(r['probability_rows'] for r in rows) == 8058340
    args.output.mkdir(parents=True)
    with (args.output / 'discrepancies.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    write_json(args.output / 'receipt.json', dict(status='complete_replay_discrepancy_summary',
        replay_receipt_sha256=sha(root / 'receipt.json'), terminal_state=state, source_pins=pins,
        nonempty_inputs=153, empty_inputs=receipt['empty_inputs'], probability_rows=8058340,
        maximum_probability_difference=max(r['maximum_probability_difference'] for r in rows),
        likelihood_absolute_difference_above_001=sum(abs(r['likelihood_difference']) > .01 for r in rows),
        differences_beyond_reported_likelihood_rounding_only=sum(r['exceeds_reported_rounding_only'] for r in rows),
        script_sha256=sha(__file__), artifacts={'discrepancies.tsv': sha(args.output / 'discrepancies.tsv')},
        scope='Output rounding half-unit assumes rounding to the displayed final digit; it does not '
              'bound effects of rounded model parameters or branches. Neither exceeding nor falling '
              'within that interval certifies correctness, optimization or model adequacy.'))


if __name__ == '__main__':
    main()
