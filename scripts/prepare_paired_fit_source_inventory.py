#!/usr/bin/env python3
"""Check input/settings compatibility before combining unchanged and refitted markers.

This inventories source runs; it neither copies fits nor declares them complete.
Both source fit audits must pass before downstream integration.
"""
import argparse
import csv
import json
from pathlib import Path

from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import sha


def ready_markers(root):
    with (root / 'marker_summary.tsv').open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if len({r['marker'] for r in rows}) != len(rows):
        raise ValueError('Duplicate marker summary rows')
    return {r['marker'] for r in rows if r['status'] == 'ready_for_inference'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['previous', 'expanded', 'changed', 'previous-fits',
                 'changed-fits', 'models', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new immutable inventory')
    receipts = {name: checked_receipt(getattr(args, name))
                for name in ['previous', 'expanded', 'changed', 'models']}
    pins = {str(getattr(args, name) / 'receipt.json'):
            sha(getattr(args, name) / 'receipt.json') for name in receipts}
    if any(receipts['previous'][key] != receipts['expanded'][key]
           for key in ['mask', 'eligibility']):
        raise ValueError('Expanded inputs changed filtering')
    for name in ['previous', 'expanded']:
        if receipts['changed'][name + '_receipt_sha256'] != pins[str(getattr(args, name) / 'receipt.json')]:
            raise ValueError('Changed selection has different parent inputs')
    configs = []
    for inputs, fits in [(args.previous, args.previous_fits),
                         (args.changed, args.changed_fits)]:
        path = fits / 'config.json'
        config = json.loads(path.read_text())
        pins[str(path)] = sha(path)
        if config['input_receipt_sha256'] != sha(inputs / 'receipt.json'):
            raise ValueError('Fit run has different inputs')
        if config['model_receipt_sha256'] != sha(args.models / 'receipt.json'):
            raise ValueError('Fit run has different models')
        if sha(Path(config['executable'])) != config['executable_sha256']:
            raise ValueError('Fit executable changed')
        if sha(Path('scripts/run_paired_marker_fits.py')) != config['script_sha256']:
            raise ValueError('Fit producer changed')
        configs.append({k: v for k, v in config.items()
                        if k not in ['input_receipt_sha256', 'input_dimensions']})
    if configs[0] != configs[1]:
        raise ValueError('Source runs have different settings')
    previous, expanded, changed = [ready_markers(getattr(args, name))
                                   for name in ['previous', 'expanded', 'changed']]
    if previous != expanded or not changed <= expanded:
        raise ValueError('Incompatible ready-marker universes')
    rows = []
    observed_changed = set()
    for marker in sorted(expanded):
        names = ['aa.faa', '3di.faa', 'columns.tsv']
        target = {name: sha(args.expanded / marker / name) for name in names}
        unchanged = all(sha(args.previous / marker / name) == target[name]
                        for name in names)
        if not unchanged:
            observed_changed.add(marker)
        source_inputs = args.previous if unchanged else args.changed
        source_fits = args.previous_fits if unchanged else args.changed_fits
        if not unchanged and marker not in changed:
            raise ValueError('Changed marker missing from refit inputs')
        for name in names:
            if sha(source_inputs / marker / name) != target[name]:
                raise ValueError('Source marker differs from expanded input')
        rows.append({'marker': marker,
                     'source': 'unchanged' if unchanged else 'refitted',
                     'input_directory': str(source_inputs / marker),
                     'fit_directory': str(source_fits / marker),
                     'expanded_input_sha256': target})
    if observed_changed != changed:
        raise ValueError('Refit selection is not the exact changed-marker set')
    result = {
        'status': 'passed_paired_fit_source_input_and_settings_compatibility',
        'markers': len(rows), 'unchanged_markers': len(expanded - changed),
        'refitted_markers': len(changed), 'source_sha256': pins,
        'script_sha256': sha(Path(__file__)), 'sources': rows,
        'scope': 'Exact AA, 3Di and column bytes match the expanded cohort for each selected source. Run settings match apart from cohort input hash and dimensions. This is not fit completion or an artifact/topology audit. Require completed, source-bound audits of both runs before combining results; preserve original fit paths/configurations.'}
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: result[k] for k in ['status', 'markers', 'unchanged_markers', 'refitted_markers']}))


if __name__ == '__main__':
    main()
