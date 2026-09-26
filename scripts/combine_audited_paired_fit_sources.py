#!/usr/bin/env python3
"""Combine audited marker tables while preserving original fit provenance.

This creates a table collection, not a native run_paired_marker_fits output.
"""
import argparse
import csv
import json
from pathlib import Path

from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import sha
from prepare_paired_fit_source_inventory import ready_markers

LABELS = {'aa', '3di_af', '3di_af_empirical', '3di_llm'}
TABLES = ['fit_summary.tsv', 'paired_branches.tsv', 'warnings.tsv']


def read_table(path):
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        return reader.fieldnames, list(reader)


def select_rows(rows, selected, source, fits, audit):
    """Retain source values as strings and attach explicit provenance columns."""
    return [dict(row, source_run=source,
                 source_fit_directory=str(fits / row['marker']),
                 source_audit_directory=str(audit))
            for row in rows if row['marker'] in selected]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['inventory', 'expanded', 'previous-audit', 'changed-audit', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new immutable output directory')
    inventory = json.loads(args.inventory.read_text())
    if inventory['status'] != 'passed_paired_fit_source_input_and_settings_compatibility':
        raise ValueError('Source compatibility not established')
    if inventory['script_sha256'] != sha(Path('scripts/prepare_paired_fit_source_inventory.py')):
        raise ValueError('Source inventory producer changed')
    pins = dict(inventory['source_sha256'])
    pins[str(args.inventory)] = sha(args.inventory)
    if pins.get(str(args.expanded / 'receipt.json')) != sha(args.expanded / 'receipt.json'):
        raise ValueError('Different expanded inputs')
    for filename, digest in pins.items():
        if sha(Path(filename)) != digest:
            raise ValueError('Changed source: ' + filename)
    checked_receipt(args.expanded)
    sources = inventory['sources']
    markers = {r['marker'] for r in sources}
    if len(markers) != len(sources) or markers != ready_markers(args.expanded):
        raise ValueError('Inventory does not cover each expanded marker exactly once')
    if {r['source'] for r in sources} != {'unchanged', 'refitted'}:
        raise ValueError('Expected unchanged and refitted source groups')
    tables = {name: [] for name in TABLES}
    fields = {}
    for source, audit_dir in [('unchanged', args.previous_audit), ('refitted', args.changed_audit)]:
        selected_rows = [r for r in sources if r['source'] == source]
        selected = {r['marker'] for r in selected_rows}
        roots = {Path(r['fit_directory']).parent for r in selected_rows}
        if len(roots) != 1:
            raise ValueError('Source group uses multiple fit runs')
        fits = roots.pop()
        # Require completed audits before attempting to consume source fit tables.
        audit = checked_receipt(audit_dir)
        if audit['status'] != 'complete_paired_fit_audit':
            raise ValueError('Source fit audit is incomplete')
        if audit['script_sha256'] != sha(Path('scripts/summarize_paired_marker_fits.py')):
            raise ValueError('Source audit producer changed')
        receipt_path = fits / 'receipt.json'
        receipt = json.loads(receipt_path.read_text())
        if (audit['fit_receipt_sha256'] != sha(receipt_path)
                or receipt['status'] != 'complete_matched_topology_point_estimates'
                or receipt['config_sha256'] != sha(fits / 'config.json')
                or pins.get(str(fits / 'config.json')) != sha(fits / 'config.json')):
            raise ValueError('Source audit/run provenance mismatch')
        pins[str(audit_dir / 'receipt.json')] = sha(audit_dir / 'receipt.json')
        pins[str(receipt_path)] = sha(receipt_path)
        results = {r['marker']: r for r in receipt['results']}
        if len(results) != len(receipt['results']) or not selected <= results.keys():
            raise ValueError('Duplicate or missing source marker results')
        loaded = {}
        for name in TABLES:
            path = audit_dir / name
            if not path.exists() and name == 'warnings.tsv':
                continue
            header, rows = read_table(path)
            if name in fields and fields[name] != header:
                raise ValueError('Source table schemas differ')
            fields[name] = header
            loaded[name] = rows
            tables[name].extend(select_rows(rows, selected, source, fits, audit_dir))
        fit_rows = loaded['fit_summary.tsv']
        keys = {(r['marker'], r['fit']) for r in fit_rows}
        if len(keys) != len(fit_rows) or keys != {(m, label) for m in results for label in LABELS}:
            raise ValueError('Source audit fit grid differs from completed run')
        for row in selected_rows:
            marker = row['marker']
            folder = Path(row['fit_directory'])
            for name, digest in row['expanded_input_sha256'].items():
                if sha(args.expanded / marker / name) != digest or sha(Path(row['input_directory']) / name) != digest:
                    raise ValueError('Selected marker input changed')
            path = folder / 'paired_branches.tsv'
            if sha(path) != results[marker]['paired_branches_sha256']:
                raise ValueError('Source branch artifact changed')
            _, native = read_table(path)
            audited = [r for r in loaded['paired_branches.tsv'] if r['marker'] == marker]
            if native != audited:
                raise ValueError('Audit branch rows differ from native source')
            for label in LABELS:
                fit = checked_receipt_for_label(folder, label)
                config_path = folder / (label + '.config.json')
                config = json.loads(config_path.read_text())
                alignment = 'aa.faa' if label == 'aa' else '3di.faa'
                if (fit['config_sha256'] != sha(config_path)
                        or config['parent_config_sha256'] != sha(fits / 'config.json')
                        or config['alignment_sha256'] != row['expanded_input_sha256'][alignment]):
                    raise ValueError('Selected fit provenance mismatch')
                if label != 'aa' and config['topology_sha256'] != sha(folder / 'aa.treefile'):
                    raise ValueError('Selected fixed topology changed')
    if len(tables['fit_summary.tsv']) != 4 * len(markers):
        raise ValueError('Combined fit grid incomplete')
    branch_keys = [(r['marker'], r['taxon_set_sha256'], r['split_taxa'])
                   for r in tables['paired_branches.tsv']]
    if len(set(branch_keys)) != len(branch_keys):
        raise ValueError('Duplicate combined branch rows')
    args.output.mkdir(parents=True)
    for name, rows in tables.items():
        if not rows:
            continue
        rows.sort(key=lambda r: tuple(r[k] for k in fields[name]))
        with (args.output / name).open('x') as handle:
            writer = csv.DictWriter(handle, fields[name] + ['source_run', 'source_fit_directory', 'source_audit_directory'],
                                    delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)
    result = {'status': 'complete_source_preserving_paired_fit_table_collection',
              'markers': len(markers), 'fits': len(tables['fit_summary.tsv']),
              'paired_branches': len(branch_keys), 'warning_rows': len(tables['warnings.tsv']),
              'source_sha256': pins, 'script_sha256': sha(Path(__file__)),
              'artifacts': {p.name: sha(p) for p in args.output.iterdir()},
              'scope': 'Union of audited point-estimate tables for the expanded inputs, preserving original source runs. Not a native single-run fit output, uncertainty analysis, or test of biological acceleration. Source fits and audits remain required for reproduction.'}
    with (args.output / 'receipt.json').open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: result[k] for k in ['status', 'markers', 'fits', 'paired_branches']}))


def checked_receipt_for_label(folder, label):
    receipt = json.loads((folder / (label + '.receipt.json')).read_text())
    if receipt['status'] != 'completed_point_estimate':
        raise ValueError('Individual fit incomplete')
    for filename, digest in receipt['artifacts'].items():
        if sha(folder / filename) != digest:
            raise ValueError('Individual fit artifact changed')
    return receipt


if __name__ == '__main__':
    main()
