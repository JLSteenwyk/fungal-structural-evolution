#!/usr/bin/env python3
"""Exact external-sort audit of native directed ortholog pair multiplicities.

Protein integers are zero-based line ordinals in the pinned SequenceIDs.txt.
Records contain two canonical six-hex-digit IDs and one direction digit.
This checks output consistency, not biological orthology or missing pairs.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

from readback_grouped_ortholog_counts import check, digest, require


def encode_pair(left, right):
    require(0 <= left < 2**24 and 0 <= right < 2**24 and left != right,
            'Invalid pair IDs')
    a, b = sorted((left, right))
    return f'{a:06x}{b:06x}{int(left > right)}\n'.encode('ascii')


def summarize_sorted(path, expected, sample_limit=1000):
    totals = dict(directed_incidences=0, unique_unordered_pairs=0,
                  unique_directed_pairs=0, duplicate_directed_incidences=0,
                  pairs_with_repeated_direction=0, pairs_missing_reverse=0,
                  pairs_with_unequal_multiplicity=0, pairs_exactly_once_each_direction=0)
    examples = []
    previous = None
    current = None
    counts = [0, 0]
    h = hashlib.sha256()

    def finish():
        if current is None:
            return
        totals['unique_unordered_pairs'] += 1
        totals['unique_directed_pairs'] += sum(n > 0 for n in counts)
        totals['duplicate_directed_incidences'] += sum(max(0, n - 1) for n in counts)
        totals['pairs_with_repeated_direction'] += int(max(counts) > 1)
        totals['pairs_missing_reverse'] += int(min(counts) == 0)
        totals['pairs_with_unequal_multiplicity'] += int(counts[0] != counts[1])
        totals['pairs_exactly_once_each_direction'] += int(counts == [1, 1])
        if counts != [1, 1] and len(examples) < sample_limit:
            examples.append(dict(left_id=int(current[:6], 16), right_id=int(current[6:], 16),
                                 forward=counts[0], reverse=counts[1]))

    with path.open('rb') as handle:
        for raw in handle:
            h.update(raw)
            require(len(raw) == 14 and raw[-1:] == b'\n' and raw[12:13] in (b'0', b'1')
                    and all(c in b'0123456789abcdef' for c in raw[:12]), 'Malformed encoded pair')
            require(int(raw[:6], 16) < int(raw[6:12], 16), 'Noncanonical encoded pair')
            require(previous is None or previous <= raw, 'Pair stream is not sorted')
            key = raw[:12]
            if key != current:
                finish()
                current = key
                counts = [0, 0]
            counts[raw[12] - 48] += 1
            totals['directed_incidences'] += 1
            previous = raw
    finish()
    require(totals['directed_incidences'] == expected, 'Encoded incidence total differs')
    totals['sorted_pairs_sha256'] = h.hexdigest()
    totals['violation_examples'] = examples
    totals['example_limit'] = sample_limit
    return totals


def protein_ids(snapshot):
    files = {Path(p).name: Path(p) for p in snapshot['input_hashes']}
    require(len(files) == len(snapshot['input_hashes']), 'Ambiguous source files')
    species = {}
    for line in files['SpeciesIDs.txt'].read_text().splitlines():
        native, label = line.split(': ', 1)
        require(native not in species, 'Duplicate native species')
        species[native] = label.rsplit('.', 1)[0]
    require(len(set(species.values())) == len(species), 'Duplicate species label')
    lookup = {label: {} for label in species.values()}
    total = 0
    with files['SequenceIDs.txt'].open() as handle:
        for index, line in enumerate(handle):
            native, protein = line.rstrip('\n').split(': ', 1)
            taxon = species[native.split('_')[0]]
            require(protein not in lookup[taxon] and index < 2**24, 'Duplicate or excessive protein IDs')
            lookup[taxon][protein] = index
            total += 1
    return lookup, total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    for path, expected in plan['pins'].items():
        require(digest(path) == expected, 'Changed plan pin: ' + path)
    identity_plan = Path(plan['identity_plan'])
    require(str(identity_plan) in plan['pins'], 'Identity plan must be pinned')
    validated = check(identity_plan)
    ip = json.loads(identity_plan.read_text())
    identity_receipt = str(Path(ip['output']) / 'receipt.json')
    require(plan['pins'].get(identity_receipt) == digest(identity_receipt), 'Unpinned identity completion')
    snapshot = json.loads(Path(ip['supplement_receipt']).read_text())
    output = Path(plan['output'])
    require(not output.exists(), 'Use a new output directory')
    require(shutil.disk_usage(output.parent).free >= plan['minimum_free_disk_gib'] * 2**30,
            'Insufficient free disk')
    require(plan['pins'].get(plan['sort_executable']) == digest(plan['sort_executable']),
            'Sort executable must be pinned')
    lookup, proteins = protein_ids(snapshot)
    require(proteins == validated['proteins'], 'Protein count differs')
    output.mkdir()
    scratch = output / 'sort-temp'
    scratch.mkdir()
    sorted_path = output / 'sorted_pairs.hex'
    command = [plan['sort_executable'], '--parallel=1', '-S', plan['sort_memory'],
               '-T', str(scratch), '-o', str(sorted_path)]
    csv.field_size_limit(16 * 1024 * 1024)
    total_rows = total_pairs = 0
    with (output / 'sort.stderr').open('wb') as errors:
        process = subprocess.Popen(command, stdin=subprocess.PIPE, stderr=errors,
                                   env=dict(os.environ, LC_ALL='C'))
        try:
            for index, (filename, expected_hash) in enumerate(sorted(snapshot['native_ortholog_hashes'].items()), 1):
                require(shutil.disk_usage(output).free >= plan['emergency_free_disk_gib'] * 2**30,
                        'Emergency disk reserve reached')
                path = Path(filename)
                taxon = path.stem
                h = hashlib.sha256()
                with path.open('rb') as handle:
                    header = handle.readline()
                    h.update(header)
                    require(header.decode().rstrip('\r\n').split('\t') ==
                            ['Orthogroup', 'Species', taxon, 'Orthologs'], 'Unexpected native header')
                    for raw in handle:
                        h.update(raw)
                        row = next(csv.reader([raw.decode()], delimiter='\t'))
                        require(len(row) == 4 and row[1] != taxon, 'Malformed native row')
                        left = [lookup[taxon][x] for x in row[2].split(', ')]
                        right = [lookup[row[1]][x] for x in row[3].split(', ')]
                        for a in left:
                            # Bound temporary allocations for large grouped sides.
                            for start in range(0, len(right), 4096):
                                process.stdin.write(b''.join(encode_pair(a, b) for b in right[start:start+4096]))
                        total_rows += 1
                        total_pairs += len(left) * len(right)
                require(h.hexdigest() == expected_hash, 'Native table differs from audited snapshot')
                state = dict(stage='encoding', tables=index, rows=total_rows, directed_incidences=total_pairs)
                (output / 'state.json').write_text(json.dumps(state) + '\n')
                print(json.dumps(state), flush=True)
            process.stdin.close()
            require(process.wait() == 0, 'External sort failed')
        except BaseException:
            if process.poll() is None:
                process.terminate()
            process.wait()
            raise
    require(total_rows == validated['rows'] and total_pairs == validated['directed_incidences'],
            'Full encoding counts differ from identity audit')
    (output / 'state.json').write_text(json.dumps(dict(stage='summarizing', directed_incidences=total_pairs)) + '\n')
    result = summarize_sorted(sorted_path, total_pairs)
    for path, expected in {**validated['source_sha256'], **plan['pins']}.items():
        require(digest(path) == expected, 'Source changed during execution: ' + path)
    result.update(status='complete_native_ortholog_pair_multiplicity_audit',
                  plan_sha256=digest(args.plan), script_sha256=digest(__file__),
                  identity_receipt_sha256=digest(identity_receipt), rows=total_rows,
                  proteins=proteins, sort_command=command,
                  scope='Exact counts of reciprocal presence and repeated directed pairs across all native grouped rows, using source-line protein IDs. Native snapshot hashes checked during encoding. Does not detect pairs absent in both directions, validate reconciliation semantics or biological orthology, or integrate small-family supplements.')
    (output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'violation_examples'}), flush=True)


if __name__ == '__main__':
    main()
