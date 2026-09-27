"""Inspect six-residue contexts for every jointly observed functional residue."""
import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    source = Path('results/functional_sites/prediction-source-comparison-20260927-v1')
    receipt = json.loads((source / 'receipt.json').read_text())
    assert receipt['status'] == 'complete_full_functional_prediction_source_comparison_with_readback'
    table = source / 'residue_source_comparison.tsv'
    assert sha(table) == receipt['artifacts'][table.name]
    rows = [r for r in csv.DictReader(table.open(), delimiter='\t') if r['paired_coverage'] == 'both']
    assert len(rows) == receipt['residue_coverage_counts']['both'] == 150
    roots = dict(afdb=Path('results/structural_markers/afdb-recovered-union-20260925-v1/encodings'),
                 esmfold=Path('results/structural_alphabet/audited-esmfold-all-completed-20260922-v1'))
    joins = dict(afdb=Path('results/functional_sites/afdb-recovered-linked-20260927-v1'),
                 esmfold=Path('results/functional_sites/esmfold-all-completed-linked-20260923-v1'))
    pins = {str(source / 'receipt.json'): sha(source / 'receipt.json'), str(table): sha(table)}
    indices, arrays = {}, {}
    for label, root in roots.items():
        join = json.loads((joins[label] / 'receipt.json').read_text())
        assert join['source_receipts']['encodings']['sha256'] == sha(root / 'receipt.json')
        assert receipt['source_sha256'][str(joins[label] / 'receipt.json')] == sha(joins[label] / 'receipt.json')
        r = json.loads((root / 'receipt.json').read_text())
        assert r['confidence_stage'] == 'plddt_and_pae'
        path = root / 'model_summary.tsv'
        assert sha(path) == r['artifacts'][path.name]
        entries = list(csv.DictReader(path.open(), delimiter='\t'))
        indices[label] = {(e['model_id'], e['version']): e for e in entries}
        assert len(indices[label]) == len(entries)
        for p in (root / 'receipt.json', path, joins[label] / 'receipt.json'):
            pins[str(p)] = sha(p)
    output = []
    identity = ['marker', 'taxon_id', 'protein_id', 'protein_residue_1based']
    for row in rows:
        result = {k: row[k] for k in identity}
        result.update(paired_state_comparison=row['paired_state_comparison'],
                      any_conserved_candidate=row['any_conserved_candidate'], pfam_accessions=row['pfam_accessions'])
        pos = int(row['protein_residue_1based'])
        sequences = []
        contexts = []
        for label in roots:
            model_key = row[label + '_model_id'], row[label + '_model_version']
            entry = indices[label][model_key]
            path = Path(entry['encoding_path'])
            cache_key = label, model_key
            if cache_key not in arrays:
                assert sha(path) == entry['encoding_sha256']
                pins[str(path)] = sha(path)
                with np.load(path, allow_pickle=False) as values:
                    arrays[cache_key] = {k: values[k].copy() for k in values.files}
            data = arrays[cache_key]
            sequence = str(data['sequence'])
            sequences.append(sequence)
            assert hashlib.sha256(sequence.encode()).hexdigest() == entry['sequence_sha256']
            assert sequence[pos - 1] == row['observed_amino_acid']
            assert data['valid'][pos - 1] and str(data['states'])[pos - 1] == row[label + '_native_state']
            partner = int(data['partner_residue_1based'][pos - 1])
            context = [pos - 1, pos, pos + 1, partner - 1, partner, partner + 1]
            assert all(1 <= value <= len(sequence) for value in context)
            low = float(np.min(data['ca_plddt'][np.array(context) - 1]))
            maximum_pae = float(data['feature_max_pae'][pos - 1])
            assert low >= 70 and maximum_pae <= 10 and math.isfinite(maximum_pae)
            assert math.isclose(low, float(data['feature_min_plddt'][pos - 1]), abs_tol=1e-10)
            assert math.isclose(low, float(row[label + '_feature_min_plddt']), abs_tol=1e-10)
            assert math.isclose(maximum_pae, float(row[label + '_feature_max_pae']), abs_tol=1e-10)
            contexts.append(context)
            result.update({label + '_model_id': model_key[0], label + '_model_version': model_key[1],
                label + '_native_state': str(data['states'])[pos - 1],
                label + '_partner_residue_1based': partner,
                label + '_context_residues_1based': ','.join(map(str, context)),
                label + '_context_amino_acids': ''.join(sequence[i - 1] for i in context),
                label + '_feature_min_plddt': low, label + '_feature_max_pae': maximum_pae})
        assert sequences[0] == sequences[1]
        result.update(sequence_sha256=hashlib.sha256(sequences[0].encode()).hexdigest(),
                      partner_comparison='same' if contexts[0][4] == contexts[1][4] else 'different',
                      distinct_context_overlap=len(set(contexts[0]) & set(contexts[1])))
        output.append(result)
    args.output.mkdir(parents=True)
    target = args.output / 'functional_context_comparison.tsv'
    with target.open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(output)
    # Reload every exported context and check its positions/confidence against
    # the qualified arrays; this checks serialization, not fresh native encoding.
    readback = list(csv.DictReader(target.open(), delimiter='\t'))
    assert {tuple(r[k] for k in identity) for r in readback} == {tuple(r[k] for k in identity) for r in rows}
    assert len(readback) == len(rows)
    for row in readback:
        contexts = []
        for label in roots:
            data = arrays[label, (row[label + '_model_id'], row[label + '_model_version'])]
            focal = int(row['protein_residue_1based']) - 1
            partner = int(data['partner_residue_1based'][focal]) - 1
            context = [focal + offset for offset in (-1, 0, 1)] + [partner + offset for offset in (-1, 0, 1)]
            assert list(map(int, row[label + '_context_residues_1based'].split(','))) == [i + 1 for i in context]
            assert int(row[label + '_partner_residue_1based']) == partner + 1
            assert row[label + '_context_amino_acids'] == ''.join(str(data['sequence'])[i] for i in context)
            assert float(row[label + '_feature_min_plddt']) == min(float(data['ca_plddt'][i]) for i in context)
            assert float(row[label + '_feature_max_pae']) == float(data['feature_max_pae'][focal])
            contexts.append(context)
        assert int(row['distinct_context_overlap']) == len(set(contexts[0]).intersection(contexts[1]))
        assert (row['partner_comparison'] == 'same') == (contexts[0][4] == contexts[1][4])
        assert (row['paired_state_comparison'] == 'same') == (row['afdb_native_state'] == row['esmfold_native_state'])
    for path, expected in pins.items():
        assert sha(path) == expected
    counts = Counter(r['paired_state_comparison'] + '_state:' + r['partner_comparison'] + '_partner' for r in readback)
    result = dict(status='complete_full_joint_functional_context_comparison', residues=len(output),
        qualified_models_loaded=len(arrays), state_partner_counts=dict(counts),
        source_sha256=pins, script_sha256=sha(__file__), artifacts={target.name: sha(target)},
        scope='All jointly observed functional coordinates, including agreeing states. Exact sequences, partners, six-residue context and pLDDT minima checked against qualified arrays; every exported context reloaded. Directional PAE maxima read from audited arrays, not recomputed here. Same partner does not imply same geometry; changed partner does not prove causation of a state difference. No fresh coordinate reconstruction, native encoding, predictor accuracy or evolutionary inference.')
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('source_sha256', 'artifacts')}, indent=2))


if __name__ == '__main__':
    main()
