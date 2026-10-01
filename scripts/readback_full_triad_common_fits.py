#!/usr/bin/env python3
"""Reconstruct every full same-residue fit through independent quaternion rotations."""
import argparse
import csv
import gzip
import json
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import numpy as np
from duplication_alignment_numeric_readback import load_pdb
from readback_whole_protein_common_fits import core_metrics, compare_row
from full_triad_fit_sources import load_sources, iterate_maps, DEFINITIONS, METRICS, fields, triple_sha, SUMMARY_FIELDS
from full_triad_fit_preflight_gate import require_preflight
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def run(plan_path, output):
    plan = json.loads(Path(plan_path).read_text()); ph = sha(plan_path)
    triads, inputs, mapping, root, bindings = load_sources(plan, plan_path); fitted_inputs = set()
    require_preflight(plan, bindings)
    out = Path(plan['output']); rp = out / 'receipt.json'; r = json.loads(rp.read_text()); rh = sha(rp)
    assert r['status'] == 'complete_full_triad_same_residue_fits_pending_independent_readback' and r['plan_sha256'] == ph
    assert r['mapping_receipt_sha256'] == sha(root / 'receipt.json') and r['mapping_completion_sha256'] == sha(plan['mapping_completion'])
    assert r['artifacts'] == {'common_residue_fits.tsv.gz': sha(out / 'common_residue_fits.tsv.gz')}
    @lru_cache(maxsize=1024)
    def residues(key):
        row = inputs[key]; assert row['status'] == 'ready'; bind(bindings, row['path'], row['sha256']); fitted_inputs.add(key)
        letters, coords, confidence = load_pdb(row)
        return {position: (coords[i], letters[i], confidence[i]) for i, position in enumerate(row['original_positions'])}
    @lru_cache(maxsize=128)
    def recompute(models, mask, triples):
        values = [residues((*model, mask)) for model in models]
        parts = [[values[column][triple[column]] for triple in triples] for column in [0, 1, 2]]
        return core_metrics([np.array([r[0] for r in part]) for part in parts],
                            [''.join(r[1] for r in part) for part in parts], [[r[2] for r in part] for part in parts])
    counts = Counter(); passes = Counter(); core_passes = Counter(); triple_counts = Counter(); rows = states = 0; maximum = 0.
    with gzip.open(out / 'common_residue_fits.tsv.gz', 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t'); assert reader.fieldnames == fields(plan)
        for item in iterate_maps(root, triads, inputs):
            ids = tuple(tuple(v) for v in item['models']); mask = item['mask']
            lengths = [inputs[(*v, 'full')]['original_length'] for v in ids]
            retained = [len(inputs[(*v, mask)]['original_positions']) for v in ids]
            problems = sorted(set(p for edge in item['edge_provenance'] for p in edge['mapping_exclusions']))
            for definition, field in DEFINITIONS:
                triples = tuple(tuple(v) for v in item[field]); n = len(triples)
                wanted = dict(triad_id=item['triad_id'], mask=mask, order_ab=item['orders'][0], order_ar=item['orders'][1], order_br=item['orders'][2],
                              mapping_definition=definition, triples_sha256=triple_sha(triples), common_residues=n,
                              source_exclusions=';'.join(problems), mapping_disagreement_count=len(item['reference_common_triples']) - len(item['cycle_consistent_triples']))
                for column, role in enumerate(['a', 'b', 'reference']):
                    wanted.update({f'model_{role}': ids[column][0], f'version_{role}': ids[column][1], f'length_{role}': lengths[column],
                                   f'retained_{role}': retained[column], f'coverage_{role}': n / lengths[column],
                                   f'retained_coverage_{role}': n / retained[column] if retained[column] else ''})
                wanted.update({k: '' for k in METRICS})
                if problems: wanted['fit_status'] = 'source_excluded'
                elif n < 3: wanted['fit_status'] = 'fewer_than_three_common_residues'
                else: wanted.update(recompute(ids, mask, triples))
                for screen in plan['screens']:
                    sid = screen['id']; reasons = list(problems); cutoff = Fraction(str(screen['minimum_original_coverage']))
                    if n < screen['minimum_aligned_residues']: reasons.append('short_common_core')
                    if any(n * cutoff.denominator < length * cutoff.numerator for length in lengths): reasons.append('low_original_protein_coverage')
                    if wanted['fit_status'] == 'fewer_than_three_common_residues': reasons.append('fewer_than_three_common_residues')
                    if wanted['fit_status'] == 'computed_degenerate_geometry': reasons.append('degenerate_common_core_geometry')
                    inherited = []
                    for column, label in enumerate(['ab', 'ar', 'br']):
                        inherited.extend(label + ':' + reason for reason in item['edge_pair_screen_exclusions'][column][sid])
                        assert item['edge_pair_screen_pass'][column][sid] == (not item['edge_pair_screen_exclusions'][column][sid])
                    assert item['all_three_pair_screen_pass'][sid] == (not inherited)
                    core_okay = not reasons; three_okay = not inherited
                    wanted.update({sid + '_core_pass': int(core_okay), sid + '_core_exclusions': ';'.join(reasons),
                                   sid + '_three_pair_pass': int(three_okay), sid + '_three_pair_exclusions': ';'.join(inherited),
                                   sid + '_pass': int(core_okay and three_okay), sid + '_exclusions': ';'.join(reasons + inherited)})
                    key = mask + ':' + definition + ':' + sid
                    core_passes[key] += core_okay; passes[key] += core_okay and three_okay
                actual = next(reader, None); assert actual is not None, 'Missing fit row'
                maximum = max(maximum, compare_row(actual, wanted)); rows += 1
                counts[mask + ':' + definition + ':' + wanted['fit_status']] += 1; triple_counts[mask + ':' + definition] += n
            states += 1
            if states % 10000 == 0: print('Independent full same-residue fit states', states, '/', mapping['mask_order_states'], flush=True)
        assert next(reader, None) is None, 'Extra fit row'
    summary = dict(target_contexts=mapping['target_contexts'], reference_tie_records=mapping['reference_tie_records'],
                   duplicate_reference_links=mapping['duplicate_reference_links'], correspondence_work_triads=len(triads), mapping_states=states,
                   fit_rows=rows, counts=dict(counts), core_screen_pass_counts=dict(core_passes), screen_pass_counts=dict(passes),
                   triple_occurrences=dict(triple_counts), fitted_pdb_inputs=len(fitted_inputs))
    assert states == mapping['mask_order_states'] == plan['expected']['mask_order_states'] and rows == 2 * states
    assert all(r[k] == summary[k] for k in SUMMARY_FIELDS) and r['source_hashes'] == bindings
    verify(bindings); assert sha(rp) == rh and sha(out / 'common_residue_fits.tsv.gz') == r['artifacts']['common_residue_fits.tsv.gz']
    result = dict(status='passed_full_triad_same_residue_quaternion_readback', plan_sha256=ph, producer_receipt_sha256=rh,
                  **summary, maximum_absolute_rmsd_or_contrast_difference=maximum, source_hashes=bindings,
                  scientific_eligibility=False, scope='Full ordered grid and both cores rebuilt from closed mappings and hashed PDB positions. Independent quaternion eigenproblem checks proper rotations, distances, signed contrasts, geometry, confidence/sequence identities, blank exclusions and all rational original-length core/inherited three-edge coverage decisions. RMSD metric bound checked. No biological asymmetry, directional rate, prediction-error calibration or accepted phylogeny inferred.')
    with Path(output).open('x') as f: f.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); run(args.plan, args.output)


if __name__ == '__main__': main()
