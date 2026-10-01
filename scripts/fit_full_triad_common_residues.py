#!/usr/bin/env python3
"""Fit full duplicate/reference triples on the identical original residue sets."""
import argparse
import csv
import gzip
import json
import shutil
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from duplication_alignment_numeric_readback import load_pdb
from fit_whole_protein_common_residues import fit_triplet
from full_triad_fit_sources import load_sources, iterate_maps, DEFINITIONS, METRICS, fields, triple_sha
from full_triad_fit_preflight_gate import require_preflight
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def run(plan_path):
    plan = json.loads(Path(plan_path).read_text()); ph = sha(plan_path)
    assert shutil.disk_usage('.').free >= plan['resources']['minimum_free_disk_gib'] * 2 ** 30
    triads, inputs, mapping, root, bindings = load_sources(plan, plan_path); fitted_inputs = set()
    require_preflight(plan, bindings)
    @lru_cache(maxsize=1024)
    def pdb(key):
        row = inputs[key]; assert row['status'] == 'ready'; bind(bindings, row['path'], row['sha256']); fitted_inputs.add(key)
        seq, xyz, confidence = load_pdb(row)
        return seq, xyz, confidence, {p: i for i, p in enumerate(row['original_positions'])}
    @lru_cache(maxsize=128)
    def calculate(models, mask, triples):
        coords = []; letters = []; confidence = []
        for column, model in enumerate(models):
            seq, xyz, conf, index = pdb((*model, mask)); ix = [index[t[column]] for t in triples]
            coords.append(xyz[ix]); letters.append(''.join(seq[i] for i in ix)); confidence.append(conf[ix])
        return fit_triplet(coords, letters, confidence)
    out = Path(plan['output']); out.mkdir(parents=True, exist_ok=False)
    counts = Counter(); passes = Counter(); core_passes = Counter(); triple_counts = Counter(); states = rows = 0
    with gzip.open(out / 'common_residue_fits.tsv.gz', 'xt', compresslevel=1) as f:
        writer = csv.DictWriter(f, fieldnames=fields(plan), delimiter='\t', lineterminator='\n'); writer.writeheader()
        for item in iterate_maps(root, triads, inputs):
            models = tuple(tuple(m) for m in item['models']); mask = item['mask']; lengths = item['original_lengths']; retained = item['retained_input_lengths']
            problems = sorted({p for edge in item['edge_provenance'] for p in edge['mapping_exclusions']})
            for definition, field in DEFINITIONS:
                triples = tuple(tuple(t) for t in item[field]); n = len(triples)
                row = dict(triad_id=item['triad_id'], mask=mask, order_ab=item['orders'][0], order_ar=item['orders'][1], order_br=item['orders'][2],
                           mapping_definition=definition, triples_sha256=triple_sha(triples), common_residues=n,
                           source_exclusions=';'.join(problems), mapping_disagreement_count=item['common_reference_count'] - item['cycle_consistent_count'])
                for i, label in enumerate(['a', 'b', 'reference']):
                    row.update({f'model_{label}': models[i][0], f'version_{label}': models[i][1], f'length_{label}': lengths[i],
                                f'retained_{label}': retained[i], f'coverage_{label}': n / lengths[i],
                                f'retained_coverage_{label}': n / retained[i] if retained[i] else ''})
                row.update({k: '' for k in METRICS})
                if problems: row['fit_status'] = 'source_excluded'
                elif n < 3: row['fit_status'] = 'fewer_than_three_common_residues'
                else: row.update(calculate(models, mask, triples))
                for screen in plan['screens']:
                    sid = screen['id']; why = list(problems); cutoff = Fraction(str(screen['minimum_original_coverage']))
                    if n < screen['minimum_aligned_residues']: why.append('short_common_core')
                    if any(Fraction(n, length) < cutoff for length in lengths): why.append('low_original_protein_coverage')
                    if row['fit_status'] == 'fewer_than_three_common_residues': why.append('fewer_than_three_common_residues')
                    if row['fit_status'] == 'computed_degenerate_geometry': why.append('degenerate_common_core_geometry')
                    inherited = [label + ':' + reason for label, reasons in zip(['ab', 'ar', 'br'], item['edge_pair_screen_exclusions']) for reason in reasons[sid]]
                    assert item['all_three_pair_screen_pass'][sid] == (not inherited)
                    row[sid + '_core_pass'] = int(not why); row[sid + '_core_exclusions'] = ';'.join(why)
                    row[sid + '_three_pair_pass'] = int(not inherited); row[sid + '_three_pair_exclusions'] = ';'.join(inherited)
                    row[sid + '_pass'] = int(not why and not inherited); row[sid + '_exclusions'] = ';'.join(why + inherited)
                    if row[sid + '_pass']: assert row['fit_status'] == 'computed_unique_at_numeric_tolerance'
                    key = mask + ':' + definition + ':' + sid
                    passes[key] += row[sid + '_pass']; core_passes[key] += row[sid + '_core_pass']
                writer.writerow(row); rows += 1; counts[mask + ':' + definition + ':' + row['fit_status']] += 1; triple_counts[mask + ':' + definition] += n
            states += 1
            if states % 10000 == 0: print('Full same-residue fit states', states, '/', mapping['mask_order_states'], flush=True)
    assert states == mapping['mask_order_states'] == plan['expected']['mask_order_states'] and rows == 2 * states
    assert sum(triple_counts.values()) == mapping['common_reference_residue_occurrences'] + mapping['cycle_consistent_residue_occurrences']
    verify(bindings)
    result = dict(status='complete_full_triad_same_residue_fits_pending_independent_readback', plan_sha256=ph,
                  mapping_receipt_sha256=sha(root / 'receipt.json'), mapping_completion_sha256=sha(plan['mapping_completion']),
                  target_contexts=mapping['target_contexts'], reference_tie_records=mapping['reference_tie_records'],
                  duplicate_reference_links=mapping['duplicate_reference_links'], correspondence_work_triads=len(triads), mapping_states=states,
                  fit_rows=rows, counts=dict(counts), core_screen_pass_counts=dict(core_passes), screen_pass_counts=dict(passes),
                  triple_occurrences=dict(triple_counts), fitted_pdb_inputs=len(fitted_inputs),
                  source_hashes=bindings, artifacts={'common_residue_fits.tsv.gz': sha(out / 'common_residue_fits.tsv.gz')},
                  scientific_eligibility=False, scope=plan['scope'])
    with (out / 'receipt.json').open('x') as f: f.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes']}, indent=2), flush=True)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True)
    run(p.parse_args().plan)


if __name__ == '__main__': main()
