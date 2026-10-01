#!/usr/bin/env python3
"""Synthetic geometry/exclusion fixtures; not a taxon pilot or production proof."""
import argparse
import contextlib
import csv
import gzip
import hashlib
import io
import itertools
import json
import tempfile
from pathlib import Path
import numpy as np
import fit_full_triad_common_residues as producer
import readback_full_triad_common_fits as reader
from full_triad_fit_sources import DEFINITIONS
from run_ortholog_pair_guide_comparison import sha


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True); args = p.parse_args()
    with tempfile.TemporaryDirectory(prefix='full-triad-fit-fixture-') as temp:
        base = Path(temp); root = base / 'maps'; root.mkdir(); out = base / 'fits'
        screens = [dict(id=f'n{n}_c{int(c * 100)}', minimum_aligned_residues=n, minimum_original_coverage=c) for n in [30, 50] for c in [.5, .7, .9]]
        inputs = {}; triads = []; maps = []
        for number in range(3):
            models = [[f'fixture_{number}_{role}', 1] for role in ['a', 'b', 'reference']]
            tid = hashlib.sha256(json.dumps(models, separators=(',', ':')).encode()).hexdigest(); triads.append(dict(triad_id=tid, models=models))
            for column, model in enumerate(models):
                length = 80 if column == 1 else 60; original = np.arange(1, length + 1)
                coords = np.column_stack((original * .8, np.sin(original) * 2, np.cos(original * .7)))
                if number == 2: coords = np.column_stack((original * .8, np.zeros(length), np.zeros(length)))
                if column == 1: coords = coords @ np.array([[0., -1, 0], [1., 0, 0], [0, 0, 1.]]) + [4, 3, 2]
                if column == 2 and number == 1: coords[:, 0] *= -1
                for mask in ['full', 'plddt70']:
                    positions = original.tolist() if mask == 'full' else original[1::2].tolist()
                    path = base / f'{model[0]}-{mask}.pdb'; lines = []; sequence = []
                    for serial, position in enumerate(positions, 1):
                        aa = 'GLY' if column == 1 and position % 4 == 0 else 'ALA'; sequence.append('G' if aa == 'GLY' else 'A')
                        x, y, z = coords[position - 1]; confidence = 85. if mask == 'plddt70' or position % 2 == 0 else 65.
                        lines.append(f'ATOM  {serial:5d}  CA  {aa} A{position:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{confidence:6.2f}          C  \n')
                    path.write_text(''.join(lines)); inputs[(*model, mask)] = dict(path=str(path), sha256=sha(path), model_id=model[0], version=1, mask=mask,
                        sequence=''.join(sequence), original_positions=positions, original_length=length, retained_residues=len(positions), status='ready')
            for mask in ['full', 'plddt70']:
                for orders in itertools.product([0, 1], repeat=3):
                    triples = [[i, i, i] for i in (range(1, 61) if mask == 'full' else range(2, 61, 2))]
                    cycle = triples[:40] if number == 0 and mask == 'full' and orders[0] else triples
                    if number == 0 and mask == 'plddt70' and orders[1]: cycle = triples[:2]
                    problems = [['rmsd_discrepancy'] if number == 1 and mask == 'plddt70' and orders[0] else [], [], []]
                    excls = [{s['id']: (['order1_not_numerically_usable'] if number == 0 and orders[2] and column == 2 else []) for s in screens} for column in range(3)]
                    maps.append(dict(triad_id=tid, role_order=['a', 'b', 'reference'], models=models, mask=mask, orders=list(orders),
                        original_lengths=[60, 80, 60], retained_input_lengths=[inputs[(*model, mask)]['retained_residues'] for model in models],
                        edge_provenance=[dict(mapping_exclusions=x) for x in problems], edge_pair_screen_exclusions=excls,
                        edge_pair_screen_pass=[{sid: not why for sid, why in ex.items()} for ex in excls],
                        all_three_pair_screen_pass={s['id']: not any(e[s['id']] for e in excls) for s in screens},
                        reference_common_triples=triples, cycle_consistent_triples=cycle, common_reference_count=len(triples), cycle_consistent_count=len(cycle),
                        common_core_fit_input_status={d: 'source_excluded' if any(problems) else ('fewer_than_three_common_residues' if len(t) < 3 else 'pending_common_coordinate_geometry') for (d, _), t in zip(DEFINITIONS, [triples, cycle])}))
        with gzip.open(root / 'common_residue_maps.jsonl.gz', 'wt') as f:
            for row in maps: f.write(json.dumps(row) + '\n')
        mapping = dict(target_contexts=3, reference_tie_records=3, duplicate_reference_links=6, correspondence_work_triads=3, mask_order_states=48,
                       common_reference_residue_occurrences=sum(r['common_reference_count'] for r in maps), cycle_consistent_residue_occurrences=sum(r['cycle_consistent_count'] for r in maps))
        (root / 'receipt.json').write_text(json.dumps(mapping)); completion = base / 'synthetic_completion.json'; completion.write_text('{"fixture_only": true}\n')
        plan_path = base / 'plan.json'; plan = dict(output=str(out), mapping_completion=str(completion), resources=dict(minimum_free_disk_gib=0), expected=dict(mask_order_states=48), screens=screens,
                                                 scope='Synthetic fixture only; closed-source proof I/O stubbed, actual PDB parser and full fit/reader algorithms exercised. No production qualification.')
        plan_path.write_text(json.dumps(plan))
        def stub(config, path):
            bindings = {str(path): sha(path), str(root / 'common_residue_maps.jsonl.gz'): sha(root / 'common_residue_maps.jsonl.gz'), str(root / 'receipt.json'): sha(root / 'receipt.json'), str(completion): sha(completion)}
            return triads, inputs, mapping, root, bindings
        producer.load_sources = reader.load_sources = stub
        producer.require_preflight = reader.require_preflight = lambda config, bindings: None
        with contextlib.redirect_stdout(io.StringIO()):
            produced = producer.run(plan_path); checked = reader.run(plan_path, base / 'readback.json')
        assert produced['fit_rows'] == checked['fit_rows'] == 96
        table = out / 'common_residue_fits.tsv.gz'; receipt = out / 'receipt.json'; pristine_receipt = receipt.read_bytes(); pristine_table = table.read_bytes()
        with gzip.open(table, 'rt') as f: records = list(csv.DictReader(f, delimiter='\t'))
        columns = list(records[0]); normal = next(i for i, r in enumerate(records) if r['fit_status'] == 'computed_unique_at_numeric_tolerance' and r['n30_c70_pass'] == '1')
        degenerate = next(i for i, r in enumerate(records) if r['fit_status'] == 'computed_degenerate_geometry')
        short = next(i for i, r in enumerate(records) if r['fit_status'] == 'fewer_than_three_common_residues')
        excluded = next(i for i, r in enumerate(records) if r['fit_status'] == 'source_excluded')
        inherited = next(i for i, r in enumerate(records) if r['n30_c70_core_pass'] == '1' and r['n30_c70_three_pair_pass'] == '0')
        assert float(records[normal]['rmsd_ab']) < .002
        assert any(float(r['rmsd_ar']) > .1 for r in records if r['model_a'] == 'fixture_1_a' and r['fit_status'] == 'computed_unique_at_numeric_tolerance')
        assert all(records[excluded][k] == '' for k in ['rmsd_ab', 'rmsd_ar_minus_br', 'joint_plddt70_fraction'])
        changes = [('changed_rmsd', normal, 'rmsd_ab', '9'), ('changed_signed_contrast', normal, 'rmsd_ar_minus_br', '1'),
                   ('wrong_original_denominator', normal, 'coverage_b', '1'), ('changed_residue_identity', normal, 'sequence_identity_ab', '0'),
                   ('changed_confidence', normal, 'mean_plddt_a', '0'), ('changed_triple_hash', normal, 'triples_sha256', '0' * 64),
                   ('swapped_model_role', normal, 'model_a', records[normal]['model_b']), ('promoted_degenerate', degenerate, 'fit_status', 'computed_unique_at_numeric_tolerance'),
                   ('promoted_short', short, 'fit_status', 'computed_unique_at_numeric_tolerance'), ('cleared_source_exclusion', excluded, 'source_exclusions', ''),
                   ('favorable_inherited_screen', inherited, 'n30_c70_three_pair_pass', '1'), ('removed_state', None, None, None), ('duplicated_state', None, None, None)]
        rejected = []
        for name, index, field, value in changes:
            altered = [dict(r) for r in records]
            if name == 'removed_state': altered.pop()
            elif name == 'duplicated_state': altered.append(dict(altered[-1]))
            else: altered[index][field] = value
            with gzip.open(table, 'wt') as f:
                w = csv.DictWriter(f, fieldnames=columns, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(altered)
            r = json.loads(pristine_receipt); r['artifacts']['common_residue_fits.tsv.gz'] = sha(table); receipt.write_text(json.dumps(r))
            try:
                with contextlib.redirect_stdout(io.StringIO()): reader.run(plan_path, base / (name + '.json'))
            except (AssertionError, ValueError): rejected.append(name)
            else: raise AssertionError('False export accepted: ' + name)
        table.write_bytes(pristine_table); receipt.write_bytes(pristine_receipt)
        x = np.array([[0., 0, 0], [1., 0, 0], [0., 2, 0], [0., 0, 3]])
        y = x.copy(); y[0, 0] += .4; z = x.copy(); z[1, 1] += .7
        first = producer.fit_triplet([x, y, z], ['AAAA'] * 3, [np.full(4, 90.)] * 3)
        swapped = producer.fit_triplet([y, x, z], ['AAAA'] * 3, [np.full(4, 90.)] * 3)
        assert np.isclose(first['rmsd_ar_minus_br'], -swapped['rmsd_ar_minus_br']) and abs(first['rmsd_ar_minus_br']) > .01
        result = dict(status='passed_full_triad_same_residue_fit_software_fixtures', mapping_states=48, fit_rows=96, independent_quaternion_maximum_difference=checked['maximum_absolute_rmsd_or_contrast_difference'],
                      rejected_rehashed_false_exports=rejected, actual_native_pdb_parser=True, proper_rotation_reflection_case=True, sign_reversal_case=True,
                      scope='Synthetic rigid/perturbed/reflected/collinear PDB coordinates, nonconsecutive confidence masks, short2-residue core, inherited pair/source exclusions and complete8order grids. Production closed-source proof I/O is stubbed. No taxon pilot or full-production source qualification.')
        with args.output.open('x') as f: f.write(json.dumps(result, indent=2) + '\n')
        print(json.dumps(result, indent=2))


if __name__ == '__main__': main()
