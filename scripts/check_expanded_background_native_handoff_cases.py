#!/usr/bin/env python3
"""Synthetic four-collection handoff and independent false-export software checks."""
import argparse
import copy
import csv
import gzip
import hashlib
import json
import tempfile
from pathlib import Path
import prepare_expanded_background_native_handoff as producer
import readback_expanded_background_native_handoff as reader
from expanded_background_input_sources import PREFERENCE
from run_ortholog_pair_guide_comparison import sha


def save_gzip(path, records):
    with gzip.open(path, 'wt') as handle:
        for record in records: handle.write(json.dumps(record) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); assert not args.output.exists()
    with tempfile.TemporaryDirectory(prefix='expanded-background-handoff-software-', dir='results') as folder:
        temp = Path(folder); models = {}
        for name in 'abcde':
            raw = temp / (name + '.cif'); raw.write_text('synthetic raw coordinate byte fixture ' + name)
            models[name] = dict(model_id=name, version=1, path=str(raw), sha256=sha(raw), sequence_sha256=hashlib.sha256(b'AAA').hexdigest(), length=3)
        groups = dict(primary='ab', reference='c', background='cd', legacy_reference='be'); specs = {}; all_input_rows = {}
        for label, names in groups.items():
            source = temp / label; source.mkdir(); mp = source / 'models.jsonl'; manifest = source / 'inputs.jsonl'; rows = []; counts = {}
            mp.write_text(''.join(json.dumps(models[n]) + '\n' for n in names))
            for name in names:
                for mask in ['full', 'plddt70']:
                    status = 'too_few_retained_residues' if name == 'd' and mask == 'plddt70' else 'ready'
                    n = 2 if status != 'ready' else 3
                    row = dict(model_id=name, version=1, mask=mask, status=status, source_sha256=models[name]['sha256'], sequence='A' * n,
                               retained_residues=n, original_length=3, original_positions=list(range(1, n + 1)), coordinate_shard=label + '.jsonl.gz')
                    if status == 'ready':
                        pdb = source / (name + '-' + mask + '.pdb'); pdb.write_text('synthetic written PDB byte fixture ' + name + '-' + mask); row.update(path=str(pdb), sha256=sha(pdb))
                    else: row['reason'] = 'fewer_than_three_retained_residues'
                    rows.append(row); counts[mask + ':' + status] = counts.get(mask + ':' + status, 0) + 1
            manifest.write_text(''.join(json.dumps(row) + '\n' for row in rows)); all_input_rows[label] = rows
            specs[label] = dict(model_path=mp, manifest=manifest, receipt=dict(models=len(names), input_dispositions=len(rows), counts=counts))
        pairs, candidates = [], {}
        for a, b, matching in [('a', 'b', ['reference_old']), ('a', 'c', []), ('c', 'd', []), ('a', 'e', [])]:
            ends = [[a, 1], [b, 1]]; key = hashlib.sha256(json.dumps(ends, separators=(',', ':')).encode()).hexdigest()
            pair = dict(pair_key=key, model_a=a, version_a='1', model_b=b, version_b='1', work_disposition='new_model_pair'); pairs.append(pair)
            candidates[key] = dict(**pair, new_sources=json.dumps(['background_expanded']), matching_old_sources=json.dumps(matching), changed_old_sources='[]',
                                   disposition='matching_catalog_sources_pending_input_and_result_checks' if matching else 'new_pair_requires_alignment')
        active = {(m['model_id'], m['version']): m for m in models.values()}
        def sources(plan, path):
            bindings = {str(path): sha(path)}
            for s in specs.values():
                for p in [s['model_path'], s['manifest']]: bindings[str(p)] = sha(p)
            return active, pairs, candidates, specs, bindings
        producer.load_sources = reader.load_sources = sources
        plan = dict(output=str(temp / 'baseline'), resources=dict(minimum_free_disk_gib=0), expected=dict(active_models=5, full_pairs=4, new_pairs=3, pending_catalog_reuse_pairs=1),
                    scope='Synthetic four-collection software fixture; proof I/O stubbed; fake CF/PDB bytes are not geometry fixtures or a biological pilot.')
        pp = temp / 'baseline-plan.json'; pp.write_text(json.dumps(plan)); producer.run(pp); baseline = reader.run(pp, temp / 'baseline-readback.json')
        root = Path(plan['output']); content = {}
        for filename in ['active_models.jsonl.gz', 'active_inputs.jsonl.gz', 'overlapping_input_states.jsonl.gz']:
            with gzip.open(root / filename, 'rt') as handle: content[filename] = [json.loads(line) for line in handle]
        with (root / 'full_background_work_partition.tsv').open() as handle: content['full_background_work_partition.tsv'] = list(csv.DictReader(handle, delimiter='\t'))
        values = {(r['model_id'], r.get('mask')): r for r in content['active_inputs.jsonl.gz']}
        assert values['e', 'full']['collection_sources'][0]['collection'] == 'legacy_reference'
        assert [v['collection'] for v in values['b', 'full']['collection_sources']] == ['primary', 'legacy_reference']
        assert values['d', 'plddt70']['status'] == 'too_few_retained_residues' and 'sha256' not in values['d', 'plddt70']
        assert baseline['directed_native_input_counts'] == {'full:ready': 6, 'plddt70:ready': 4, 'plddt70:input_unavailable': 2}
        assert baseline['overlapping_models'] == 2 and baseline['overlapping_model_mask_states'] == 4
        def input_row(c, name, mask): return next(r for r in c['active_inputs.jsonl.gz'] if r['model_id'] == name and r['mask'] == mask)
        mutations = {
            'wrong_collection_preference': lambda c, r: input_row(c, 'b', 'full')['collection_sources'].reverse(),
            'missing_legacy_active_model': lambda c, r: c['active_models.jsonl.gz'].pop(),
            'missing_confidence_mask': lambda c, r: c['active_inputs.jsonl.gz'].pop(),
            'duplicated_active_input': lambda c, r: c['active_inputs.jsonl.gz'].insert(1, copy.deepcopy(c['active_inputs.jsonl.gz'][0])),
            'promoted_short_input': lambda c, r: input_row(c, 'd', 'plddt70').update(status='ready'),
            'lost_overlap_provenance': lambda c, r: c['overlapping_input_states.jsonl.gz'].pop(),
            'changed_overlap_position': lambda c, r: c['overlapping_input_states.jsonl.gz'][0]['input']['original_positions'].reverse(),
            'wrong_native_partition': lambda c, r: c['full_background_work_partition.tsv'][0].update(measurement_disposition='invented'),
            'changed_physical_endpoint': lambda c, r: c['full_background_work_partition.tsv'][0].update(model_a='wrong'),
            'invented_semantic_signature': lambda c, r: input_row(c, 'a', 'full').update(semantic_sha256='wrong'),
            'changed_disposition_total': lambda c, r: r.update(new_pairs=4),
            'invented_scientific_eligibility': lambda c, r: r.update(scientific_eligibility=True),
        }
        receipt = json.loads((root / 'receipt.json').read_text()); rejected = []
        for name, mutate in mutations.items():
            records = copy.deepcopy(content); r = copy.deepcopy(receipt); mutate(records, r); case = temp / name; case.mkdir()
            config = {**plan, 'output': str(case)}; cp = temp / (name + '-plan.json'); cp.write_text(json.dumps(config))
            for filename, rows in records.items():
                if filename.endswith('.gz'): save_gzip(case / filename, rows)
                else:
                    with (case / filename).open('w') as handle:
                        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t'); writer.writeheader(); writer.writerows(rows)
            r.update(plan_sha256=sha(cp), source_hashes={**r['source_hashes'], str(cp): sha(cp)})
            r['source_hashes'].pop(str(pp)); r['artifacts'] = {name: sha(case / name) for name in records}; (case / 'receipt.json').write_text(json.dumps(r))
            try: reader.run(cp, case / 'readback.json')
            except (AssertionError, StopIteration): rejected.append(name)
            else: raise AssertionError('False export accepted: ' + name)
        # A contradictory overlap in the input source itself cannot be resolved
        # by choosing whichever collection is ready or has favorable geometry.
        altered = copy.deepcopy(all_input_rows['background']); altered[0]['original_positions'].reverse()
        specs['background']['manifest'].write_text(''.join(json.dumps(row) + '\n' for row in altered))
        bad = temp / 'contradictory-source-plan.json'; bad.write_text(json.dumps({**plan, 'output': str(temp / 'contradictory-source')}))
        try: producer.run(bad)
        except AssertionError: rejected.append('producer_rejects_contradictory_source_overlap')
        else: raise AssertionError('Contradictory source overlap accepted')
    result = dict(status='passed_synthetic_expanded_background_full_handoff_checks', baseline_readback_status=baseline['status'], synthetic_active_models=5,
                  synthetic_current_collections=4, synthetic_background_pairs=4, rejected_rehashed_false_exports=rejected[:12], contradictory_source_overlap_rejected=True,
                  script_hashes={str(p): sha(p) for p in [Path(__file__), Path(producer.__file__), Path(reader.__file__), Path('scripts/expanded_background_input_sources.py')]},
                  scope='Synthetic complete four-collection/model/mask/overlap/preference/legacy-fallback/new-vs-pending work fixture. Actual synthetic byte hashes checked, independent SQL baseline and12 rehashed false exports rejected, contradictory full-residue source overlap rejected before native work. Proof I/O stubbed only here; no structural geometry, real source qualification, biological pilot or old-result reuse claim.')
    with args.output.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__': main()
