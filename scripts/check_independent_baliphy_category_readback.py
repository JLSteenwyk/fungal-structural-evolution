#!/usr/bin/env python3
"""Complete 405-group categorical workflow contracts with genuine oracle metrics."""
import argparse
from collections import Counter
import copy
import gzip
import json
import os
from pathlib import Path
import tempfile
import warnings

import numpy as np

from ancestral_categorical_diagnostics import diagnose_states as oracle
from ancestral_chain_diagnostics import diagnose as scalar_oracle
from ancestral_state_patterns import group_patterns as original_patterns
from check_independent_baliphy_scalar_readback_v2 import setup as scalar_fixture
from independent_baliphy_category_sources import ALPHABET
from prepare_independent_baliphy_category_readback import run
from prepare_independent_baliphy_scalar_readback_v2 import atomic
from readback_independent_baliphy_category_readback import run as readback
from run_ortholog_pair_guide_comparison import sha


def setup(root):
    path = scalar_fixture(root); plan = json.loads(path.read_text())
    source = root / 'source'; complete_path = Path(plan['diagnostic_completion'])
    complete = json.loads(complete_path.read_text()); archive_path = Path(complete['full_hash_archive'])
    archive = json.loads(archive_path.read_text()); bindings = archive['source_hashes']
    def create(p, value): atomic(p, value); bindings[str(p)] = sha(p); return p
    rng = np.random.default_rng(3334); bank = np.zeros((4, 4, 101), dtype=np.uint8)
    bank[1] = rng.integers(0, 2, (4, 101), dtype=np.uint8)
    bank[2] = rng.integers(20, 22, (4, 101), dtype=np.uint8)
    bank[3, :, 63:] = 1
    mapping = np.tile(np.arange(4), 4)
    values = bank[mapping].transpose(1, 2, 0).reshape(4, 101, 4, 4)
    iterations = np.arange(0, 1001, 10)
    free = np.stack([((np.arange(101)[:, None] + i) % np.arange(2, 6)).astype(np.int32) for i in range(4)])
    coords = dict(nodes=['node-' + str(i) for i in range(4)], tips=[dict(tip='tip', length=4)],
        alphabet=''.join(ALPHABET), input_alignment_sha256='a' * 64)
    templates = source / 'categorical-templates'; templates.mkdir()
    quartet = templates / 'quartet.npz'
    np.savez_compressed(quartet, values=values, iterations=iterations, unanchored_residue_counts=free)
    bindings[str(quartet)] = sha(quartet)
    state_templates = []
    for i in range(4):
        p = templates / ('states-' + str(i) + '.npz')
        counts = {f'counts_after_{cut}':np.stack([np.count_nonzero(values[i, iterations > cut] == s, axis=0)
            for s in range(len(ALPHABET))], axis=-1).astype(np.uint16) for cut in [250, 500]}
        np.savez_compressed(p, states=values[i], iterations=iterations, unanchored_residue_counts=free[i], **counts)
        bindings[str(p)] = sha(p); state_templates.append(p)
    cutoff_templates = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        for cutoff in ['250', '500']:
            retained = values[:, iterations > int(cutoff)]
            patterns, ids = original_patterns(retained); multiplicity = np.bincount(ids.ravel())
            p = templates / ('patterns-' + cutoff + '.npz')
            np.savez_compressed(p, patterns=patterns, coordinate_pattern_ids=ids); bindings[str(p)] = sha(p)
            report = templates / ('diagnostics-' + cutoff + '.jsonl.gz'); statuses = Counter(); cs = Counter()
            with gzip.open(report, 'wt') as f:
                for number, (pattern, mult) in enumerate(zip(patterns, multiplicity)):
                    diagnostic = oracle(pattern, ALPHABET)
                    statuses[diagnostic['status']] += 1; cs[diagnostic['status']] += int(mult)
                    f.write(json.dumps(dict(pattern_id=number, coordinate_multiplicity=int(mult), diagnostic=diagnostic), allow_nan=False) + '\n')
            bindings[str(report)] = sha(report)
            summary = dict(discard_through=int(cutoff), retained_samples_per_chain=retained.shape[1],
                patterns=len(patterns), coordinates=16, pattern_status_counts=dict(statuses), coordinate_status_counts=dict(cs),
                unanchored_count_screens={node:scalar_oracle(free[:, iterations > int(cutoff), i]) for i,node in enumerate(coords['nodes'])})
            cutoff_templates[cutoff] = dict(patterns=p, report=report, summary=summary)
    producer_path = Path(complete['producer_receipt']); producer = json.loads(producer_path.read_text())
    recovery = json.loads(Path(plan['recovery_completion']).read_text())
    overlay_path = Path(recovery['overlay']); overlay = json.loads(overlay_path.read_text()); states = {}
    rows = {r['chain']['chain_id']:r for r in overlay['rows']}
    pattern_counts = Counter(); coordinate_counts = Counter(); statuses = Counter(); cs = Counter()
    for group, info in producer['groups'].items():
        folder = source / group; chain_manifests = []
        for i, cid in enumerate(info['chain_ids']):
            native = rows[cid]; native['chain']['alignment_sha256'] = coords['input_alignment_sha256']
            selected = native['selected_disposition']; audit_path = Path(selected['sample_audit'])
            audit = json.loads(audit_path.read_text()); audit['candidate_samples'] = [dict(source_node=node) for node in coords['nodes']]
            create(audit_path, audit); selected['sample_audit_sha256'] = sha(audit_path)
            sf = folder / (cid + '-states'); sf.mkdir(); ap = sf / 'states.npz'; os.link(state_templates[i], ap)
            bindings[str(ap)] = sha(ap); cp = create(sf / 'coordinates.json', coords)
            rp = create(sf / 'receipt.json', dict(summaries=[dict(candidate_anchor_coordinates=16,
                source_audit_sha256=sha(audit_path), artifacts={str(ap):sha(ap), str(cp):sha(cp)})]))
            states[cid] = dict(receipt=str(rp), receipt_sha256=sha(rp))
            chain_manifests.append(dict(chain_id=cid, seed=native['chain']['seed'], model_input_identity=group,
                log=audit['scalar_log'], log_sha256=audit['scalar_log_sha256']))
        if info['status'] == 'unresolved_failed_native_chain_retained': continue
        target = folder / 'categorical'; target.mkdir(); ap = target / 'quartet.npz'; os.link(quartet, ap); bindings[str(ap)] = sha(ap)
        mp = create(target / 'manifest.json', dict(status='provenance_checked_state_quartet', chains=chain_manifests,
            coordinates=coords, expected_iterations=iterations.tolist(), arrays=str(ap), arrays_sha256=sha(ap)))
        report_root = target / 'reports'; report_root.mkdir(); outputs = {}; artifacts = {}
        for cutoff, template in cutoff_templates.items():
            cut = report_root / ('discard-' + cutoff); cut.mkdir()
            for name, original in [('patterns.npz', template['patterns']), ('diagnostics.jsonl.gz', template['report'])]:
                dest = cut / name; os.link(original, dest); bindings[str(dest)] = sha(dest); artifacts[str(dest.relative_to(report_root))] = sha(dest)
            summary = template['summary']; sp = create(cut / 'summary.json', summary); artifacts[str(sp.relative_to(report_root))] = sha(sp)
            outputs[cutoff] = {k:summary[k] for k in ['patterns', 'coordinates', 'retained_samples_per_chain']}
            pattern_counts[cutoff] += summary['patterns']; coordinate_counts[cutoff] += summary['coordinates']
            statuses.update({cutoff + ':' + k:v for k,v in summary['pattern_status_counts'].items()})
            cs.update({cutoff + ':' + k:v for k,v in summary['coordinate_status_counts'].items()})
        rp = create(report_root / 'receipt.json', dict(status='both_cutoff_categorical_reports_complete_not_posterior_qualification',
            arviz_version='0.22.0', numpy_version='2.2.6', outputs=outputs, artifacts=artifacts,
            pins={str(mp):sha(mp), str(ap):sha(ap)}))
        cp = create(target / 'receipt.json', dict(status='verified_quartet_categorical_reports_complete_not_posterior_qualification',
            group=group, manifest=str(mp), manifest_sha256=sha(mp), report=str(rp), report_sha256=sha(rp)))
        info['categorical'] = dict(receipt=str(cp), receipt_sha256=sha(cp))
    create(overlay_path, overlay); recovery['overlay_sha256'] = sha(overlay_path); create(Path(plan['recovery_completion']), recovery)
    producer['states'] = states; create(producer_path, producer)
    for k,v in [('categorical_pattern_counts',dict(pattern_counts)), ('categorical_coordinate_counts',dict(coordinate_counts)),
        ('pattern_status_counts',dict(statuses)), ('coordinate_status_counts',dict(cs))]:
        complete[k] = v; archive['summary'][k] = v
    archive['source_hashes'] = bindings; atomic(archive_path, archive)
    complete.update(full_hash_archive_sha256=sha(archive_path), bound_source_hashes=len(bindings), producer_receipt_sha256=sha(producer_path))
    atomic(complete_path, complete); plan['scope'] = 'Full synthetic405group categorical workflow; genuine locked-oracle metrics on shared four-pattern fixtures. Not production posterior or native-journal evidence.'
    atomic(path, plan); return path


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); rejected = []
    with tempfile.TemporaryDirectory(prefix='independent-baliphy-category-contract-') as directory:
        root = Path(directory); path = setup(root)
        try: run(path, stop_after_groups=3)
        except InterruptedError: pass
        else: raise AssertionError('Interruption was not exercised')
        before = {str(p):sha(p) for p in (root / 'output/groups').rglob('*') if p.is_file()}
        assert len(list((root / 'output/groups').glob('*.json'))) == 3
        result = run(path)
        assert result['pattern_cutoff_rows'] == 3224 and result['declared_indicator_rows'] == 70928
        assert result['singular_indicator_rows'] == 806 and result['unanchored_count_rows'] == 3224
        assert all(sha(p) == d for p,d in before.items())
        readback(path, root / 'output/readback.json')
        for action in ['completed_producer', 'alternate_reader']:
            try: run(path) if action == 'completed_producer' else readback(path, root / 'alternate.json')
            except AssertionError: pass
            else: raise AssertionError('Completed restart accepted')
        (root / 'output/readback.json').unlink()
        manifest = root / 'output/group_manifest.json'; receipt = root / 'output/receipt.json'
        manifest_bytes, receipt_bytes = manifest.read_bytes(), receipt.read_bytes()
        first = json.loads(manifest_bytes)[0]; cp = root / 'output' / first['path']; cp_bytes = cp.read_bytes()
        record = json.loads(cp_bytes); cut = record['result']['cutoffs']['250']
        exported = root / 'output' / cut['serialized_comparison']; original_gzip = exported.read_bytes()
        with gzip.open(exported, 'rt') as f: lines = [json.loads(line) for line in f]
        changes = ['missing_pattern', 'duplicate_pattern', 'changed_metric', 'missing_state', 'changed_unknown_count',
            'source_digest', 'false_science', 'false_failed_disposition', 'changed_multiplicity', 'hidden_cutoff', 'changed_unanchored', 'omit_failed_group']
        for change in changes:
            altered = copy.deepcopy(lines); checkpoint = json.loads(cp_bytes)
            if change == 'missing_pattern': altered.pop()
            elif change == 'duplicate_pattern': altered[-1] = copy.deepcopy(altered[0])
            elif change == 'changed_metric': altered[1]['independent']['indicators']['A']['screen']['bulk_ess'] += 1
            elif change == 'missing_state': del altered[0]['independent']['indicators']['X']
            elif change == 'changed_unknown_count': altered[2]['independent']['indicators']['X']['counts_per_chain'][0] += 1
            elif change == 'source_digest': altered[0]['source_diagnostic_sha256'] = '0' * 64
            elif change == 'false_science': altered[0]['scientific_eligibility'] = True
            elif change == 'false_failed_disposition': checkpoint['result']['status'] = 'unresolved_failed_native_chain_retained'; checkpoint['result']['cutoffs'] = {}
            elif change == 'changed_multiplicity': altered[0]['coordinate_multiplicity'] += 1
            elif change == 'hidden_cutoff': del checkpoint['result']['cutoffs']['500']
            elif change == 'changed_unanchored': checkpoint['result']['cutoffs']['250']['unanchored_count_rows']['node-0']['independent']['bulk_ess'] += 1
            with gzip.open(exported, 'wt') as f:
                for row in altered: f.write(json.dumps(row, allow_nan=False) + '\n')
            if '250' in checkpoint['result']['cutoffs']:
                checkpoint['result']['cutoffs']['250']['serialized_comparison_sha256'] = sha(exported)
            atomic(cp, checkpoint)
            entries = json.loads(manifest_bytes); entries[0]['sha256'] = sha(cp); atomic(manifest, entries)
            if change == 'omit_failed_group':
                assert entries[-1]['group'] == 'group-0404'
                entries.pop(); atomic(manifest, entries)
            document = json.loads(receipt_bytes)
            for p in [cp, exported, manifest]: document['artifacts'][str(p.relative_to(root / 'output'))] = sha(p)
            atomic(receipt, document)
            try: readback(path, root / 'malformed.json')
            except AssertionError: rejected.append(change)
            else: raise AssertionError(('Rehashed false export accepted', change))
            exported.write_bytes(original_gzip); cp.write_bytes(cp_bytes); manifest.write_bytes(manifest_bytes); receipt.write_bytes(receipt_bytes)
    paths = [__file__, *['scripts/' + name + '.py' for name in ['independent_baliphy_category_sources',
        'independent_baliphy_category_replay', 'prepare_independent_baliphy_category_readback',
        'readback_independent_baliphy_category_readback', 'independent_ancestral_categorical_diagnostics',
        'independent_ancestral_scalar_diagnostics_v2', 'check_independent_baliphy_scalar_readback_v2']]]
    result = dict(status='passed_full_independent_baliphy_categorical_grid_and_checkpoint_contracts',
        full_quartets=405, complete_quartets=403, unresolved_quartets=2,
        pattern_cutoff_rows=3224, declared_indicator_rows=70928, unanchored_count_rows=3224,
        singular_indicator_rows_retained=806, full_alphabet=ALPHABET,
        unchanged_interrupted_checkpoints_reused=True, completed_producer_and_alternate_reader_restart_refused=True,
        rehashed_false_exports_rejected=rejected, source_hashes={str(p):sha(p) for p in paths},
        scientific_eligibility=False, production_patterns_replayed=False,
        scope='Complete synthetic405group source/count/checkpoint/serialized workflow with two retained failures; '
              'four temporal patterns use genuine locked-oracle nonconstant, unknown/gap, constant and singular '
              'metrics. Synthetic source journals are not actual provenance. No native-parser or posterior qualification.')
    with args.output.open('x') as f: f.write(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(result))


if __name__ == '__main__': main()
