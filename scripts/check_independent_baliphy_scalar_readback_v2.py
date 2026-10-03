#!/usr/bin/env python3
"""Check the complete 405-group serialized comparison and checkpoint contracts."""
import argparse
from collections import Counter
import copy
import json
import os
from pathlib import Path
import tempfile

from prepare_independent_baliphy_scalar_readback_v2 import run, atomic
from readback_independent_baliphy_scalar_readback_v2 import run as readback
from run_ortholog_pair_guide_comparison import sha


def setup(root):
    source = root / 'source'; source.mkdir()
    native_log = source / 'native.tsv'; length_log = source / 'length.tsv'
    variables = ['variable-' + str(i) for i in range(37)]
    lengths = ['length_level' + str(i) for i in range(4)]
    with native_log.open('w') as f:
        f.write('iter\t' + '\t'.join(variables) + '\n')
        for i in range(1001):
            f.write(str(i) + '\t' + '\t'.join(['1'] * 37) + '\n')
    with length_log.open('w') as f:
        f.write('iter\t' + '\t'.join(lengths) + '\n')
        for i in range(0, 1001, 10):
            f.write(str(i) + '\t' + '\t'.join(['3'] * 4) + '\n')
    bindings = {str(native_log): sha(native_log), str(length_log): sha(length_log)}
    def create(path, value):
        atomic(path, value); bindings[str(path)] = sha(path); return path
    groups = {}; overlay = []; statuses = Counter()
    for number in range(405):
        group = 'group-' + str(number).zfill(4); folder = source / group; folder.mkdir()
        ids = [group + '-chain' + str(i) for i in range(4)]; failed = number >= 403
        chains = {kind: [] for kind in ['scalar', 'length']}
        for i, cid in enumerate(ids):
            chain = dict(chain_id=cid, seed=number * 4 + i + 1)
            disposition = dict(status='failed' if failed and i == 0 else 'all_saved_alignments_and_candidate_nodes_checked')
            for kind, template in [('scalar', native_log), ('length', length_log)]:
                log = folder / (cid + '-' + kind + '.tsv'); os.link(template, log); bindings[str(log)] = sha(log)
                chains[kind].append(dict(chain_id=cid, seed=chain['seed'], model_input_identity=group,
                    log=str(log), log_sha256=sha(log)))
            audit = create(folder / (cid + '-audit.json'), dict(scalar_log=chains['scalar'][-1]['log'],
                scalar_log_sha256=chains['scalar'][-1]['log_sha256']))
            disposition.update(sample_audit=str(audit), sample_audit_sha256=sha(audit))
            overlay.append(dict(chain=chain, model_input_identity=group, selected_disposition=disposition))
        info = dict(chain_ids=ids, scientific_eligibility=False,
            status='unresolved_failed_native_chain_retained' if failed else 'complete_scalar_length_category_screens_not_posterior_qualification')
        if not failed:
            for kind, names, grid in [('scalar', variables, list(range(1001))), ('length', lengths, list(range(0, 1001, 10)))]:
                target = folder / kind; target.mkdir(); outputs = {}
                for cutoff in ['250', '500']:
                    mp = create(target / ('manifest-' + cutoff + '.json'), dict(chains=chains[kind],
                        variables=names, expected_iterations=grid, discard_through_iteration=int(cutoff)))
                    draws = (750 if cutoff == '250' else 500) if kind == 'scalar' else (75 if cutoff == '250' else 50)
                    reports = {name: dict(status='constant_chain_requires_review', draws_per_chain=draws) for name in names}
                    fp = create(target / ('report-' + cutoff + '.json'), dict(status='scalar_diagnostics_complete_not_posterior_qualification',
                        arviz_version='0.22.0', numpy_version='2.2.6', discard_through_iteration=int(cutoff), variables=reports,
                        thresholds=dict(rhat_strict_upper=1.01, bulk_ess_minimum=400, tail_ess_minimum=400),
                        pins={str(mp): sha(mp), **{r['log']: r['log_sha256'] for r in chains[kind]}}))
                    outputs[cutoff] = dict(path=str(fp), sha256=sha(fp))
                    statuses[kind + ':' + cutoff + ':constant_chain_requires_review'] += len(names)
                rp = create(target / 'receipt.json', dict(outputs=outputs))
                info[kind] = dict(receipt=str(rp), receipt_sha256=sha(rp))
        groups[group] = info
    op = create(source / 'overlay.json', dict(rows=overlay, original_samples_concatenated=False, scientific_eligibility=False))
    recovery = create(source / 'recovery-completed.json', dict(overlay=str(op), overlay_sha256=sha(op)))
    producer = create(source / 'receipt.json', dict(groups=groups))
    summary = dict(full_quartets=405, complete_quartets=403, unresolved_quartets=2,
        scalar_status_counts={k.split(':', 1)[1]: v for k, v in statuses.items() if k.startswith('scalar:')},
        length_status_counts={k.split(':', 1)[1]: v for k, v in statuses.items() if k.startswith('length:')},
        quartets_passing_every_scalar={'250': 0, '500': 0}, quartets_passing_every_length_scalar={'250': 0, '500': 0})
    archive = source / 'archive.json'
    atomic(archive, dict(status='complete_verified_full_baliphy_recovery_diagnostics_archive', services=[{}, {}],
        summary=summary, source_hashes=bindings))
    complete = source / 'diagnostic-completed.json'
    atomic(complete, dict(status='complete_verified_full_baliphy_recovery_diagnostics', **summary,
        exact_process_journals_checked=2, bound_source_hashes=len(bindings), full_hash_archive=str(archive),
        full_hash_archive_sha256=sha(archive), producer_receipt=str(producer), producer_receipt_sha256=sha(producer)))
    path = root / 'plan.json'; atomic(path, dict(pins={}, diagnostic_completion=str(complete),
        recovery_completion=str(recovery), output=str(root / 'output'), comparison_tolerances=dict(atol=1e-8, rtol=1e-8),
        resources=dict(minimum_free_disk_gib=0), scope='Complete synthetic source/count/provenance contracts; source journal fixtures are synthetic.'))
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); rejected = []
    with tempfile.TemporaryDirectory(prefix='independent-baliphy-scalar-contract-') as directory:
        root = Path(directory); path = setup(root)
        try:
            run(path, stop_after_groups=3)
        except InterruptedError:
            pass
        else:
            raise AssertionError('Interruption not exercised')
        before = {str(p): sha(p) for p in (root / 'output/groups').glob('*.json')}
        assert len(before) == 3
        result = run(path); assert result['variable_cutoff_rows'] == 33046
        assert all(sha(p) == d for p, d in before.items())
        readback(path, root / 'output/readback.json')
        for action in ['producer', 'alternate_reader']:
            try:
                run(path) if action == 'producer' else readback(path, root / 'alternate-reader.json')
            except AssertionError:
                pass
            else:
                raise AssertionError('Completed stage restarted: ' + action)
        (root / 'output/readback.json').unlink()
        manifest = root / 'output/group_manifest.json'; receipt = root / 'output/receipt.json'
        original_manifest = manifest.read_bytes(); original_receipt = receipt.read_bytes()
        first = json.loads(original_manifest)[0]; cp = root / 'output' / first['path']; original = cp.read_bytes()
        for change in ['remove_variable', 'duplicate_variable', 'change_metric', 'change_source_disposition', 'foreign_report', 'promote_science', 'hide_unresolved_group']:
            altered = json.loads(original); row = altered['result']['rows'][0]
            if change == 'remove_variable':
                altered['result']['rows'].pop()
            elif change == 'duplicate_variable':
                altered['result']['rows'][-1] = copy.deepcopy(row)
            elif change == 'change_metric':
                row['independent']['rhat'] = .5
            elif change == 'change_source_disposition':
                row['original']['status'] = 'passes_scalar_screen_only'
            elif change == 'foreign_report':
                row['original_report_sha256'] = '0' * 64
            elif change == 'promote_science':
                row['scientific_eligibility'] = True
            else:
                altered['result']['status'] = 'unresolved_failed_native_chain_retained'; altered['result']['rows'] = []
            atomic(cp, altered)
            parts = json.loads(original_manifest); parts[0]['sha256'] = sha(cp); atomic(manifest, parts)
            document = json.loads(original_receipt)
            document['artifacts'][first['path']] = sha(cp); document['artifacts']['group_manifest.json'] = sha(manifest)
            atomic(receipt, document)
            try:
                readback(path, root / 'malformed-readback.json')
            except AssertionError:
                rejected.append(change)
            else:
                raise AssertionError('Rehashed false comparison accepted: ' + change)
            cp.write_bytes(original); manifest.write_bytes(original_manifest); receipt.write_bytes(original_receipt)
    paths = [Path(__file__), *[Path('scripts/' + name + '.py') for name in ['independent_ancestral_scalar_diagnostics_v2',
        'independent_baliphy_scalar_sources', 'prepare_independent_baliphy_scalar_readback_v2', 'readback_independent_baliphy_scalar_readback_v2']]]
    result = dict(status='passed_full_independent_baliphy_scalar_grid_and_checkpoint_contracts',
        synthetic_full_quartets=405, complete_quartets=403, unresolved_quartets=2,
        variable_cutoff_rows=33046, unchanged_interrupted_checkpoints_reused=True,
        completed_producer_and_alternate_reader_restart_refused=True, rehashed_false_exports_rejected=rejected,
        source_hashes={str(p): sha(p) for p in paths}, scientific_eligibility=False,
        scope='Full synthetic constant-trace grid/count/serialized/provenance/checkpoint contracts, retaining every405group/two cutoffs/37scalar+4lengthvariables. Source journals and constant diagnostics are fixtures; separate106case oracle contracts test actual nonconstant numerical methods. No production posterior or biological pilot.')
    with args.output.open('x') as f:
        f.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
