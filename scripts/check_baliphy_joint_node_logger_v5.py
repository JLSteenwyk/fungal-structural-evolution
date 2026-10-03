#!/usr/bin/env python3
"""Full source-grid and paired native CJSON logger qualification, including small alpha."""
import argparse
import copy
from datetime import datetime, timezone
import json
import math
from pathlib import Path

import numpy as np

from ancestral_chain_attempt import run_attempt, sha, write_json
from baliphy_joint_node_logger_v5 import transform, restore, numeric_transform, fresh_seeds, validate_frame
from independent_native_ancestral_alignment import fasta_records, native_alignments, tree_labels
from independent_joint_ancestral_frames import decode, write_arrays, verify_arrays
from independent_short_sampler_outputs import strict_json
from run_baliphy_reference_preflight import verify


def rejected(action):
    try:
        action()
    except (AssertionError, KeyError, ValueError):
        return
    raise AssertionError('Invalid logger or frame accepted')


def instrument(source):
    """Test-only separate rate encoder oracle; never included in future models."""
    replacements = [
        ('import qualified Data.JSON as J', 'import qualified SModel.Property as Property\nimport qualified Data.JSON as J'),
        ('[logA] [logCatStates] = do', '[logA] [logCatStates] [logEncoderChecks] = do'),
        (';return (parameterLogValues loggerValues ++',
         ';addLogger $ ((every 10) $ (logEncoderChecks (((J.toJSONKey "roundtrips") .= '\
         '[x == (read (Text.unpack (J.cjsonToText (J.toCJSON x))) :: Double) | '\
         'Property.StateProperties row <- Property.getComponentStateProperties '\
         '((prop_smodel_properties properties) Map.! (Text.pack "rate")), x <- row]) <> '\
         '((J.toJSONKey "originalProperties") .= (prop_smodel_properties properties)))))\n'\
         ';return (parameterLogValues loggerValues ++'),
        (';unless isTest (do', ';logEncoderChecks <- if loggingEnabled then ejsonLogger '\
         '(outputDirectory </> "encoder-checks.jsonl") else return noLogger\n;unless isTest (do'),
        ('[logA] [logCatStates])', '[logA] [logCatStates] [logEncoderChecks])'),
    ]
    for before, after in replacements:
        assert source.count(before) == 1 and after not in source
        source = source.replace(before, after)
    return source


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    a = p.parse_args()
    assert not a.receipt.exists()
    root = a.output.resolve(); root.mkdir(exist_ok=False)
    source = Path('metadata/baliphy_reference_sampler_qualification_plan_20261003.json')
    plan = json.loads(source.read_text()); verify(plan)
    jobs = json.loads(Path(plan['jobs']).read_text()); assert len(jobs) == 1620
    old_gate = Path('metadata/baliphy_joint_node_logger_software_validation_20261003_v4.json')
    gate = json.loads(old_gate.read_text()); verify({'pins': gate['source_hashes']})
    bindings = dict(plan['pins']); bindings.update(gate['source_hashes'])
    for path in [source, old_gate, Path(plan['jobs']), Path(__file__),
                 Path('scripts/baliphy_joint_node_logger_v5.py'),
                 Path('scripts/baliphy_joint_node_logger_v3.py'),
                 Path('scripts/baliphy_full_node_logger.py'),
                 Path('scripts/independent_joint_ancestral_frames.py'),
                 Path('scripts/independent_short_sampler_outputs_v2.py'),
                 Path('scripts/independent_native_ancestral_alignment.py')]:
        bindings[str(path)] = sha(path)
    programs = {j['chain']['program']: j['chain']['program_sha256'] for j in jobs}
    assert len(programs) == 405
    for name, digest in programs.items():
        assert sha(name) == digest
        text = Path(name).read_text()
        assert restore(transform(text)) == text
    first = Path(jobs[0]['chain']['program']).read_text()
    rejected(lambda: transform(transform(first)))
    from baliphy_joint_node_logger_v5 import PROPERTY_BEFORE
    rejected(lambda: transform(first.replace(PROPERTY_BEFORE, '')))
    rejected(lambda: transform(first+PROPERTY_BEFORE))
    used = {j['chain']['seed'] for j in jobs} | {j['source_seed'] for j in jobs}
    old_future = Path('metadata/baliphy_joint_node_logger_future_models_20261003_v4.json')
    future = json.loads(old_future.read_text()); verify(future)
    old_roles = Path(future['future_roles'])
    used.update(r['chain']['seed'] for r in json.loads(old_roles.read_text()))
    used.update({20267001, 20267002, 20267003, 20268001, 20268002, 20268003})
    namespace = 'fungal-joint-node-cjson-future-20261003-v5'
    ids = [j['chain']['chain_id'] for j in jobs]
    seeds = fresh_seeds(ids, used, namespace)
    assert seeds == fresh_seeds(list(reversed(ids)), used, namespace)
    assert len(seeds) == 1620 and set(seeds.values()).isdisjoint(used)
    for path in [old_future, old_roles]: bindings[str(path)] = sha(path)
    baseline = Path('data/software_audits/baliphy-joint-node-logger-20261003-v4').resolve()
    observed = {k: v.replace('-', '') for k, v in fasta_records((baseline/'alignment.faa').read_text().splitlines()).items()}
    binary = Path(jobs[0]['config']['command'][5]); prlimit = Path('/usr/bin/prlimit')
    api = binary.parent.parent/'lib/bali-phy/haskell'
    apis = [api/name for name in ['Graph.hs', 'Tree.hs', 'Bio/Alignment.hs', 'Bio/Alphabet.hs',
        'Probability/Distribution/PhyloCTMC/Properties.hs', 'SModel/ASRV.hs', 'SModel/MixtureModel.hs',
        'SModel/Property.hs', 'Probability/Logger.hs', 'Data/Text.hs', 'Data/OldList.hs',
        'Data/JSON.hs', 'Data/JSON/Encoding.hs', 'Data/JSON/Types/ToJSON.hs', 'Data/JSON/Types/Foreign.hs']]
    samples = []; pairs = []; failures = []; oracle_frames = []; corrupt_original = []

    def native(name, text, seed):
        program = root/(name+'.hs'); program.write_text(text)
        paths = [prlimit, binary, program, baseline/'alignment.faa', baseline/'tree.nwk', *apis]
        config = dict(command=[str(prlimit), '--as='+str(12*2**30), '--cpu=300',
            '--fsize='+str(64*2**20), '--', str(binary), '--seed', str(seed), 'run', str(program),
            '--iterations', '20', '--log-format', 'json,tsv', '--name', 'independent-chain'],
            timeout_seconds=300, pins={str(path): sha(path) for path in paths})
        receipt = run_attempt(root/'native'/name, config)
        result = json.loads(receipt.read_text()); assert result['exit_code'] == 0, str(receipt)
        folders = list(receipt.parent.glob('independent-chain-*')); assert len(folders) == 1
        bindings.update(config['pins']); bindings[str(receipt)] = sha(receipt)
        bindings[str(receipt.parent.parent/'configuration.json')] = sha(receipt.parent.parent/'configuration.json')
        for path, digest in result['artifacts'].items(): bindings[str(receipt.parent/path)] = digest
        print('cjson_logger_native', name, 'exited_zero', flush=True)
        return folders[0], str(receipt)

    def frames(directory, name):
        labels, tips = tree_labels((directory/'runtime-tree.nwk').read_text())
        records = [strict_json(line) for line in (directory/'C1.P1.site-property-samples.jsonl').read_text().splitlines()]
        assert [f['iter'] for f in records] == [0, 10, 20]
        legacy = list(native_alignments((directory/'C1.P1.fastas').read_text().splitlines()))
        for frame, (iteration, separate) in zip(records, legacy):
            sequence = fasta_records(frame['alignmentLines'])
            assert set(sequence) == labels and len(labels-tips) == 4
            assert {k: [c == '-' for c in v] for k, v in sequence.items()} == {k: [c == '-' for c in v] for k, v in separate.items()}
            checked = validate_frame(frame, iteration, sequence, tips)
            candidates = {node: node for node in sorted(labels-tips)}
            independent, arrays = decode(frame, iteration, observed, sorted(labels), candidates)
            assert checked['tip_pairs'] == independent['tip_pairs']
            assert checked['ancestral_pairs'] == independent['ancestral_pairs']
            path = root/(name+'-'+str(iteration)+'.npz')
            write_arrays(path, arrays); verify_arrays(path, arrays); bindings[str(path)] = sha(path)
            samples.append(dict(test=name, iteration=iteration, joint_check=checked, independent_check=independent))
            for case in ['rate_mean', 'rate_nonfinite', 'missing_ancestor', 'category_support', 'state_letter', 'coordinate_alignment']:
                bad = copy.deepcopy(frame); node = sorted(labels-tips)[0]
                if case == 'rate_mean': bad['properties']['rate'][0] = [2.]*20
                elif case == 'rate_nonfinite': bad['properties']['rate'][0] = [float('inf')]*20
                elif case == 'missing_ancestor': bad['catStates'].pop(node)
                elif case == 'category_support': bad['catStates'][node]['categories'][0] = 4
                elif case == 'state_letter': bad['catStates'][node]['states'][0] = (bad['catStates'][node]['states'][0]+1)%20
                else:
                    index = bad['alignmentLines'].index('>'+node)+1
                    value = bad['alignmentLines'][index]; offset = next(i for i, c in enumerate(value) if c != '-')
                    bad['alignmentLines'][index] = value[:offset]+('C' if value[offset] != 'C' else 'A')+value[offset+1:]
                rejected(lambda: validate_frame(bad, iteration, sequence, tips))
                rejected(lambda: decode(bad, iteration, observed, sorted(labels), candidates))
                failures.append(name+':'+str(iteration)+':'+case)
        return records

    for index, prior in enumerate(['broad', 'centered', 'package']):
        old_program = baseline/(prior+'-joint_node_logger.hs')
        old = next((baseline/'native'/prior/'joint_node_logger').glob('attempt-*/independent-chain-*'))
        text = numeric_transform(old_program.read_text())
        directory, receipt = native(prior+'-numeric-only', text, 20267001+index)
        for name in ['C1.log', 'C1.log.column-map.json', 'C1.log.json', 'runtime-tree.nwk', 'C1.P1.fastas']:
            assert (old/name).read_bytes() == (directory/name).read_bytes(), (prior, name)
        new = frames(directory, prior+'-numeric-only')
        before = [strict_json(line) for line in (old/'C1.P1.site-property-samples.jsonl').read_text().splitlines()]
        for left, right in zip(before, new):
            assert {k: v for k, v in left.items() if k != 'properties'} == {k: v for k, v in right.items() if k != 'properties'}
            assert np.allclose(left['properties']['rate'], right['properties']['rate'], rtol=2e-14, atol=0)
        pairs.append(dict(prior=prior, native_receipt=receipt, retained_baseline=str(old),
            exact_scalar_alignment_and_nonproperty_frame_identity=True,
            property_float_comparison='rtol2e-14_abs0_for_these_retained_nondamaged_fixture_values_only'))
        alpha = [0.06590801999395074, 0.03, 0.01][index]
        import re
        fixed, substitutions = re.subn(r';alpha_2 <- sample \([^\n]+\)', ';let {alpha_2 = '+repr(alpha)+'}', text)
        assert substitutions == 1
        small, small_receipt = native(prior+'-small-alpha', fixed, 20268001+index)
        instrumented, oracle_receipt = native(prior+'-small-alpha-oracle', instrument(fixed), 20268001+index)
        for name in ['C1.log', 'C1.log.column-map.json', 'C1.log.json', 'runtime-tree.nwk',
                     'C1.P1.fastas', 'C1.P1.site-property-samples.jsonl']:
            assert (small/name).read_bytes() == (instrumented/name).read_bytes(), ('Oracle changed output', prior, name)
        corrected = frames(small, prior+'-small-alpha')
        oracle = [strict_json(line) for line in (instrumented/'encoder-checks.jsonl').read_text().splitlines()]
        assert [f['iter'] for f in oracle] == [0, 10, 20]
        for frame, probe in zip(corrected, oracle):
            assert set(probe) == {'iter', 'roundtrips', 'originalProperties'}
            assert probe['roundtrips'] == [True]*80
            original_values = [row[0] for row in probe['originalProperties']['rate']]
            corrected_values = [row[0] for row in frame['properties']['rate']]
            mean = math.fsum(original_values)/4
            if not math.isclose(mean, 1., rel_tol=1e-10, abs_tol=1e-10):
                corrupt_original.append(dict(prior=prior, alpha=alpha, iteration=frame['iter'],
                    original_rates=original_values, original_mean=mean, corrected_rates=corrected_values))
            oracle_frames.append(dict(prior=prior, alpha=alpha, iteration=frame['iter'],
                exact_native_roundtrips=80, oracle_does_not_change_saved_outputs=True,
                native_receipt=small_receipt, oracle_receipt=oracle_receipt))
    assert len(samples) == 18 and len(oracle_frames) == 9 and corrupt_original
    assert sum(r['joint_check']['native_ancestors'] for r in samples) == 72
    assert len(failures) == 108
    for path, digest in bindings.items(): assert sha(path) == digest, path
    proof = dict(status='passed_joint_node_logger_v5_cjson_full_source_and_native_qualification',
        checked_utc=datetime.now(timezone.utc).isoformat(), full_programs_checked=405, full_roles_checked=1620,
        deterministic_disjoint_future_seeds=1620, forbidden_seed_count=len(used), future_seed_namespace=namespace,
        native_runs=9, paired_retained_prior_checks=pairs, same_record_frames_checked=18,
        ancestral_records_checked=72, independent_all_node_sequence_category_coordinate_checks=True,
        serialized_candidate_array_checks=18, strict_mean_one_assertion_unchanged=True,
        small_alpha_native_oracle_frames=oracle_frames, native_exact_double_roundtrips=720,
        original_encoder_normalization_failures_reproduced=corrupt_original,
        altered_frames_rejected_by_both_decoders=failures, source_hashes=bindings,
        test_instrumentation_in_future_models=False, scientific_eligibility=False, posterior_qualified=False,
        scope='All405programs exactly reversible,1620roles and new seeds retained; three corrected pure-rendering paired runs match retained scalar/alignment/nonproperty JSON bytes. Three fixed-small-alpha software runs plus same-seed separate native-oracle counterparts expose old rate exponent corruption and pass strict corrected normalization and720 exact native-double roundtrips. Two independently implemented frame checks plus full NPZ coordinate/dtype/value readback cover18frames/72ancestors;108 malformed cases rejected by both. Fixed alpha/oracle changes exist only in software fixtures. Full-grid native execution/startup, historical corruption handling, longer posterior/root/model/likelihood/memory adequacy remain unqualified. Existing sources/runs unchanged; no GPU/new charges.')
    with a.receipt.open('x') as handle: handle.write(json.dumps(proof, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: v for k, v in proof.items() if k != 'source_hashes'}, indent=2), flush=True)


if __name__ == '__main__': main()
