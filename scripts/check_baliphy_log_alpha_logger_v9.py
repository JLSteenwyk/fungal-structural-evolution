#!/usr/bin/env python3
"""Qualify captured latent alpha without changing existing short-run inference outputs."""
import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path

from ancestral_chain_attempt import run_attempt, sha
from baliphy_joint_fasta_lines_v7 import transform as joint_v7
from baliphy_log_alpha_logger_v9 import transform, reverse, fields
from read_baliphy_scalar_json_v6b import load, read_record, compare_tsv
from reference_measurement_union_sources import bind, verify


PRIORS = {'broad': (0, 2), 'centered': (0, 1), 'package': (6, 2)}
GRID = [-20., -10., -6., -3., 0., 6., 20., 50., 100., 500., 700., 709., 709.7,
        709.78, 709.79, 710., 750., 1000., 5000.]


def check_diagnostic(record, mu, scale, original=None):
    audit = read_record(record)
    state = record['parameters//']['S1/']
    assert set(state) == {'latentLogAlpha', 'derivedAlpha', 'categoryRates',
                         'laplaceLocation', 'laplaceScale', 'latentLogDensity'}
    x = state['latentLogAlpha']
    assert type(x) in (float, int) and math.isfinite(x)
    assert state['laplaceLocation'] == mu and state['laplaceScale'] == scale
    density = -math.log(2*scale)-abs(x-mu)/scale
    assert math.isclose(state['latentLogDensity'], density, rel_tol=2e-13, abs_tol=2e-13)
    try:
        alpha = math.exp(x)
    except OverflowError:
        alpha = math.inf
    if math.isinf(alpha):
        assert state['derivedAlpha'] == '__project_scalar_v6__:positive_infinity'
        assert audit['nonfinite_reviews'] == [dict(section='parameters',
            path=['S1/', 'derivedAlpha'], kind='positive_infinity')]
    else:
        assert math.isclose(state['derivedAlpha'], alpha, rel_tol=2e-13, abs_tol=0)
        assert not audit['nonfinite_reviews']
    rates = state['categoryRates']
    assert len(rates) == 4 and all(type(v) in (int, float) and math.isfinite(v) and v >= 0 for v in rates)
    assert math.isclose(math.fsum(rates)/4, 1., rel_tol=1e-10, abs_tol=1e-10)
    if x >= 100:
        assert rates == [1., 1., 1., 1.]
    if original is not None:
        assert record['iter'] == original['iter']
        assert record['statistics//'] == original['statistics//']
        assert state['derivedAlpha'] == original['parameters//']['S1/']['ASRV.Gamma:alpha']
    return bool(math.isinf(alpha))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists() and not args.receipt.exists()
    root = args.output.resolve()
    root.mkdir(parents=True)
    gate = Path('metadata/baliphy_joint_fasta_v7_software_validation_20261004_v1.json')
    qualified = json.loads(gate.read_text())
    assert qualified['status'] == 'passed_reversible_rowwise_joint_fasta_v7_paired_native_controls'
    pins = dict(qualified['source_hashes'])
    verify(pins)
    for p in [gate, Path(__file__), Path('scripts/baliphy_log_alpha_logger_v9.py')]:
        bind(pins, p)
    programs = sorted(Path('results/ancestral/prepared-full-scalar-json-v6-20261004-v1/models').glob('*.hs'))
    assert len(programs) == 405
    for program in programs:
        source = joint_v7(program.read_text())
        candidate = transform(source)
        assert reverse(candidate) == source
        assert candidate.count(';addLogger') == source.count(';addLogger')
        bind(pins, program)
    unchanged = ['C1.log', 'C1.log.json', 'C1.log.column-map.json', 'runtime-tree.nwk',
                 'C1.P1.fastas', 'C1.P1.site-property-samples.jsonl']
    outcomes = []
    native_library = None
    for item in qualified['outcomes']:
        previous_receipt = Path(item['native_receipt'])
        original_config = previous_receipt.parent.parent/'configuration.json'
        config = json.loads(original_config.read_text())
        verify(config['pins'])
        command = list(config['command'])
        index = command.index('run') + 1
        source = Path(command[index]).read_text()
        model = root/(item['prior']+'-log-alpha-v9.hs')
        model.write_text(transform(source))
        assert reverse(model.read_text()) == source
        command[index] = str(model)
        config = dict(command=command, timeout_seconds=config['timeout_seconds'],
                      pins=dict(config['pins'], **{str(model):sha(model)}))
        native_library = Path(command[command.index('--')+1]).parent.parent/'lib/bali-phy/haskell'
        transform_source = native_library/'Probability/Distribution/Transform.hs'
        text = transform_source.read_text()
        assert 'sample (ExpTransform dist) = exp <$> sample dist' in text
        assert 'logLaplace m s = ExpTransform $ laplace m s' in text
        target = root/'native'/item['prior']
        receipt_path = run_attempt(target, config)
        receipt = json.loads(receipt_path.read_text())
        assert receipt['exit_code'] == 0, receipt_path
        before = previous_receipt.parent/'independent-chain-1'
        after = receipt_path.parent/'independent-chain-1'
        for name in unchanged:
            assert (before/name).read_bytes() == (after/name).read_bytes(), (item['prior'], name)
        scalar = compare_tsv(after)
        assert scalar['rows'] == 21 and scalar['mapped_values_compared'] == 903
        old = [load(line) for line in (after/'C1.log.json').read_text().splitlines()[1:]]
        diagnostics = [load(line) for line in (after/'C1.P1.log-alpha-samples.jsonl').read_text().splitlines()[1:]]
        assert len(diagnostics) == 21
        mu, scale = PRIORS[item['prior']]
        for row, reference in zip(diagnostics, old):
            assert not check_diagnostic(row, mu, scale, reference)
        for p, digest in config['pins'].items():
            bind(pins, p, digest)
        for p in [original_config, target/'configuration.json', receipt_path]:
            bind(pins, p)
        for p, digest in receipt['artifacts'].items():
            bind(pins, receipt_path.parent/p, digest)
        outcomes.append(dict(prior=item['prior'], seed=item['seed'], native_receipt=str(receipt_path),
                             byte_identical_files=unchanged, diagnostic_rows=21, **scalar))
        print('log_alpha_v9_native_control', item['prior'], 'passed', flush=True)
    # A deterministic diagnostic writer exercises the same generated fields,
    # explicit special-value scalar encoding and actual native Gamma routine.
    source = (root/'broad-log-alpha-v9.hs').read_text()
    prefix = source[:source.index('sample_smodel alpha  =')]
    probe = prefix + '\nprobe priorNumber mu scale x iteration = do\n'
    probe += '  let projectLatentLogAlphaV9 = x :: Double\n      alpha_2 = exp x :: Double\n'
    # The field expression is taken directly from the qualified transformation.
    expression = fields('0', '2').removeprefix(';let {projectAlphaFieldsV9 = ').removesuffix('}\n')
    expression = expression.replace('(0 :: Double)', 'mu').replace('(2 :: Double)', 'scale')
    expression = expression.replace('(laplace 0 2)', '(laplace mu scale)')
    probe += '      projectAlphaFieldsV9 = ' + expression + '\n'
    probe += '      logged = ["S1" %>% projectAlphaFieldsV9]\n'
    probe += '  T.putStrLn (J.fromEncoding (projectV6EncodeRecord iteration (projectV6ContextValue []) logged))\n'
    probe += 'main = sequence_ [probe p mu scale x (p*19+n) | (p,mu,scale) <- [(0,0,2),(1,0,1),(2,6,2)], (n,x) <- zip [0..] ' + repr(GRID) + ']\n'
    probe_path = root/'deterministic-log-alpha-v9.hs'
    probe_path.write_text(probe)
    binary = native_library.parents[2]/'bin/bali-phy'
    paths = [Path('/usr/bin/prlimit'), binary, probe_path, *sorted(native_library.rglob('*.hs'))]
    config = dict(command=['/usr/bin/prlimit', '--as='+str(12*2**30), '--cpu=300', '--fsize='+str(64*2**20),
                           '--', str(binary), 'run', str(probe_path)], timeout_seconds=300,
                  pins={str(p):sha(p) for p in paths})
    receipt_path = run_attempt(root/'deterministic', config)
    receipt = json.loads(receipt_path.read_text())
    assert receipt['exit_code'] == 0, receipt_path
    records = [load(line) for line in (receipt_path.parent/'stdout.log').read_text().splitlines()]
    assert len(records) == 57 and [r['iter'] for r in records] == list(range(57))
    overflows = 0
    for label, offset in [('broad', 0), ('centered', 19), ('package', 38)]:
        mu, scale = PRIORS[label]
        subset = records[offset:offset+19]
        assert [row['parameters//']['S1/']['latentLogAlpha'] for row in subset] == GRID
        overflows += sum(check_diagnostic(row, mu, scale) for row in subset)
    assert overflows == 15
    for p, digest in config['pins'].items():
        bind(pins, p, digest)
    bind(pins, root/'deterministic/configuration.json')
    bind(pins, receipt_path)
    for p, digest in receipt['artifacts'].items():
        bind(pins, receipt_path.parent/p, digest)
    verify(pins)
    result = dict(status='passed_reversible_latent_log_alpha_v9_full_source_and_paired_native_controls',
        checked_utc=datetime.now(timezone.utc).isoformat(), full_programs_reversibly_checked=405,
        paired_prior_controls=3, scalar_rows=63, unchanged_scalar_comparisons=2709,
        new_diagnostic_rows=63, deterministic_diagnostic_rows=57, explicit_overflow_rows=15,
        outcomes=outcomes, source_hashes=pins, installed_software_changed=False,
        original_jobs_restarted=False, prior_or_likelihood_changed=False,
        scientific_eligibility=False, posterior_qualified=False, original_chain_overflow_cause_proven=False,
        scope='Expose the same latent Laplace draw used by the installed exp-transform, retain derived '
              'alpha and native category rates in a separate corrected scalar file. All405sources reverse '
              'exactly; three paired short native controls preserve six original scientific files byte for '
              'byte.57deterministic rows retain15overflow states. No historical latent value is reconstructed, '
              'existing review array admitted, prior truncated, native binary changed or posterior accepted.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes', 'outcomes']}, indent=2))


if __name__ == '__main__':
    main()
