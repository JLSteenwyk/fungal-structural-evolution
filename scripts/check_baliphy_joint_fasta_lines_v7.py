#!/usr/bin/env python3
"""Compare V7 row-wise joint logging against saved native results, without changing inference."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import run_attempt, sha
from baliphy_joint_fasta_lines_v7 import reverse, transform
from read_baliphy_scalar_json_v6b import compare_tsv
from reference_measurement_union_sources import bind, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    assert not args.receipt.exists()
    root = args.output.resolve()
    root.mkdir(exist_ok=False)
    qualified_path = Path('metadata/baliphy_scalar_json_v6_software_validation_20261004_v9.json')
    qualified = json.loads(qualified_path.read_text())
    pins = dict(qualified['source_hashes'])
    verify(pins)
    for path in [Path(__file__), Path('scripts/baliphy_joint_fasta_lines_v7.py'), qualified_path]:
        bind(pins, path)
    programs = sorted(Path('results/ancestral/prepared-full-scalar-json-v6-20261004-v1/models').glob('*.hs'))
    assert len(programs) == 405
    for program in programs:
        text = program.read_text()
        assert reverse(transform(text)) == text
        bind(pins, program)
    outcomes = []
    unchanged = ['C1.log', 'C1.log.json', 'C1.log.column-map.json', 'runtime-tree.nwk',
                 'C1.P1.fastas', 'C1.P1.site-property-samples.jsonl']
    for item in qualified['paired_prior_checks']:
        before = Path(item['new_directory'])
        old_config_path = before.parent.parent / 'configuration.json'
        old_config = json.loads(old_config_path.read_text())
        verify(old_config['pins'])
        command = list(old_config['command'])
        model_index = command.index('run') + 1
        original = Path(command[model_index])
        model = root / (item['prior'] + '-joint-fasta-v7.hs')
        source = original.read_text()
        model.write_text(transform(source))
        assert reverse(model.read_text()) == source
        command[model_index] = str(model)
        config = dict(command=command, timeout_seconds=old_config['timeout_seconds'],
                      pins={**old_config['pins'], str(model): sha(model)})
        target = root / 'native' / item['prior']
        assert not target.exists()
        receipt_path = run_attempt(target, config)
        receipt = json.loads(receipt_path.read_text())
        assert receipt['exit_code'] == 0, receipt_path
        after = receipt_path.parent / 'independent-chain-1'
        for name in unchanged:
            assert (after / name).read_bytes() == (before / name).read_bytes(), (item['prior'], name)
        comparison = compare_tsv(after)
        assert comparison['rows'] == 21 and comparison['mapped_values_compared'] == 903
        frames = [json.loads(line) for line in (after/'C1.P1.site-property-samples.jsonl').read_text().splitlines()]
        assert [frame['iter'] for frame in frames] == [0,10,20]
        pins.update(config['pins'])
        for path in [old_config_path, target/'configuration.json', receipt_path]:
            bind(pins, path)
        for name, digest in receipt['artifacts'].items():
            pins[str(receipt_path.parent/name)] = digest
        outcomes.append(dict(prior=item['prior'],seed=item['seed'],native_receipt=str(receipt_path),
                             byte_identical_files=unchanged,joint_frames=len(frames),**comparison))
        print('joint_fasta_v7_native_control',item['prior'],'passed',flush=True)
    verify(pins)
    result = dict(status='passed_reversible_rowwise_joint_fasta_v7_paired_native_controls',
                  checked_utc=datetime.now(timezone.utc).isoformat(),full_programs_reversibly_checked=405,
                  native_prior_controls=3,native_scalar_rows=63,native_mapped_comparisons=2709,
                  native_joint_frames=9,outcomes=outcomes,source_hashes=pins,
                  probability_model_changed=False,installed_software_changed=False,
                  original_jobs_restarted=False,scientific_eligibility=False,posterior_qualified=False,
                  scope='Pure logger construction replaces the whole-alignment Text.unpack/lines/pack '
                        'roundtrip with two native Text rows per sequence, sharing the same ancestral data. '
                        'Three prior controls preserve all six scientific files byte-for-byte; '
                        'all405full models reverse exactly to original source. This is a candidate '
                        'software correction, not a validated full-input crash repair or posterior.')
    with args.receipt.open('x') as handle:
        json.dump(result,handle,indent=2,allow_nan=False)
        handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','outcomes']},indent=2))


if __name__ == '__main__':
    main()
