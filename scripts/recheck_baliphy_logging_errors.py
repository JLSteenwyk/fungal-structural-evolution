#!/usr/bin/env python3
"""Recheck corrected retained frames and freshly reproduce native numeric encoding."""
import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import resource

from ancestral_chain_attempt import run_attempt, sha
from baliphy_joint_node_logger_v3 import validate_frame
from independent_native_ancestral_alignment import fasta_records, tree_labels
from independent_short_sampler_outputs import strict_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    assert not args.receipt.exists()
    root = args.output.resolve()
    root.mkdir(exist_ok=False)
    software = Path('metadata/baliphy_joint_node_logger_software_validation_20261003_v4.json')
    previous = json.loads(software.read_text())
    bindings = {str(software): sha(software), str(Path(__file__)): sha(__file__)}
    for path, digest in previous['source_hashes'].items():
        assert sha(path) == digest, path
        bindings[path] = digest
    frames = []
    for name in previous['source_hashes']:
        if '/joint_node_logger/' not in name or not name.endswith('site-property-samples.jsonl'):
            continue
        path = Path(name)
        labels, tips = tree_labels((path.parent/'runtime-tree.nwk').read_text())
        records = [strict_json(line) for line in path.read_text().splitlines()]
        assert [f['iter'] for f in records] == [0, 10, 20]
        for frame in records:
            sequences = fasta_records(frame['alignmentLines'])
            assert set(sequences) == labels
            frames.append(validate_frame(frame, frame['iter'], sequences, tips))
    assert len(frames) == 9
    assert sum(f['native_ancestors'] for f in frames) == 36
    assert sum(f['ancestral_pairs'] for f in frames) == 233
    binary = Path('data/software_audits/baliphy-4.3-20260927/install/bali-phy-4.3/bin/bali-phy').resolve()
    api = binary.parent.parent/'lib/bali-phy/haskell'
    library_paths = [api/name for name in ['Data/JSON.hs', 'Data/JSON/Encoding.hs',
        'Data/JSON/Types/ToJSON.hs', 'Data/JSON/Types/Foreign.hs', 'Data/Text/Display.hs', 'Data/Text.hs']]
    probes = {}
    for version, source in [('original', 'original.hs'), ('cjson', 'cjson-native-roundtrip-v2.hs')]:
        source = Path('config/native-format-probes')/source
        program = root/(version+'.hs')
        program.write_bytes(source.read_bytes())
        paths = [Path('/usr/bin/prlimit'), binary, program, *library_paths]
        config = dict(command=[str(paths[0]), '--as='+str(8*2**30), '--cpu=20',
            '--fsize='+str(16*2**20), '--', str(binary), 'run', str(program)],
            timeout_seconds=30, pins={str(p): sha(p) for p in paths})
        receipt = run_attempt(root/version, config)
        result = json.loads(receipt.read_text())
        assert result['exit_code'] == 0
        lines = (receipt.parent/'stdout.log').read_text().splitlines()
        values = json.loads(lines[0])
        expected = [2.34e-10, 2.34e-11, 2.34e-20, 2.34e-50, 2.34e-100,
                    2.34e-101, 2.34e10, 2.34e20, 0., 1., .25, 4.]
        assert len(values) == 12
        if version == 'original':
            assert len(lines) == 1 and values[0] == .234 and values[2] == .0234
            assert values[7] == 234 and values[8:] == expected[8:]
            wrong = [i for i, (v, e) in enumerate(zip(values, expected)) if v != e]
            assert len(wrong) == 5
            probes[version] = dict(error_reproduced=True, incorrect_constant_indices=wrong, values=values)
        else:
            assert len(lines) == 2 and json.loads(lines[1]) == [True]*12
            assert all(math.isclose(v, e, rel_tol=2e-15, abs_tol=0) for v, e in zip(values, expected))
            probes[version] = dict(native_exact_roundtrips=12, values=values,
                                  full_model_logger_fix_qualified=False)
        probes[version]['native_receipt'] = str(receipt)
        bindings.update(config['pins'])
        bindings[str(source)] = sha(source)
        bindings[str(receipt)] = sha(receipt)
        for name, digest in result['artifacts'].items():
            bindings[str(receipt.parent/name)] = digest
    for path, digest in bindings.items():
        assert sha(path) == digest, path
    proof = dict(status='completed_fresh_native_and_retained_logger_error_recheck',
        checked_utc=datetime.now(timezone.utc).isoformat(),
        corrected_fixture_bindings_checked=len(previous['source_hashes']),
        corrected_saved_frames_checked=9, corrected_ancestral_records_checked=36,
        corrected_ancestral_residue_category_pairs=233,
        corrected_same_record_sequence_state_mismatches=0, numeric_probes=probes,
        source_hashes=bindings, new_mcmc_runs=0, existing_jobs_restarted=False,
        child_cpu_seconds=resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime +
            resource.getrusage(resource.RUSAGE_CHILDREN).ru_stime,
        scope='Retained corrected sequence/state frames pass. Installed original numeric formatter still fails; CJSON exact roundtrip passes only the twelve fixed native constants. Full future logger correction, historical float assessment and posterior qualification remain pending. No installed binary/library or old output changes.')
    with args.receipt.open('x') as handle:
        handle.write(json.dumps(proof, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: v for k, v in proof.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
