#!/usr/bin/env python3
"""Merge independently closed full-genome evidence without a third corpus scan."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    paths = [Path('metadata/' + name + '_20261005_v1.json') for name in
             ['full_assembly_dna', 'full_assembly_dna_transport',
              'full_assembly_dna_readback', 'full_assembly_dna_readback_transport']]
    producer, producer_transport, reader, reader_transport = [json.loads(p.read_text()) for p in paths]
    assert producer['status'] == 'completed_all_selected_assembly_dna_attempts_pending_independent_readback'
    assert reader['status'] == 'passed_full_selected_genome_disposition_and_original_fasta_reconstruction'
    assert producer['expected_taxa'] == reader['taxa'] == 526
    assert producer['counts'] == {'verified_publisher_bound_genomic_dna': 526}
    assert reader['counts'] == {'independently_reconstructed_original_genome': 526}
    assert producer['verified_genome_bases'] == reader['verified_genome_bases'] == 26645610508
    assert len(reader['proofs']) == len({r['taxon_id'] for r in reader['proofs']}) == 526
    pins = {}
    original_waits = []
    for receipt_path, receipt, transport in [(paths[0], producer, producer_transport),
                                             (paths[2], reader, reader_transport)]:
        assert transport['status'] == 'verified_original_software_wait_and_whole_wrapper_payloads'
        assert transport['original_tool_terminal_exit_code'] == 0
        assert transport['validation_sha256'] == hashlib.sha256(receipt_path.read_bytes()).hexdigest()
        assert transport['whole_wrapper_initial_and_terminal_payloads_matched']
        assert transport['manager_start_records'] == transport['manager_completion_records'] == 1
        for mapping in [receipt['source_hashes'], transport['source_hashes']]:
            for path, digest in mapping.items():
                assert path not in pins or pins[path] == digest, path
                pins[path] = digest
        original_waits.append(dict(unit=transport['unit'], invocation_id=transport['invocation_id'],
                                   original_tool_session_id=transport['original_tool_session_id'],
                                   actual_terminal_exit_code=transport['original_tool_terminal_exit_code']))
    for path in [*paths, Path(__file__)]:
        name = str(path)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        assert name not in pins or pins[name] == digest, name
        pins[name] = digest
    result = dict(status='complete_verified_full_selected_assembly_dna',
        checked_utc=datetime.now(timezone.utc).isoformat(), taxa=526, verified_genomes=526,
        genome_acquisition_errors=0, verified_genome_bases=reader['verified_genome_bases'],
        original_transports=original_waits, complete_bound_files=len(pins), source_hashes=pins,
        scientific_eligibility=False, gpu=False, new_predictions=0, all_eight_aims_incomplete=True,
        scope='Compact merger of already independently reconstructed whole-genome source evidence and actual '
              'original API/native/whole-journal closures. No third corpus scan. Every original record, case '
              'and ambiguity is retained. This does not qualify genome-CDS agreement, translation, gene copies, '
              'species identity, contamination, selection or evolutionary results.')
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
