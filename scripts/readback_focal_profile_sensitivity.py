#!/usr/bin/env python3
"""Read focal HMMER domain tables with Bio.SearchIO and compare every recorded hit."""
import argparse
import json
from pathlib import Path
from Bio import SearchIO
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--search', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new immutable readback')
    receipt = checked_receipt(args.search)
    config_path = args.search / 'config.json'
    config = json.loads(config_path.read_text())
    if receipt['config_sha256'] != sha(config_path):
        raise ValueError('Changed search configuration')
    counts = {}
    for run in receipt['runs']:
        reconstructed = []
        for query in SearchIO.parse(args.search / (run['label'] + '.domtblout'), 'hmmsearch3-domtab'):
            ga = config['gathering_thresholds'][query.accession]
            for hit in query:
                for domain in hit.hsps:
                    reconstructed.append({'protein': hit.id, 'profile': query.accession,
                        'sequence_score': hit.bitscore, 'domain_score': domain.bitscore,
                        'sequence_evalue': hit.evalue, 'domain_independent_evalue': domain.evalue,
                        'hmm_start': domain.query_start + 1, 'hmm_end': domain.query_end,
                        'alignment_start': domain.hit_start + 1, 'alignment_end': domain.hit_end,
                        'envelope_start': domain.env_start + 1, 'envelope_end': domain.env_end,
                        'passes_sequence_and_domain_ga': hit.bitscore >= ga['sequence'] and domain.bitscore >= ga['domain']})
        if reconstructed != run['hits']:
            raise ValueError('Independent parser disagrees: ' + run['label'])
        counts[run['label']] = len(reconstructed)
    by_label = {r['label']: r['hits'] for r in receipt['runs']}
    same = by_label['standard_ga'] == by_label['max_ga']
    permissive_pass = [r for r in by_label['max_permissive'] if r['passes_sequence_and_domain_ga']]
    same_pass = permissive_pass == by_label['max_ga']
    result = {'status': 'passed_focal_profile_hit_readback', 'hit_rows_checked': counts,
              'standard_and_max_ga_identical': same,
              'permissive_hits_passing_ga_equal_max_ga': same_pass,
              'search_receipt_sha256': sha(args.search / 'receipt.json'),
              'script_sha256': sha(Path(__file__)),
              'scope': 'All producer-recorded hit fields checked with Bio.SearchIO, including coordinate conventions and gathering-threshold flags. Output artifacts rehashed. Same HMMER output, not an independent profile search or biological homology validation.'}
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
