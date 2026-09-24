#!/usr/bin/env python3
"""Compare focal domain coordinate mappings under two MAFFT settings."""
import argparse
import json
from pathlib import Path
import subprocess
from Bio import SeqIO
from inspect_focal_domain_architectures import sha


def read_fasta(path):
    rows = list(SeqIO.parse(path, 'fasta'))
    assert len({r.id for r in rows}) == len(rows)
    return {r.id: str(r.seq).upper() for r in rows}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ['sequences', 'architectures', 'context', 'mafft_root', 'output']:
        ap.add_argument('--' + name.replace('_', '-'), type=Path, required=True)
    a = ap.parse_args()
    executable = a.mafft_root / 'bin/mafft'
    paths = [a.sequences, a.sequences.parent / 'receipt.json', a.architectures,
             a.context, Path(__file__), Path(__file__).with_name('inspect_focal_domain_architectures.py')]
    for subdir in ['bin', 'libexec', 'binaries']:
        paths.extend(p for p in (a.mafft_root / subdir).rglob('*') if p.is_file())
    pins = {str(p): sha(p) for p in paths}
    assert json.loads(paths[1].read_text())['artifacts'][a.sequences.name] == pins[str(a.sequences)]
    sequences = read_fasta(a.sequences)
    architectures = json.loads(a.architectures.read_text())
    context = json.loads(a.context.read_text())
    pairs = sorted({(f['label'], sister) for r in context['results']
                    for f in r['focal_context'] for sister in f['sister_labels']})
    domains = {}
    for record in architectures['records']:
        domains[record['tree_label']] = {h['hit_id']: h for alt in record['candidate_architectures']['alternatives']
                                         for h in alt['annotations'] if h['pfam_accession'].split('.')[0] == 'PF26973'}
    a.output.mkdir(parents=True, exist_ok=False)
    results = []; commands = []
    for mode in ['localpair', 'globalpair']:
        target = a.output / (mode + '.faa')
        cmd = [str(executable.resolve()), '--thread', '1', '--' + mode,
               '--maxiterate', '1000', str(a.sequences.resolve())]
        with target.open('w') as out, (a.output / (mode + '.log')).open('w') as log:
            subprocess.run(cmd, stdout=out, stderr=log, check=True, timeout=600)
        commands.append(cmd)
        alignment = read_fasta(target)
        assert set(alignment) == set(sequences)
        assert len({len(s) for s in alignment.values()}) == 1
        assert all(s.replace('-', '') == sequences[k] for k, s in alignment.items())
        positions = {}
        for label, seq in alignment.items():
            n = 0; pos = []
            for letter in seq:
                n += letter != '-'
                pos.append(n if letter != '-' else None)
            positions[label] = pos
        for focal, sister in pairs:
            for hit in domains[focal].values():
                for boundary in ['alignment', 'envelope']:
                    start, end = [hit[boundary + '_' + x] for x in ['start', 'end']]
                    matched = [(p, q) for p, q in zip(positions[focal], positions[sister])
                               if p is not None and start <= p <= end and q is not None]
                    identical = sum(sequences[focal][p-1] == sequences[sister][q-1] for p, q in matched)
                    results.append(dict(mode=mode, focal=focal, sister=sister, hit_id=hit['hit_id'],
                                        boundary=boundary, start=start, end=end, domain_residues=end-start+1,
                                        matched_residues=len(matched), identical_residues=identical,
                                        residue_pairs=matched))
    assert all(sha(p) == digest for p, digest in pins.items())
    result = dict(status='complete_focal_alignment_sensitivity_mapping', source_hashes=pins,
                  commands=commands, results=results,
                  artifacts={p.name: sha(p) for p in a.output.iterdir() if p.is_file()},
                  scope='Five-protein alignment sensitivity and explicit residue correspondences. '
                  'Not domain homology validation, structural comparison, or evidence of gain/loss. '
                  'Whole-protein alignments can misalign repeated domains.')
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps([{k: v for k, v in r.items() if k != 'residue_pairs'} for r in results], indent=2))


if __name__ == '__main__':
    main()
