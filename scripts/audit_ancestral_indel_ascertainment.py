#!/usr/bin/env python3
"""Inventory every SIC character before selecting an ascertainment model."""
import csv
import hashlib
import json
from pathlib import Path
from Bio import SeqIO


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    root = Path('results/ancestral/indel-coding-20260927-v1')
    receipt = json.loads((root / 'receipt.json').read_text())
    assert receipt['encodings'] == 156
    out = Path('results/ancestral/indel-ascertainment-audit-20260927-v1')
    out.mkdir(parents=True, exist_ok=False)
    summaries, characters, pins = [], [], {}
    for rel, digest in sorted(receipt['job_receipts'].items()):
        rp = root / rel
        assert sha(rp) == digest
        job = json.loads(rp.read_text())
        fp = rp.parent / 'characters.faa'
        assert sha(fp) == job['artifacts']['characters.faa']
        pins[str(rp)] = digest
        pins[str(fp)] = sha(fp)
        records = list(SeqIO.parse(fp, 'fasta'))
        assert len({r.id for r in records}) == len(records)
        seqs = [str(r.seq) for r in records]
        assert len(set(map(len, seqs))) == 1
        assert set(''.join(seqs)) <= set('01?')
        assert len(seqs) == job['summary']['proteins']
        assert len(seqs[0]) == job['summary']['characters']
        rows = []
        for i, states in enumerate(zip(*seqs), 1):
            mask = ''.join('?' if s == '?' else '.' for s in states)
            row = dict(encoding=rp.parent.name, character=i,
                       observed_zero=states.count('0'), observed_one=states.count('1'),
                       unknown=states.count('?'),
                       missing_mask_sha256=hashlib.sha256(mask.encode()).hexdigest())
            rows.append(row)
        characters.extend(rows)
        summaries.append(dict(encoding=rp.parent.name, proteins=len(seqs),
            characters=len(rows), characters_with_unknown=sum(r['unknown'] > 0 for r in rows),
            characters_without_observed_one=sum(r['observed_one'] == 0 for r in rows),
            characters_without_observed_zero=sum(r['observed_zero'] == 0 for r in rows),
            distinct_unknown_masks=len({r['missing_mask_sha256'] for r in rows}),
            disposition='no_coded_characters' if not rows else 'requires_explicit_ascertainment_assumption'))
    assert len(summaries) == 156 and len(characters) == 12957
    for filename, rows in [('encoding_summary.tsv', summaries), ('character_observability.tsv', characters)]:
        with (out / filename).open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)
    source = Path('data/software_audits/fastml-3.11/source/FastML.v3.11')
    source_files = ['libs/phylogeny/sequenceContainer.cpp', 'libs/phylogeny/unObservableData.cpp',
                    'libs/phylogeny/likelihoodComputation.cpp', 'programs/gainLoss/gainLoss.cpp']
    result = dict(status='complete_all_156_character_observability_audit',
        encodings=len(summaries), characters=len(characters),
        characters_with_unknown=sum(r['unknown'] > 0 for r in characters),
        characters_without_observed_one=sum(r['observed_one'] == 0 for r in characters),
        empty_encodings=[r['encoding'] for r in summaries if not r['characters']],
        source_receipt_sha256=sha(root / 'receipt.json'), input_hashes=pins,
        software_source_hashes={str(source / p):sha(source / p) for p in source_files},
        script_sha256=sha(__file__),
        artifacts={p.name:sha(p) for p in out.iterdir() if p.is_file()},
        scope='Character observability only. No fitted indel model or ancestral posterior. Missing states are coding-dependent; a mask-conditional correction alone does not establish an adequate generative model.')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['input_hashes','software_source_hashes','artifacts']}))


if __name__ == '__main__':
    main()
