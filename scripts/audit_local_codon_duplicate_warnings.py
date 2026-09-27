#!/usr/bin/env python3
"""Resolve all fit warnings and inventory identical codon sequences by taxon."""
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    fits = Path('results/cds/local-mg94-diagnostics-20260927-v1')
    audit = Path('results/cds/local-mg94-audit-20260927-v1')
    receipt = json.loads((fits/'receipt.json').read_text())
    ar = json.loads((audit/'receipt.json').read_text())
    assert ar['status'] == 'passed_saved_fit_readback' and ar['audited_cases'] == 1632
    with (audit/'case_audit.tsv').open() as f: summary = {r['case_id']:r for r in csv.DictReader(f,delimiter='\t')}
    assert sha(audit/'case_audit.tsv') == ar['artifacts']['case_audit.tsv']
    assert set(summary) == {r['case_id'] for r in receipt['case_receipts']}
    rows, groups = [], []
    for item in receipt['case_receipts']:
        case = item['case_id'];folder = fits/case
        assert sha(folder/'receipt.json') == item['receipt_sha256']
        cr = json.loads((folder/'receipt.json').read_text())
        assert sha(folder/'fit.log') == cr['artifacts']['fit.log']
        assert sha(folder/'config.json') == cr['config_sha256']
        config = json.loads((folder/'config.json').read_text())
        command = config['command'];alignment = Path(command[command.index('--alignment')+1])
        assert sha(alignment) == config['alignment_sha256']
        records = {};key = None
        for line in alignment.read_text().splitlines():
            if line.startswith('>'):
                key = line[1:].split()[0];assert key not in records;records[key]=''
            else:
                assert key is not None;records[key] += line.strip().upper()
        by_sequence = defaultdict(list)
        for taxon, sequence in records.items():by_sequence[sequence].append(taxon)
        duplicate_count = len(records)-len(by_sequence)
        log = (folder/'fit.log').read_text()
        warning_lines = [line for line in log.splitlines() if 'warning' in line.lower()]
        blocks = re.findall(r'>\[WARNING\]\n(.*?)(?=\n-------|\Z)',log,re.S)
        assert len(blocks) == len(warning_lines) == int(summary[case]['warning_lines'])
        reported = 0
        for block in blocks:
            matches = re.findall(r'contains (\d+) duplicate sequences?\.',block)
            assert len(matches)==1 and 'Identical sequences do not contribute' in block and 'prior to running selection analyses' in block, block
            reported += int(matches[0])
        assert reported == duplicate_count, (case,reported,duplicate_count)
        assert len(records) == int(summary[case]['taxa'])
        for sequence,taxa in by_sequence.items():
            if len(taxa)>1:
                groups.append(dict(case_id=case,taxa=';'.join(sorted(taxa)),taxon_count=len(taxa),aligned_sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),aligned_nucleotides=len(sequence)))
        rows.append(dict(case_id=case,taxa=len(records),distinct_aligned_sequences=len(by_sequence),redundant_aligned_sequences=duplicate_count,warning_blocks=len(blocks),warning_class='identical_aligned_sequence_advisory' if blocks else 'none',alignment_sha256=sha(alignment),fit_log_sha256=sha(folder/'fit.log')))
    root=Path('results/cds/local-codon-duplicate-warning-audit-20260927-v1');root.mkdir(parents=True,exist_ok=False)
    for name,data in [('cases.tsv',rows),('identical_sequence_groups.tsv',groups)]:
        with (root/name).open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');w.writeheader();w.writerows(data)
    result=dict(status='complete_full_local_duplicate_warning_audit',cases=len(rows),warning_cases=sum(r['warning_blocks']>0 for r in rows),warning_blocks=sum(r['warning_blocks'] for r in rows),redundant_aligned_sequences=sum(r['redundant_aligned_sequences'] for r in rows),identical_sequence_groups=len(groups),taxon_memberships_in_identical_groups=sum(r['taxon_count'] for r in groups),unclassified_warnings=0,source_fit_receipt_sha256=sha(fits/'receipt.json'),source_audit_receipt_sha256=sha(audit/'receipt.json'),script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in root.iterdir()},scope='Every fit warning matched to the identical-sequence advisory and its count independently reconstructed from checksum-bound aligned FASTA strings. Identical alignment strings can include missing characters and need not imply identical full proteins or genomes. Species are retained; this does not clear near-zero branches, independence concerns, optimization, or selection eligibility.')
    (root/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
