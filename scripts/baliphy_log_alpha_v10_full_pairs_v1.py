"""Inventory and compare every original scientific file for the complete V10 grid."""
from pathlib import Path

from ancestral_chain_attempt import sha


FILES = ['C1.log', 'C1.log.json', 'C1.log.column-map.json', 'runtime-tree.nwk',
         'C1.P1.fastas', 'C1.P1.site-property-samples.jsonl']


def same_bytes(first, second):
    with first.open('rb') as a, second.open('rb') as b:
        while True:
            left, right = a.read(65536), b.read(65536)
            if left != right:
                return False
            if not left:
                return True


def inspect_pairs(jobs, rows):
    by_id = {r['chain_id']: r for r in rows}
    assert len(jobs) == len(rows) == len(by_id) == 24
    assert set(by_id) == {j['chain']['chain_id'] for j in jobs}
    results = []
    for job in jobs:
        cid = job['chain']['chain_id']
        row = by_id[cid]
        old_receipt = Path(job['paired_v7_native_receipt'])
        assert sha(old_receipt) == job['paired_v7_native_receipt_sha256']
        old = old_receipt.parent / 'independent-chain-1'
        current = Path(row['native_receipt']).parent / 'independent-chain-1'
        files = []
        for name in FILES:
            first, second = old / name, current / name
            assert first.is_file(), first
            files.append(dict(name=name, original_path=str(first), current_path=str(second),
                original_sha256=sha(first), current_sha256=sha(second) if second.is_file() else None,
                current_present=second.is_file(), byte_identical=second.is_file() and same_bytes(first, second)))
        diagnostic = current / 'C1.P1.log-alpha-samples.jsonl'
        count = 0
        if diagnostic.is_file():
            with diagnostic.open('rb') as handle:
                count = max(0, sum(1 for _ in handle) - 1)
        identical = row['exit_code'] == 0 and all(f['byte_identical'] for f in files)
        results.append(dict(chain_id=cid, source_v7_chain_id=job['source_v7_chain_id'],
            native_exit_code=row['exit_code'], files=files,
            all_original_files_byte_identical=identical,
            diagnostic_path=str(diagnostic), diagnostic_present=diagnostic.is_file(),
            diagnostic_sha256=sha(diagnostic) if diagnostic.is_file() else None,
            diagnostic_rows_present=count,
            status='all_original_scientific_files_byte_identical_pending_diagnostic_readback' if identical else
                   'native_or_paired_scientific_output_review_required',
            scientific_eligibility=False, posterior_qualified=False))
    return sorted(results, key=lambda r: r['chain_id'])
