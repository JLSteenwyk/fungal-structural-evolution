"""Record terminal production and complete disposition coverage, not numeric validation."""
from collections import Counter
import csv
import json
from pathlib import Path
import subprocess
from screen_duplication_domain_alignment_coverage import sha

unit = 'fungal-duplication-alignments-20260926.service'
state = dict(line.split('=', 1) for line in subprocess.check_output(
    ['systemctl', '--user', 'show', unit, '-p', 'ActiveState', '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0')
planpath = Path('metadata/duplication_alignment_plan_20260926.json')
plan = json.loads(planpath.read_text())
root = Path(plan['output'])
receiptpath = root/'receipt.json'
receipt = json.loads(receiptpath.read_text())
assert receipt['status'] == 'complete_duplication_alignment_dispositions_pending_readback'
assert receipt['plan_sha256'] == sha(planpath)
manifest = root/'checkpoint_manifest.tsv'
assert sha(manifest) == receipt['artifacts'][manifest.name]
queue = Path(plan['queue'])
qr = json.loads((queue/'receipt.json').read_text())
assert sha(queue/'model_pairs.tsv') == qr['artifacts']['model_pairs.tsv']
with (queue/'model_pairs.tsv').open() as handle:
    pairs = [r['pair_key'] for r in csv.DictReader(handle, delimiter='\t')]
assert len(pairs) == len(set(pairs)) == receipt['distinct_model_pairs'] == 103200
expected = {(p, m, o) for p in pairs for m in ['full', 'plddt70'] for o in ['0', '1']}
seen, counts = set(), Counter()
with manifest.open() as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        p, m, o = Path(row['path']).stem.rsplit('-', 2)
        key = p, m, o
        assert key in expected and key not in seen
        assert row['path'] == f'pairs/{p[:2]}/{p}-{m}-{o}.json'
        assert len(row['sha256']) == 64 and all(c in '0123456789abcdef' for c in row['sha256'])
        seen.add(key)
        counts[m+':'+row['status']] += 1
assert seen == expected and len(seen) == receipt['directed_dispositions'] == 412800
assert dict(counts) == receipt['counts']
proof = dict(status='complete_primary_alignment_production_numeric_audit_pending',
             unit=unit, terminal_state=state, distinct_pairs=len(pairs), dispositions=len(seen), counts=dict(counts),
             receipt_sha256=sha(receiptpath), manifest_sha256=sha(manifest), script_sha256=sha(__file__),
             scope='Complete pair/mask/order coverage and receipt/manifest binding verified after terminal producer success. Individual checkpoint hashes, PDB mappings and numeric values remain the responsibility of the active full auditor; not certified here. No biological duplication effect claimed.')
Path('metadata/primary_alignment_production_closed_20260927.json').write_text(json.dumps(proof, indent=2)+'\n')
print(json.dumps(proof, indent=2))
