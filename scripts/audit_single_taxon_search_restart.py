#!/usr/bin/env python3
"""Check exact synthetic single-taxon family membership after native restart."""
import argparse
import hashlib
import json
from pathlib import Path
from audit_orthology_family_universe import groups


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    root=args.input;r=json.loads((root/'receipt.json').read_text())
    if r['timed_out'] or not r['original_inputs_unchanged'] or sha(root/'stdout.log')!=r['log_sha256']:
        raise ValueError('Incomplete fixture')
    clusters=list(root.glob('WorkingDirectory/OrthoFinder/Results_*/WorkingDirectory/clusters*_id_pairs.txt'))
    if len(clusters)!=1:raise ValueError('Expected one native output cluster file')
    wd=clusters[0].parent
    ids=dict(line.split(': ',1) for line in (root/'WorkingDirectory/SequenceIDs.txt').read_text().splitlines())
    observed=[frozenset(ids[x] for x in genes) for _,genes in groups(clusters[0])]
    expected={frozenset([f'100_{2*f}',f'100_{2*f+1}']) for f in range(24)}
    if len(observed)!=24 or set(observed)!=expected:raise ValueError('Synthetic partition differs')
    if list(wd.glob('Trees_ids/*.txt')):raise ValueError('Unexpected tree inference')
    for name,h in r['original_input_sha256'].items():
        if sha(root/'WorkingDirectory'/name)!=h:raise ValueError('Source changed')
    result=dict(status='passed_single_taxon_native_search_restart_readback',taxa=1,proteins=48,families=24,
                exact_synthetic_family_membership=True,gene_trees=0,original_inputs_unchanged=True,
                receipt_sha256=sha(root/'receipt.json'),cluster_sha256=sha(clusters[0]),script_sha256=sha(Path(__file__)),
                scope='All 48 original protein IDs and all 24 synthetic duplicate-pair families checked through native -b self-search restart. No fictitious taxa or biological subset. Production self-search generation and full-scale normalization remain to be validated.')
    args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
