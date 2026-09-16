#!/usr/bin/env python3
"""Validate synthetic native discovery results independently of exit status."""
import argparse
import hashlib
import json
from pathlib import Path
from audit_orthology_family_universe import groups


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    receipt=json.loads((args.input/'receipt.json').read_text()); reports=[]
    for case in receipt['cases']:
        folder=args.input/case['case']; log=folder/'stdout.log'
        if sha(log)!=case['log_sha256'] or case['timed_out'] or not case['source_unchanged']:
            raise ValueError('Fixture did not finish cleanly')
        for name,h in case['source_sha256'].items():
            if sha(folder/'inputs'/name)!=h:raise ValueError('Source modified')
        clusters=list((folder/'native-results').glob('Results_*/WorkingDirectory/clusters*_id_pairs.txt'))
        trees=list((folder/'native-results').glob('Results_*/WorkingDirectory/Trees_ids/*.txt'))
        if case['taxa']==1:
            if clusters or 'At least two species are required' not in log.read_text():
                raise ValueError('Unexpected single-taxon behavior')
            reports.append(dict(case=case['case'],behavior='rejected_single_taxon_without_clusters',returncode=case['returncode']))
            continue
        if len(clusters)!=1:raise ValueError('Expected one completed cluster file')
        wd=clusters[0].parent
        ids=dict(line.strip().split(': ',1) for line in (wd/'SequenceIDs.txt').read_text().splitlines())
        observed=[frozenset(ids[g] for g in genes) for _,genes in groups(clusters[0])]
        expected={frozenset(f'{100+s}_{2*f+c}' for s in range(case['taxa']) for c in range(2)) for f in range(24)}
        if len(observed)!=24 or set(observed)!=expected:
            raise ValueError('Synthetic family partition differs')
        if case['no_fix_files'] and trees:raise ValueError('Groups-only run inferred gene trees')
        if not case['no_fix_files'] and len(trees)!=24:raise ValueError('Default downstream behavior differs')
        reports.append(dict(case=case['case'],behavior='exact_synthetic_partition',families=24,proteins=sum(map(len,observed)),gene_trees=len(trees),returncode=case['returncode'],cluster_sha256=sha(clusters[0])))
    output=dict(status='passed_native_discovery_software_contract_readback',cases=reports,
                receipt_sha256=sha(args.input/'receipt.json'),script_sha256=sha(Path(__file__)),
                conclusion='For >=2 taxa use --only-groups --no-fix-files to avoid unintended downstream inference. Single-taxon -f input is rejected with exit 0 and requires a separately validated path. Validate artifact contents rather than exit status.',
                limitations='Synthetic exact-copy grouping and software control flow only, not biological accuracy or scale validation.')
    args.output.write_text(json.dumps(output,indent=2)+'\n');print(json.dumps(output,indent=2))


if __name__=='__main__':main()
