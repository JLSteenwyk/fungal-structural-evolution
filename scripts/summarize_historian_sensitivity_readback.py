#!/usr/bin/env python3
"""Summarize an immutable Historian audit without pooling distinct sensitivities."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    receipt_path = args.audit/'receipt.json'
    receipt = json.loads(receipt_path.read_text())
    for name, digest in receipt['artifacts'].items():
        assert sha(args.audit/name) == digest
    dispositions = json.loads((args.audit/'dispositions.json').read_text())
    assert len(dispositions) == receipt['expected_jobs'] == 324
    assert dict(Counter(r['status'] for r in dispositions)) == receipt['dispositions']
    by_id = {r['job_id']: r for r in dispositions}
    assert len(by_id) == 324
    with (args.audit/'paired_sequences.tsv').open() as handle:
        pairs = list(csv.DictReader(handle, delimiter='\t'))
    assert len(pairs) == receipt['paired_node_comparisons']
    groups = defaultdict(list)
    for row in pairs:
        for key in ['job_a', 'job_b']:
            assert by_id[row[key]]['status'] == 'independent_tip_tree_and_candidate_readback_passed'
        assert (row['same_sequence']=='True') == (int(row['edit_distance'])==0)
        assert (row['root_assumption']=='True') == (int(row['level'])==3)
        family = row['input_id'].split('-')[0]
        dataset = 'whole' if '-whole-' in row['input_id'] else 'domain'
        groups[family, dataset, row['comparison'], row['root_assumption']].append(row)
    summaries = []
    for (family,dataset,comparison,assumed_root), rows in sorted(groups.items()):
        summaries.append(dict(family=family,dataset=dataset,comparison=comparison,
                              assumed_root=assumed_root,dependent_node_comparisons=len(rows),
                              changed=sum(int(r['edit_distance'])>0 for r in rows),
                              maximum_edit_distance=max(int(r['edit_distance']) for r in rows)))
    # This polytomy has six alignment inputs, two floors and all 15 resolutions.
    star = [r for r in dispositions if r['job_id'].startswith('OG0000230-')]
    assert len(star) == 180 and all(r['status']=='independent_tip_tree_and_candidate_readback_passed' for r in star)
    star_groups = defaultdict(set)
    for row in star:
        name = row['job_id']
        prefix, resolution = name.rsplit('-resolution',1)
        star_groups[prefix].add(int(resolution))
    assert len(star_groups)==12 and all(v==set(range(15)) for v in star_groups.values())
    star_pairs = [r for r in pairs if r['input_id'].startswith('OG0000230-')]
    star_types = Counter(r['comparison'] for r in star_pairs)
    assert star_types == {'star_resolution':5040,'branch_floor':360}
    args.output.mkdir(parents=True,exist_ok=False)
    with (args.output/'family_sensitivity.tsv').open('w') as handle:
        writer=csv.DictWriter(handle,list(summaries[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(summaries)
    result=dict(status='completed_summary_of_frozen_partial_audit' if receipt['dispositions'].get('pending_receipt') else 'completed_summary_of_full_audit',
                input_receipt=str(receipt_path),input_receipt_sha256=sha(receipt_path),
                script_sha256=sha(__file__),dispositions=receipt['dispositions'],
                dependent_node_comparisons=len(pairs),changed_node_comparisons=sum(r['changed'] for r in summaries),
                complete_polytomy_family=dict(family='OG0000230',audited_jobs=180,resolutions_per_input_and_floor=15,
                    comparisons_by_factor=dict(star_types),changed=sum(int(r['edit_distance'])>0 for r in star_pairs)),
                artifacts={'family_sensitivity.tsv':sha(args.output/'family_sensitivity.tsv')},
                scope='Original capacity-grid audit only. Higher-memory retries and cross-alignment comparisons are separate and must be reported alongside this summary. Zero edits in these fixed-parameter profile diagnostics do not imply posterior concentration, biological invariance, convergence, or a uniquely identified root. Pairwise alternatives share data and are not independent replicates.')
    (args.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':
    main()
