#!/usr/bin/env python3
"""Summarize audited genus-tree warnings without treating trees as selection-ready."""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
from audit_genus_codon_trees import sha, table


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('Use a new summary directory')
    receipt=json.loads((args.audit/'receipt.json').read_text())
    if receipt['status'] not in ['passed_completed_case_snapshot','passed_full_genus_tree_audit']:
        raise ValueError('Passed source audit required')
    for name,digest in receipt['artifacts'].items():
        if sha(args.audit/name)!=digest:raise ValueError('Changed source artifact')
    cases=table(args.audit/'case_audit.tsv');branches=table(args.audit/'branch_support.tsv')
    counts=defaultdict(Counter)
    for row in table(args.audit/'warnings.tsv'):
        message=row['warning'].lower()
        if 'too low -mem' in message:category='memory_adjustment'
        elif 'saturated' in message:category='saturated_pairwise_distance'
        elif 'boundar' in message:category='parameter_boundary'
        elif 'nni search needs unusual' in message:category='nni_convergence'
        else:category='other'
        counts[row['case_id']][category]+=1
    categories=['memory_adjustment','saturated_pairwise_distance','parameter_boundary','nni_convergence','other']
    output=[]
    for row in cases:
        count=counts[row['case_id']]
        if sum(count.values())!=int(row['warning_lines']):raise ValueError('Warning tally differs')
        output.append(dict(row,**{category+'_warning_lines':count[category] for category in categories},
            review_required_before_codon_interpretation=bool(count['saturated_pairwise_distance'] or count['parameter_boundary'] or count['nni_convergence'] or count['other']),
            selection_eligibility='not_established'))
    result={'status':'complete_descriptive_audit_summary','source_audit_sha256':sha(args.audit/'receipt.json'),
        'source_audit_status':receipt['status'],'cases':len(cases),'bootstrap_trees':sum(int(r['bootstrap_trees']) for r in cases),
        'internal_edges':len(branches),'cases_by_code':dict(Counter(r['translation_table'] for r in cases)),
        'cases_with_warning_category':{c:sum(counts[r['case_id']][c]>0 for r in cases) for c in categories},
        'branches_ufb_lt95':sum(float(r['reported_ufb_percent'])<95 for r in branches),
        'branches_sh_alrt_lt80':sum(float(r['sh_alrt_percent'])<80 for r in branches),
        'cases_with_near_zero_edges':sum(int(r['near_zero_edges_le1e_5'])>0 for r in cases),
        'script_sha256':sha(Path(__file__)),
        'interpretation':'Descriptive completed-case snapshot, potentially biased by finishing time. Warning flags require follow-up, not automatic exclusion. IQ-TREE saturated-distance warnings are not a codon-specific synonymous-saturation test. Memory-adjustment warnings are retained separately; their cause is unresolved. Selection eligibility is not established even for unflagged cases.'}
    args.output.mkdir(parents=True)
    with (args.output/'case_review.tsv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,list(output[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(output)
    result['artifacts']={'case_review.tsv':sha(args.output/'case_review.tsv')}
    (args.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
