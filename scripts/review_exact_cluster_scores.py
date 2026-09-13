#!/usr/bin/env python3
"""Review exact-score conversion against all original directed cluster alignments."""
import argparse,csv,json,math
from collections import Counter
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha,read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['validation','exact','output']:p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable review')
    old=checked_receipt(a.validation);completion=json.loads((a.exact/'completion.json').read_text());config=json.loads((a.exact/'config.json').read_text())
    if completion['config_sha256']!=sha(a.exact/'config.json') or completion['alignment_table_sha256']!=sha(a.exact/'alignments.tsv') or config['source_validation_receipt_sha256']!=sha(a.validation/'receipt.json'):raise ValueError('Exact conversion lineage mismatch')
    command=config['command']
    if command[command.index('--exact-tmscore')+1]!='1':raise ValueError('Exact-score option not selected')
    fields=command[command.index('--format-output')+1].split(',');before={(r['query'],r['target']):r for r in read_table(a.validation/'directed_alignment_review.tsv')};after={};counts=Counter()
    fixed=['evalue','qstart','qend','qlen','tstart','tend','tlen','alnlen','qcov','tcov']
    with (a.exact/'alignments.tsv').open() as f:
        for values in csv.reader(f,delimiter='\t'):
            r=dict(zip(fields,values));key=r['query'],r['target']
            if len(values)!=len(fields) or key not in before or key in after:raise ValueError('Unexpected/duplicate exact pair')
            if any(r[k]!=before[key][k] for k in fixed):raise ValueError('Underlying alignment fields changed')
            nums={k:float(r[k]) for k in fields[2:]}
            if any(not math.isfinite(v) or v<0 for v in nums.values()):raise ValueError('Nonfinite/negative exact output')
            invalid=any(nums[k]>1 for k in ['qcov','tcov','alntmscore','qtmscore','ttmscore','lddt'])
            passed=not invalid and nums['evalue']<=1e-5 and nums['qcov']>=.8 and nums['tcov']>=.8 and nums['alntmscore']>=.5
            prior=before[key]['passes_reported_cluster_criteria']=='True'
            r.update(exact_option_score_out_of_bounds=invalid,passes_exact_reported_criteria=passed,passes_approximate_reported_criteria=prior,decision_changed=passed!=prior)
            after[key]=r;counts['invalid_score_rows']+=invalid;counts['criterion_decision_changes']+=passed!=prior
    if set(after)!=set(before) or len(after)!=old['returned_directed_pairs']:raise ValueError('Exact pair grid incomplete')
    edges=[]
    for row in read_table(a.validation/'edge_review.tsv'):
        rep,member=row['representative_input_id'],row['member_input_id'];forward=after[rep,member];reverse=after[member,rep];passes=[x['passes_exact_reported_criteria'] for x in [forward,reverse]]
        status='both_directions_pass' if all(passes) else 'one_direction_passes' if any(passes) else 'neither_direction_passes'
        edges.append({'representative_input_id':rep,'member_input_id':member,'approximate_status':row['status'],'exact_status':status,'forward_exact_passes':passes[0],'reverse_exact_passes':passes[1],'any_exact_score_out_of_bounds':forward['exact_option_score_out_of_bounds'] or reverse['exact_option_score_out_of_bounds']})
    a.output.mkdir(parents=True);write_table(a.output/'directed_score_review.tsv',list(after.values()));write_table(a.output/'edge_score_review.tsv',edges)
    result={'status':'complete_exact_option_cluster_edge_review','validation_receipt_sha256':sha(a.validation/'receipt.json'),'exact_completion_sha256':sha(a.exact/'completion.json'),'script_sha256':sha(Path(__file__)),'directed_pairs':len(after),'edge_pairs':len(edges),'checks':dict(counts),'exact_edge_status_counts':dict(Counter(r['exact_status'] for r in edges)),'edge_decision_changes':sum(r['approximate_status']!=r['exact_status'] for r in edges),'interpretation':'Exact-TM option of the same Foldseek implementation on unchanged saved alignments. Identity, fixed alignment fields, ranges and threshold decisions checked. Out-of-range scores are explicit failures of this review. Finite printed precision and fixed-alignment dependence remain. Original memberships remain immutable; no all-pairs, orthology, novelty or experimental validity claim.','artifacts':{p.name:sha(p) for p in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
