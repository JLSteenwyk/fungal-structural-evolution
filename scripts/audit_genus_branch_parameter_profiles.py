#!/usr/bin/env python3
"""Audit every saved MG94 profile point and summarize finite-grid diagnostics."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def table(path):
    with path.open() as handle:return list(csv.DictReader(handle,delimiter='\t'))


def declarations(text):
    values={};fixed=set()
    for line in text.splitlines():
        if '=' not in line:continue
        left,right=line.split('=',1);name=left.removeprefix('global ').rstrip(':')
        if not (line.startswith('global ') and not left.endswith(':') and '.model_MGREV.' in name or '.tree_0.' in name and name.endswith('.t')):continue
        try:value=float(right.split(';',1)[0])
        except ValueError:continue
        if name in values:raise ValueError('Duplicate fitted declaration')
        values[name]=value
        if left.endswith(':'):fixed.add(name)
    return values,fixed


def unchanged_model_text(text,parameters):
    lines=[]
    for line in text.splitlines():
        left,sep,right=line.partition('=')
        name=left.removeprefix('global ').rstrip(':')
        if sep and name in parameters:
            line='PROFILE_PARAMETER '+name+';'+right.partition(';')[2]
        lines.append(line)
    return '\n'.join(lines)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['profiles','fits','slices','plan','output']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    root=json.loads((a.profiles/'receipt.json').read_text());config=json.loads((a.profiles/'config.json').read_text());plan=json.loads(a.plan.read_text())
    if root['status']!='complete_full_case_branch_parameter_profiles_pending_summary_audit' or sha(a.profiles/'config.json')!=root['config_sha256'] or sha(a.plan)!=config['plan_sha256']:raise ValueError('Incomplete or changed profile batch')
    for folder,key in [(a.fits,'source_fit_receipt_sha256'),(a.slices,'source_slice_receipt_sha256')]:
        if sha(folder/'receipt.json')!=config[key]:raise ValueError('Source receipt differs')
    for name,digest in root['artifacts'].items():
        if sha(a.profiles/name)!=digest:raise ValueError('Profile table changed')
    sr=json.loads((a.slices/'receipt.json').read_text())
    for name,digest in sr['artifacts'].items():
        if sha(a.slices/name)!=digest:raise ValueError('Slice artifact changed')
    original_cases={r['case_id']:r['receipt_sha256'] for r in json.loads((a.fits/'receipt.json').read_text())['case_receipts']}
    exposure=Path('results/qc/fcs-marker-analysis-exposure-v1')
    er=json.loads((exposure/'receipt.json').read_text())
    if er['source_fit_receipt_sha256']!=sha(a.fits/'receipt.json') or sha(exposure/'codon_case_exposure.tsv')!=er['artifacts']['codon_case_exposure.tsv']:raise ValueError('FCS exposure source differs')
    fcs={}
    for row in table(exposure/'codon_case_exposure.tsv'):
        if row['fit_receipt_sha256']!=original_cases[row['case_id']]:raise ValueError('FCS case receipt differs')
        fcs.setdefault(row['case_id'],[]).append(row['taxon_id'])
    reviews={r['case_id']:r for r in table(a.slices/'case_review.tsv')}
    slices={(r['case_id'],float(r['branch_multiplier'])):float(r['log_likelihood']) for r in table(a.slices/'likelihood_slices.tsv')}
    points=table(a.profiles/'profile_points.tsv')
    point_key=lambda r:(r['case_id'],r['fit_kind'],str(r['branch_parameter_multiplier']))
    grid={point_key(r):r for r in points}
    if len(grid)!=len(points) or root['cases']!=plan['cases'] or len(root['case_receipts'])!=len(original_cases):raise ValueError('Incomplete case/point grid')
    summaries=[];seen=set();hashes=0;parameters_checked=0;max_readback=0.0;all_rows=0
    for proof in root['case_receipts']:
        case=proof['case_id'];folder=a.profiles/case;source=a.fits/case
        if case in seen or case not in reviews:raise ValueError('Invalid case identity')
        seen.add(case)
        if sha(folder/'receipt.json')!=proof['receipt_sha256'] or sha(source/'receipt.json')!=original_cases[case]:raise ValueError('Case receipt changed')
        r=json.loads((folder/'receipt.json').read_text());baseline=json.loads((source/'receipt.json').read_text())
        if r['status']!='complete_branch_parameter_profile_with_fresh_fit_readbacks' or r['source_fit_receipt_sha256']!=original_cases[case]:raise ValueError('Case source mismatch')
        for name,digest in baseline['artifacts'].items():
            if sha(source/name)!=digest:raise ValueError('Original fit artifact changed')
        original_text=(source/'fit.bf').read_text();original,_=declarations(original_text);target_names=[n for n in original if n.endswith('.tree_0.'+reviews[case]['target_node']+'.t')]
        if len([n for n in original if '.model_MGREV.' in n])!=6:raise ValueError('Expected five free exchangeabilities and omega')
        expected_nodes=set(json.loads((source/'fit.json').read_text())['branch attributes']['0'])
        if {n.split('.tree_0.')[1][:-2] for n in original if '.tree_0.' in n}!=expected_nodes or len(target_names)!=1:raise ValueError('Incomplete original parameter grid')
        target=target_names[0];rows=r['rows'];factors={float(x['branch_parameter_multiplier']) for x in rows if x['fit_kind']=='fixed_branch_parameter_profile'}
        if factors!=set(plan['branch_parameter_multipliers']) or len(rows)!=8 or len(r['proofs'])!=8:raise ValueError('Incomplete multiplier grid')
        bylabel={x['fit_kind']:x for x in r['proofs']}
        best=max(x['log_likelihood'] for x in rows)
        for row in rows:
            if point_key(row) not in grid or any(str(v)!=grid[point_key(row)][k] for k,v in row.items()):raise ValueError('Aggregate point table differs')
            factor=None if row['fit_kind']=='unconstrained_reoptimization' else row['branch_parameter_multiplier']
            label='unconstrained' if factor is None else 'factor_'+str(factor).replace('.','p');work=folder/label
            for name,digest in bylabel[label]['artifacts'].items():
                if sha(work/name)!=digest:raise ValueError('Profile point artifact changed')
                hashes+=1
            log=(work/'optimize.log').read_text().splitlines();opt=[x.split('\t') for x in log if x.startswith('OPTIMUM\t')];pars=[x.split('\t')[1:] for x in log if x.startswith('PARAM\t')]
            values={n:float(v) for n,v in pars};export_text=(work/'profile_fit.bf').read_text();exported,constrained=declarations(export_text)
            if unchanged_model_text(original_text,original)!=unchanged_model_text(export_text,original):raise ValueError('Non-profiled model content changed')
            replay=[float(x.partition('=')[2]) for x in (work/'readback.log').read_text().splitlines() if x.startswith('READBACK=')]
            if len(opt)!=1 or len(replay)!=1 or len(values)!=len(pars) or set(values)!=set(original) or values!=exported:raise ValueError('Incomplete or mismatched saved fitted parameters')
            if any(not math.isfinite(v) or v<0 for v in values.values()):raise ValueError('Invalid parameter')
            ll=row['log_likelihood'];start=baseline['log_likelihood'] if factor is None else slices[case,float(factor)]
            if not math.isfinite(ll) or max(abs(float(opt[0][1])-ll),abs(float(opt[0][2])-ll),abs(replay[0]-ll))>1e-6:raise ValueError('Likelihood readback mismatch')
            if constrained!=(set() if factor is None else {target}) or factor is not None and not math.isclose(values[target],original[target]*factor,rel_tol=1e-10,abs_tol=1e-10):raise ValueError('Fixed-target constraint differs')
            if ll<start-1e-5 or abs(row['improvement_over_fixed_nuisance_slice']-(ll-start))>1e-8 or abs(row['delta_ll_from_best_evaluated_fit']-(ll-best))>1e-8:raise ValueError('Likelihood diagnostic differs')
            max_readback=max(max_readback,abs(replay[0]-ll));parameters_checked+=len(values);all_rows+=1
        unconstrained=next(x for x in rows if x['fit_kind']=='unconstrained_reoptimization')
        fixed={float(x['branch_parameter_multiplier']):x for x in rows if x['fit_kind']=='fixed_branch_parameter_profile'}
        summaries.append({'case_id':case,'target_node':reviews[case]['target_node'],'original_dS_review_label':reviews[case]['target_dS'],
                          'original_t':original[target],'unconstrained_t':unconstrained['target_parameter'],
                          'unconstrained_ll_gain':unconstrained['delta_ll_from_original_fit'],
                          'best_grid_exceeds_unconstrained_by':max(x['log_likelihood'] for x in fixed.values())-unconstrained['log_likelihood'],
                          'factor_0p1_ll_deficit_from_best':best-fixed[.1]['log_likelihood'],
                          'factor_10_ll_deficit_from_best':best-fixed[10]['log_likelihood'],
                          'maximum_nuisance_reoptimization_gain':max(x['improvement_over_fixed_nuisance_slice'] for x in fixed.values()),
                          'marker_copy_caveat':reviews[case]['marker_copy_caveat'],
                          'fcs_exposed_taxa':';'.join(sorted(fcs.get(case,[]))),
                          'fcs_exposed_case':case in fcs})
    if seen!=set(original_cases) or all_rows!=root['optimized_fits'] or all_rows!=plan['total_optimized_fits']:raise ValueError('Incomplete final audit grid')
    a.output.mkdir(parents=True)
    with (a.output/'case_summary.tsv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,list(summaries[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(summaries)
    result={'status':'passed_full_branch_parameter_profile_artifact_and_numeric_audit','cases':len(seen),'optimized_points':all_rows,'saved_point_artifact_hashes_checked':hashes,'saved_parameter_values_checked':parameters_checked,'maximum_fresh_readback_likelihood_error':max_readback,
            'cases_grid_beats_unconstrained_by_gt_1e_5':sum(x['best_grid_exceeds_unconstrained_by']>1e-5 for x in summaries),
            'source_receipt_sha256':sha(a.profiles/'receipt.json'),'script_sha256':sha(Path(__file__)),
            'fcs_exposure_receipt_sha256':sha(exposure/'receipt.json'),'fcs_exposed_cases':len(fcs),
            'artifacts':{'case_summary.tsv':sha(a.output/'case_summary.tsv')},
            'interpretation':'Full finite-grid artifact and numerical readback, not a fresh rerun of optimizations. Profiles vary branch t while nuisance exchangeabilities and omega change; original dS is only a review label. Single-start finite grids can miss optima. Grid beating unconstrained fit identifies residual optimization concern. No calibrated intervals, absence-of-saturation certificate or selection conclusion. FCS and biological eligibility review remain separate.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
