#!/usr/bin/env python3
"""Inventory the full expanded fixed-effect design grid without running fits."""
import argparse
from collections import Counter, defaultdict
import csv
import fcntl
import gzip
import itertools
import json
from pathlib import Path
import shutil
import numpy as np
from full_expanded_model_design_sources import load, array_digest, cohort_id, design_id, fit_id, AXES, OUTCOMES, DEGREES, SETTING_FIELDS, SUMMARY_FIELDS
from full_expanded_model_input_sources import ORDERS, NUISANCE, STRATUM
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def rank_audit(matrix, degree):
    n, p = matrix.shape
    assert np.isfinite(matrix).all()
    maxabs = np.max(np.abs(matrix), axis=0, initial=0)
    safe = np.where(maxabs > 0, maxabs, 1.)
    rescaled = matrix / safe
    lengths = np.sqrt(np.sum(rescaled**2, axis=0))
    active = list(range(p)) if not n else np.flatnonzero(maxabs > 0).tolist()
    normalized = (rescaled / np.where(lengths > 0, lengths, 1.))[:,active]
    singular = np.linalg.svd(normalized, full_matrices=False, compute_uv=False)
    tolerance = max(n, p) * np.finfo(float).eps * (float(singular[0]) if len(singular) else 0.)
    ranks = [int((singular > tolerance * f).sum()) for f in [.1, 1., 10.]]
    q = len(active)
    status = 'empty_setting' if not n else 'sequence_axis_uninformative_requires_review' if not any(k in active for k in range(1,degree+1)) else 'insufficient_residual_dimension' if n <= q else 'rank_boundary_requires_review' if ranks[0] != ranks[2] else 'rank_deficient_requires_review' if ranks[1] < q else 'full_rank_design'
    return dict(records=n, columns=p, active_columns=q, active_column_indices=active, exactly_zero_column_indices=[k for k in range(p) if k not in active], column_maxabs=maxabs.tolist(), column_l2_after_maxabs=lengths.tolist(),
        singular_values=singular.tolist(), rank_tolerance=tolerance, rank_at_one_tenth_tolerance=ranks[0],
        rank=ranks[1], rank_at_ten_times_tolerance=ranks[2],
        normalized_condition_number=float(singular[0]/singular[-1]) if len(singular) == q and ranks[1] == q and len(singular) else None,
        scaled_design_sha256=array_digest(normalized, '<f8'), disposition=status)


def memberships(source):
    lookup = {c['case_id']: i for i, c in enumerate(source['cases'])}; groups = defaultdict(list); reuse = Counter()
    with gzip.open(source['selections'], 'rt') as f:
        for number, row in enumerate(csv.DictReader(f, delimiter='\t'), 1):
            assert int(row['source_row_ordinal']) == number
            i = lookup[row['case_id']]; case = source['cases'][i]
            for k in ['guide', 'physical_case_id', 'target_id', 'background_id']: assert row[k] == case[k]
            groups[row['guide'], row['policy'], row['scenario_id']].append(i); reuse[row['case_id']] += 1
    assert number == source['config']['expected']['selected_records']
    assert reuse == {c['case_id']: int(c['selection_records']) for c in source['cases']}
    return groups


def run(path, stop_after_cohorts=None):
    path = Path(path); plan = json.loads(path.read_text()); source, bindings = load(plan, path)
    root = Path(plan['output']); assert shutil.disk_usage(root.parent).free >= plan['resources']['minimum_free_disk_gib']*2**30
    root.mkdir(exist_ok=True); lock = (root/'stage.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
    assert not (root/'receipt.json').exists(), 'Completed design inventory cannot restart'
    state = dict(plan_sha256=sha(path), schema='expanded-model-design-v1'); marker = root/'stage_plan.json'
    if marker.exists(): assert json.loads(marker.read_text()) == state
    else: marker.write_text(json.dumps(state, indent=2)+'\n')
    groups = memberships(source); cases = source['cases']; cohorts = {}; mapping = {}
    bits = {m: np.asarray([int(c['joint_'+m+'_pass_bits']) for c in cases],dtype=np.int64) for m in ['full','plddt70','both']}
    screen_bits = {s['id']: i for i,s in enumerate(source['config']['screens'])}
    for row in source['counts']:
        key = tuple(row[k] for k in STRATUM[:-1])
        if key in mapping: continue
        guide, policy, scenario, mask, gate, screen = key
        candidates = groups.get((guide,policy,scenario), [])
        assert len({cases[i]['target_id'] for i in candidates}) == len(candidates)
        b = bits[mask if gate == 'mask' else 'both']; flag = 1 << screen_bits[screen]
        members = np.asarray(sorted([i for i in candidates if b[i] & flag], key=lambda i: cases[i]['case_id']),dtype=np.int64)
        assert len(members) == int(row['quality_retained_numeric_records'])
        ids_sha = array_digest([cases[i]['case_id'] for i in members], 'S64'); identifier = cohort_id(guide,mask,ids_sha)
        if identifier not in cohorts: cohorts[identifier] = dict(cohort_id=identifier, guide=guide, mask=mask, members=members, ordered_case_ids_sha256=ids_sha, membership_occurrences=0)
        else: assert np.array_equal(cohorts[identifier]['members'],members)
        cohorts[identifier]['membership_occurrences'] += 1; mapping[key] = identifier
    folder = root/'cohorts'; folder.mkdir(exist_ok=True); manifest=[]; model_map={}; ds=Counter(); fs=Counter(); ss=Counter()
    design_count=fit_count=0
    with (root/'unique_designs.jsonl').open('w') as df, (root/'unique_fit_inputs.jsonl').open('w') as ff:
        for number,(cid,cohort) in enumerate(sorted(cohorts.items()),1):
            members=cohort['members']; fp=folder/(cid+'.npz'); np.savez_compressed(fp,case_rows=members)
            info={k:v for k,v in cohort.items() if k!='members'}
            info.update(path=str(fp.relative_to(root)),sha256=sha(fp),records=len(members),case_rows_sha256=array_digest(members,'<i8'))
            manifest.append(info)
            dep = dict(species_pattern_rows_sha256=array_digest([int(source['cov'][cases[i]['case_id']]['species_pattern_row']) for i in members],'<i8'),
                family_components_sha256=array_digest([source['cov'][cases[i]['case_id']]['family_component'] for i in members],'S64'),
                background_nodes_sha256=array_digest([cases[i]['background_id'] for i in members],'S64'),
                background_physical_pairs_sha256=array_digest([cases[i]['background_pair_key'] for i in members],'S64'),
                unique_family_components=len({source['cov'][cases[i]['case_id']]['family_component'] for i in members}),
                unique_species_patterns=len({source['cov'][cases[i]['case_id']]['species_pattern_id'] for i in members}),
                largest_family_component=max(Counter(source['cov'][cases[i]['case_id']]['family_component'] for i in members).values(),default=0))
            for order in ORDERS:
                values=source['arrays'][cohort['mask'],order]; assert np.all(values['numerical_usable'][members]==1)
                ordered_ids=array_digest(values['input_id'][members],'S64')
                for axis,degree in itertools.product(AXES,DEGREES):
                    columns=['intercept',*AXES[axis][:degree],*NUISANCE]
                    matrix=np.column_stack([np.ones(len(members))]+[values[k][members] for k in columns[1:]])
                    audit=rank_audit(matrix,degree); identifier=design_id(cid,order,axis,degree,source['contract'])
                    design=dict(design_id=identifier,cohort_id=cid,mask=cohort['mask'],order_contrast=order,sequence_axis=axis,degree=degree,
                        predictor_columns=columns,ordered_input_ids_sha256=ordered_ids,raw_design_sha256=array_digest(matrix,'<f8'),**dep,**audit)
                    df.write(json.dumps(design,sort_keys=True,allow_nan=False)+'\n');design_count+=1;ds[design['disposition']]+=1
                    for outcome in OUTCOMES:
                        y=values[outcome][members];assert np.isfinite(y).all()
                        response=array_digest(y,'<f8');fid=fit_id(identifier,outcome,response)
                        status=design['disposition'] if design['disposition']!='full_rank_design' else 'constant_response_requires_review' if np.ptp(y)==0 else 'ready_for_working_covariance_fit'
                        fit=dict(fit_input_id=fid,design_id=identifier,cohort_id=cid,outcome=outcome,response_sha256=response,
                            records=len(y),response_min=float(y.min()) if len(y) else None,response_max=float(y.max()) if len(y) else None,
                            disposition=status,trees=plan['trees'],source_contract=source['contract'])
                        ff.write(json.dumps(fit,sort_keys=True,allow_nan=False)+'\n');fit_count+=1;fs[status]+=1
                        model_map[cid,order,axis,degree,outcome]=(identifier,fid,status)
            if number%25==0: print('full_expanded_design_cohorts',number,'/',len(cohorts),flush=True)
            if stop_after_cohorts==number: raise InterruptedError('Software interruption contract')
    (root/'cohort_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    settings=0
    with gzip.open(root/'model_settings.tsv.gz','wt') as f:
        w=csv.DictWriter(f,SETTING_FIELDS,delimiter='\t',lineterminator='\n');w.writeheader()
        for row in source['counts']:
            cid=mapping[tuple(row[k] for k in STRATUM[:-1])]
            for outcome,axis,degree in itertools.product(OUTCOMES,AXES,DEGREES):
                identifier,fid,status=model_map[cid,row['order_contrast'],axis,degree,outcome]
                w.writerow({**{k:row[k] for k in STRATUM},'outcome':outcome,'sequence_axis':axis,'degree':degree,
                    'cohort_id':cid,'design_id':identifier,'fit_input_id':fid,'records':len(cohorts[cid]['members']),
                    'disposition':status,'nominal_tree_fits':len(plan['trees'])})
                settings+=1;ss[status]+=1
    assert settings==plan['expected']['model_setting_rows']==len(source['counts'])*12
    assert design_count==len(cohorts)*30 and fit_count==design_count*2
    summary=dict(logical_cases=len(cases),selected_records=source['config']['expected']['selected_records'],input_setting_rows=len(source['counts']),
        model_setting_rows=settings,unique_cohorts=len(cohorts),unique_designs=design_count,unique_fit_inputs=fit_count,
        nominal_tree_setting_fits=settings*len(plan['trees']),unique_tree_fit_inputs=fit_count*len(plan['trees']),
        design_status_counts=dict(ds),fit_input_status_counts=dict(fs),setting_status_counts=dict(ss),
        largest_cohort=max(c['records'] for c in manifest),cohort_member_occurrences=sum(c['records'] for c in manifest),trees=plan['trees'])
    verify(bindings)
    names=['stage_plan.json','cohort_manifest.json','unique_designs.jsonl','unique_fit_inputs.jsonl','model_settings.tsv.gz']+[r['path'] for r in manifest]
    receipt=dict(status='complete_full_expanded_model_designs_pending_independent_readback',plan_sha256=sha(path),**summary,
        artifacts={n:sha(root/n) for n in names},source_hashes=bindings,scientific_eligibility=False,scope=plan['scope'])
    with (root/'receipt.json').open('x') as f:f.write(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(summary),flush=True);return receipt


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);run(p.parse_args().plan)
