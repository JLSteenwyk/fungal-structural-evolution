#!/usr/bin/env python3
"""Summarize every paired predictor draw without accepting biological intervals."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
from scipy.stats import beta

from ancestral_chain_attempt import sha
from matched_predictor_branch_fits import ROLES
from matched_predictor_resampling import read_arrays
from readback_matched_predictor_resampling import independent_draw
from reference_measurement_union_sources import bind, verify

MODELS = ['af', 'af_empirical', 'llm']


def statistics(aa, af, esm, point_af, point_esm):
    assert aa.shape == af.shape == esm.shape == (200,)
    assert all(np.isfinite(x).all() and np.all(x >= 0) for x in [aa, af, esm])
    delta = esm - af
    assert np.isfinite(delta).all()
    positive, negative, zero = [int(np.sum(op(delta, 0))) for op in [np.greater, np.less, np.equal]]
    assert positive + negative + zero == 200
    q = np.quantile(delta, [.025, .5, .975], method='linear')
    lo = 0.0 if positive == 0 else float(beta.ppf(.025, positive, 201-positive))
    hi = 1.0 if positive == 200 else float(beta.ppf(.975, positive+1, 200-positive))
    result = dict(point_AlphaFold=float(point_af), point_ESMFold=float(point_esm),
        point_ESMFold_minus_AlphaFold=float(point_esm-point_af),
        delta_mean=float(np.mean(delta)), delta_sd=float(np.std(delta, ddof=1)),
        delta_q025=float(q[0]), delta_median=float(q[1]), delta_q975=float(q[2]),
        positive_draws=positive, negative_draws=negative, zero_draws=zero,
        positive_fraction=positive/200, positive_fraction_mc95_low=lo,
        positive_fraction_mc95_high=hi,
        empirical_central_range_excludes_zero=bool(q[0] > 0 or q[2] < 0),
        empirical_direction='positive' if q[0] > 0 else 'negative' if q[2] < 0 else 'includes_zero')
    for label, vector in [('aa', aa), ('AlphaFold', af), ('ESMFold', esm)]:
        values = np.quantile(vector, [.025, .5, .975], method='linear')
        result.update({label+'_q025':float(values[0]), label+'_median':float(values[1]),
                       label+'_q975':float(values[2]), label+'_zero_draws':int(np.sum(vector==0))})
    assert all(np.isfinite(v) for v in result.values() if isinstance(v, float))
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True); p.add_argument('--receipt', type=Path, required=True)
    a = p.parse_args(); assert not a.receipt.exists() and not a.output.exists()
    gate_path = Path('metadata/matched_predictor_resampling_full_closure_verified_20261004_v1.json')
    gate = json.loads(gate_path.read_text()); verify(gate['source_hashes'])
    assert gate['paired_uncertainty_stage_numerically_closed'] and gate['unresolved_cases'] == 0
    cp = Path('metadata/matched_predictor_resampling_completed_20261003_v1.json')
    completed = json.loads(cp.read_text()); archive_path = Path(completed['full_hash_archive'])
    assert sha(archive_path) == completed['full_hash_archive_sha256']
    archive = json.loads(archive_path.read_text()); expected_pins = archive['source_hashes']
    plan_path = Path('metadata/matched_predictor_resampling_plan_20261003_v1.json')
    plan = json.loads(plan_path.read_text()); root = Path(plan['output'])
    inputs_plan_path = Path(plan['input_plan']); inputs_plan = json.loads(inputs_plan_path.read_text())
    input_root = Path(inputs_plan['output']); point_root = Path(plan['point_output'])
    axes_path = root/'axes.json'; manifest_path = root/'case_manifest.json'
    axes = json.loads(axes_path.read_text()); manifest = json.loads(manifest_path.read_text())
    assert axes['native_roles'] == [r[0] for r in ROLES] and axes['replicates_per_mode'] == 200
    keys = axes['original_inputs']; assert len(keys) == len(set(keys)) == 133
    assert [(r['original_input_id'],r['mode'],r['replicate']) for r in manifest] == [
        (key,mode,i) for key in sorted(keys) for mode in sorted(axes['modes']) for i in range(200)]
    assert len({r['case_id'] for r in manifest}) == len(manifest) == 53200
    pins = {}; role_index = {name:i for i,(name,_,_) in enumerate(ROLES)}
    for path in [gate_path,cp,archive_path,plan_path,inputs_plan_path,axes_path,manifest_path,Path(__file__),
                 Path('scripts/matched_predictor_resampling.py'),Path('scripts/readback_matched_predictor_resampling.py')]:
        bind(pins,path)
    def source(path):
        digest = expected_pins.get(str(path),expected_pins.get(str(path.resolve())))
        assert digest is not None, ('Source outside closed archive',str(path))
        assert sha(path) == digest; bind(pins,path,digest)
    for path in [axes_path,manifest_path]: source(path)
    a.output.mkdir(); table = a.output/'paired_branch_distributions.tsv'
    counts=Counter(); marker_ids=set(); taxa=set(); rows=0; arrays=0; branch_slots=0
    with table.open('x') as f:
        writer=None
        for ki,key in enumerate(sorted(keys)):
            config_path=input_root/'inputs'/key/'config.json'; source(config_path)
            config=json.loads(config_path.read_text()); marker_ids.add(config['marker']);taxa.update(config['taxa'])
            splits=axes['splits'][key]; assert len(splits)==len(set(splits))==2*len(config['taxa'])-3
            points={}
            for role,_,_ in ROLES:
                path=point_root/'roles'/(key+'-'+role+'.json');source(path);record=json.loads(path.read_text())
                assert record['status']=='native_point_output_integrity_checked_not_model_qualified'
                points[role]={v['split_mask_hex']:v['length'] for v in record['point_estimate']['branches']}
                assert set(points[role])==set(splits)
            for mode in sorted(axes['modes']):
                values=[]
                for replicate in range(200):
                    path=root/'cases'/key/mode/(f'{replicate:04d}.npz');source(path)
                    draw=independent_draw(len(config['columns']),axes['resampling_groups'][key],mode,replicate)
                    array=read_arrays(path,draw,splits)
                    assert np.all(array['role_status']==0), 'No successful-only summary after unresolved source'
                    values.append(array['branch_lengths']);arrays+=1;branch_slots+=array['branch_lengths'].size
                matrix=np.stack(values)
                assert matrix.shape==(200,7,len(splits))
                for model in MODELS:
                    af,esm='AlphaFold_'+model,'ESMFold_'+model
                    for bi,split in enumerate(splits):
                        result=statistics(matrix[:,0,bi],matrix[:,role_index[af],bi],matrix[:,role_index[esm],bi],
                                          points[af][split],points[esm][split])
                        row=dict(input_id=key,marker=config['marker'],resampling_group=axes['resampling_groups'][key],
                            topology_sha256=config['topology_sha256'],taxa=len(config['taxa']),retained_columns=len(config['columns']),
                            split_mask_hex=split,internal=split in config['internal_splits'],mode=mode,model=model,
                            replicates=200,unit='expected_state_substitutions_per_site',**result,
                            calibrated_effect_interval=False,scientific_eligibility=False)
                        if writer is None:writer=csv.DictWriter(f,fieldnames=list(row),delimiter='\t');writer.writeheader()
                        writer.writerow(row);rows+=1;counts[mode+':'+model+':'+result['empirical_direction']]+=1
            print('paired_predictor_summary_inputs',ki+1,'/133',flush=True)
    assert arrays==53200 and branch_slots==6146000 and len(marker_ids)==71 and len(taxa)==21
    verify(pins)
    receipt=dict(status='complete_full_paired_predictor_empirical_branch_distributions_pending_independent_summary_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(),inputs=133,markers=71,fungal_taxa=21,
        modes=axes['modes'],models=MODELS,paired_draw_arrays=arrays,branch_value_slots=branch_slots,
        table_rows=rows,empirical_direction_counts=dict(counts),table=str(table),table_sha256=sha(table),
        source_hashes=pins,scientific_eligibility=False,calibrated_effect_intervals=False,
        all_eight_aims_incomplete=True,gpu=False,new_native_fits=0,
        scope='Every completed original input/draw/role/branch retained; paired differences are ESMFold minus '
              'AlphaFold under the same state model, topology and sampled columns. Quantiles are descriptive '
              'conditional resampling ranges, not calibrated confidence intervals or significance tests. '
              'Binomial MC bounds quantify draw-count precision for the conditional positive fraction only. '
              'Alternative topologies share resampling groups and branches are dependent; no pooled independent '
              'replicate count or multiple-test discoveries are inferred. Original seven-role draws and '
              'point-fit sources remain immutable. Independent summary replay, coverage and direct-structure '
              'checks plus broad fungal replication remain required.')
    with a.receipt.open('x') as f:json.dump(receipt,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='source_hashes'},indent=2))


if __name__=='__main__': main()
