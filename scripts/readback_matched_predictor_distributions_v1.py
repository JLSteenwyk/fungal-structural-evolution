#!/usr/bin/env python3
"""Independently replay all paired summaries with sorted ranks and scalar arithmetic."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import binom

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import verify


def rank_quantile(vector, p):
    ordered = sorted(map(float, vector)); rank = 199*p
    lo = math.floor(rank); hi = math.ceil(rank)
    return ordered[lo] + (rank-lo)*(ordered[hi]-ordered[lo])


def numeric(row, name, expected, scale=1.0):
    value = float(row[name]); assert math.isfinite(value)
    assert abs(value-expected) <= 8*np.finfo(float).eps*max(1.0,scale,abs(expected)), (name,value,expected)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['source','transport','receipt']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();assert not a.receipt.exists()
    source=json.loads(a.source.read_text());transport=json.loads(a.transport.read_text())
    assert source['status']=='complete_full_paired_predictor_empirical_branch_distributions_pending_independent_summary_readback'
    assert transport['status']=='verified_original_software_wait_and_whole_wrapper_payloads'
    assert source['scientific_eligibility'] is source['calibrated_effect_intervals'] is False
    assert source['table_rows']==13170 and source['paired_draw_arrays']==53200
    assert transport['validation_sha256']==sha(a.source) and transport['original_tool_terminal_exit_code']==0
    assert transport['whole_wrapper_initial_and_terminal_payloads_matched']
    verify(transport['source_hashes']);verify(source['source_hashes'])
    table=Path(source['table']);assert sha(table)==source['table_sha256']
    plan_path=Path('metadata/matched_predictor_resampling_plan_20261003_v1.json')
    plan=json.loads(plan_path.read_text());root=Path(plan['output'])
    axes=json.loads((root/'axes.json').read_text())
    inputs_plan=json.loads(Path(plan['input_plan']).read_text());input_root=Path(inputs_plan['output'])
    role_indices={role:i for i,role in enumerate(axes['native_roles'])}
    counts=Counter();rows=arrays=slots=0
    with table.open() as f:
        reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames
        for ki,key in enumerate(sorted(axes['original_inputs'])):
            config=json.loads((input_root/'inputs'/key/'config.json').read_text());splits=axes['splits'][key]
            points={}
            for role in axes['native_roles']:
                record=json.loads((Path(plan['point_output'])/'roles'/(key+'-'+role+'.json')).read_text())
                points[role]={b['split_mask_hex']:b['length'] for b in record['point_estimate']['branches']}
            for mode in sorted(axes['modes']):
                draw_values=[]
                for replicate in range(200):
                    path=root/'cases'/key/mode/(f'{replicate:04d}.npz')
                    assert sha(path)==source['source_hashes'][str(path)]
                    with np.load(path,allow_pickle=False) as data:
                        assert data['role_status'].tolist()==[0]*7
                        values=data['branch_lengths'].copy()
                        assert values.dtype==np.dtype('<f8') and values.shape==(7,len(splits))
                        assert np.isfinite(values).all() and (values>=0).all()
                    draw_values.append(values);arrays+=1;slots+=values.size
                data=np.stack(draw_values)
                for model in ['af','af_empirical','llm']:
                    afrole,esmrole='AlphaFold_'+model,'ESMFold_'+model
                    for bi,split in enumerate(splits):
                        row=next(reader);rows+=1
                        metadata=dict(input_id=key,marker=config['marker'],resampling_group=axes['resampling_groups'][key],
                            topology_sha256=config['topology_sha256'],taxa=len(config['taxa']),retained_columns=len(config['columns']),
                            split_mask_hex=split,internal=split in config['internal_splits'],mode=mode,model=model,
                            replicates=200,unit='expected_state_substitutions_per_site',
                            calibrated_effect_interval=False,scientific_eligibility=False)
                        assert all(row[k]==str(v) for k,v in metadata.items())
                        aa,af,esm=data[:,0,bi],data[:,role_indices[afrole],bi],data[:,role_indices[esmrole],bi]
                        delta=[float(e)-float(v) for e,v in zip(esm,af)]
                        scale=max(map(abs,delta),default=1.0);mean=math.fsum(delta)/200
                        numeric(row,'delta_mean',mean,scale)
                        numeric(row,'delta_sd',math.sqrt(math.fsum((v-mean)**2 for v in delta)/199),scale)
                        q=[rank_quantile(delta,v) for v in [.025,.5,.975]]
                        for name,value in zip(['delta_q025','delta_median','delta_q975'],q):numeric(row,name,value,scale)
                        pos=sum(v>0 for v in delta);neg=sum(v<0 for v in delta);zero=sum(v==0 for v in delta)
                        assert [int(row[n]) for n in ['positive_draws','negative_draws','zero_draws']]==[pos,neg,zero]
                        numeric(row,'positive_fraction',pos/200)
                        lo,hi=float(row['positive_fraction_mc95_low']),float(row['positive_fraction_mc95_high'])
                        assert 0<=lo<=pos/200<=hi<=1
                        assert lo==0 if pos==0 else abs(binom.sf(pos-1,200,lo)-.025)<1e-10
                        assert hi==1 if pos==200 else abs(binom.cdf(pos,200,hi)-.025)<1e-10
                        direction='positive' if q[0]>0 else 'negative' if q[2]<0 else 'includes_zero'
                        assert row['empirical_direction']==direction
                        assert row['empirical_central_range_excludes_zero']==str(direction!='includes_zero')
                        counts[mode+':'+model+':'+direction]+=1
                        for name,vector in [('aa',aa),('AlphaFold',af),('ESMFold',esm)]:
                            for suffix,prob in [('q025',.025),('median',.5),('q975',.975)]:
                                numeric(row,name+'_'+suffix,rank_quantile(vector,prob),max(vector))
                            assert int(row[name+'_zero_draws'])==sum(v==0 for v in vector)
                        numeric(row,'point_AlphaFold',points[afrole][split])
                        numeric(row,'point_ESMFold',points[esmrole][split])
                        numeric(row,'point_ESMFold_minus_AlphaFold',points[esmrole][split]-points[afrole][split])
                        expected_fields=set(metadata)|{'point_AlphaFold','point_ESMFold','point_ESMFold_minus_AlphaFold',
                            'delta_mean','delta_sd','delta_q025','delta_median','delta_q975','positive_draws','negative_draws',
                            'zero_draws','positive_fraction','positive_fraction_mc95_low','positive_fraction_mc95_high',
                            'empirical_central_range_excludes_zero','empirical_direction'}|{
                                name+'_'+suffix for name in ['aa','AlphaFold','ESMFold'] for suffix in ['q025','median','q975','zero_draws']}
                        assert set(fields)==expected_fields and len(fields)==len(expected_fields)
            print('independent_paired_summary_inputs',ki+1,'/133',flush=True)
        assert next(reader,None) is None
    assert rows==13170 and arrays==53200 and slots==6146000
    assert dict(counts)==source['empirical_direction_counts']
    verify(source['source_hashes']);assert sha(table)==source['table_sha256']
    result=dict(status='passed_full_independent_paired_predictor_distribution_summary_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(),source_receipt_sha256=sha(a.source),
        summary_rows_checked=rows,original_draw_arrays_reread=arrays,branch_value_slots=slots,
        empirical_direction_counts=dict(counts),scientific_eligibility=False,calibrated_effect_intervals=False,
        all_eight_aims_incomplete=True,new_native_fits=0,gpu=False,
        source_hashes={str(q):sha(q) for q in [Path(__file__),a.source,a.transport,table,plan_path]},
        scope='Every summary row independently reconstructed from original paired branch arrays '
              'using sorted linear order statistics and math.fsum/scalar sums. Exact binomial tails '
              'check MC-fraction bounds. All source digests rechecked before/after. Eight machine '
              'epsilons times scale checks arithmetic roundoff only; no original numerical or '
              'biological tolerance is changed. Descriptive conditional ranges remain uncalibrated '
              'for effects, and dependent topology/branch/model summaries are not discoveries.')
    with a.receipt.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
