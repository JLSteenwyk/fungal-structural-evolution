#!/usr/bin/env python3
"""Inventory both nonlinear grids; deduplicate exact expanded ordered inputs."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from screen_duplication_domain_alignment_coverage import sha


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args(); plan=json.loads(args.plan.read_text()); ph=sha(args.plan)
    def verify():
        assert sha(args.plan)==ph
        for path,digest in plan['pins'].items():assert sha(path)==digest,path
    verify()
    nonlinear=Path(plan['nonlinear'])
    nr=json.loads((nonlinear/'receipt.json').read_text())
    assert nr['status']=='complete_polynomial_identity_contrasts_full_dual_implementation_check'
    for name,digest in nr['artifacts'].items():assert sha(nonlinear/name)==digest
    source=Path(plan['summaries'])
    partitions=json.loads((source/'partition_manifest.json').read_text())
    nodes=pd.read_csv(plan['nodes'],sep='\t')
    targets=nodes[nodes.role.eq('target')][['node_id','guide','family_component']].rename(columns={'node_id':'target_id'})
    pair_patterns=pd.read_csv(plan['pairs'],sep='\t',usecols=['target_id','background_id','species_pattern_id'])
    pair_patterns=pair_patterns.merge(targets,on='target_id',validate='many_to_one')
    assert len(pair_patterns)==52675
    # Keys retain exact target/control identities, including potential future
    # dependence associated with these nodes; do not merge just equal values.
    pair_patterns['row_identity']=[hashlib.sha256(json.dumps(list(row),separators=(',',':')).encode()).hexdigest()
                                  for row in pair_patterns[['target_id','background_id','family_component','species_pattern_id']].itertuples(index=False,name=None)]
    selected=pd.read_csv(plan['selections'],sep='\t',usecols=['target_id','background_id','policy','scenario_id','domain_config_id'])
    selected=selected.merge(pair_patterns,on=['target_id','background_id'],validate='many_to_one')
    assert len(selected)==2786912
    assert not selected.duplicated(['target_id','policy','scenario_id']).any()
    selected=selected.sort_values('target_id',kind='stable')
    summary=pd.read_csv(source/'record_summary.tsv',sep='\t',dtype={'target_order':str,'background_order':str})
    keys=['guide','policy','scenario_id']; setting=['boundary','mask','cohort','screen','target_order','background_order']
    all_numeric=['rmsd_difference','identity_difference','original_coverage_difference','log_aligned_length_ratio','confidence_fraction_difference','identity_power_2_difference','identity_power_3_difference']
    out=Path(plan['output']);out.mkdir(exist_ok=False)
    manifest=[]; recipes={}; signatures={}
    for part in partitions:
        path=source/part['path'];assert sha(path)==part['sha256']
        values=pd.read_parquet(path,columns=['domain_config_id']+all_numeric[:5])
        extra=pd.read_parquet(nonlinear/f"{part['index']:03d}.parquet")
        assert set(values.domain_config_id)==set(extra.domain_config_id)
        values=values.merge(extra[['domain_config_id']+all_numeric[5:]],on='domain_config_id',validate='one_to_one')
        records=selected.merge(values,on='domain_config_id',validate='many_to_one',sort=False)
        expected=summary
        for key in setting:expected=expected[expected[key].eq(str(part[key]))]
        counts=expected.set_index(keys).matched_records.to_dict()
        assert len(counts)==432
        for group,frame in records.groupby(keys,sort=True):
            assert len(frame)==counts[group] and frame.target_id.is_monotonic_increasing
            for degree in [2,3]:
                numeric=all_numeric[:degree+4]
                matrix=np.ascontiguousarray(frame[numeric].to_numpy(),dtype='<f8')
                assert np.isfinite(matrix).all()
                # Canonicalize signed zeros only; no decimal rounding or tolerance-based merging.
                matrix[matrix==0]=0.
                identities=np.asarray(frame.row_identity,dtype='S64')
                value_hash=hashlib.sha256(matrix.tobytes()).hexdigest()
                identity_hash=hashlib.sha256(identities.tobytes()).hexdigest()
                specification=dict(polynomial_degree=degree,records=len(frame),numeric_columns=numeric,values_sha256=value_hash,
                                   ordered_identity_sha256=identity_hash,model_specification=plan['model_specification'])
                serialized=json.dumps(specification,sort_keys=True,separators=(',',':'))
                fit_id=hashlib.sha256(serialized.encode()).hexdigest()
                if fit_id in signatures:assert signatures[fit_id]==serialized
                else:
                    signatures[fit_id]=serialized
                    recipes[fit_id]=dict(fit_input_id=fit_id,partition_index=part['index'],
                                        **dict(zip(keys,group)),**specification)
                manifest.append(dict(polynomial_degree=degree,fit_input_id=fit_id,partition_index=part['index'],
                                     **dict(zip(keys,group)),**{key:part[key] for key in setting},records=len(frame)))
        print('Inventoried full fit setting',part['index']+1,'/ 192',flush=True)
    assert len(manifest)==165888
    table=pd.DataFrame(manifest)
    table.to_csv(out/'full_setting_fit_map.tsv',sep='\t',index=False)
    with (out/'unique_fit_recipes.jsonl').open('w') as handle:
        for key in sorted(recipes):handle.write(json.dumps(recipes[key],sort_keys=True,separators=(',',':'))+'\n')
    verify()
    receipt=dict(status='complete_full_nonlinear_fit_inventory_pending_independent_readback',plan_sha256=ph,
                 full_settings=len(table),unique_record_inputs=len(recipes),tree_alternatives=5,
                 full_tree_fits=len(table)*5,unique_tree_fits=len(recipes)*5,
                 record_occurrences=int(table.records.sum()),
                 unique_input_record_occurrences=sum(r['records'] for r in recipes.values()),
                 artifacts={p.name:sha(p) for p in out.iterdir()},
                 scope='All settings retained; exact numeric and ordered node/family/species identities define reusable likelihood inputs. Only signed zero is canonicalized. Quadratic and cubic specifications and five trees remain separate. Tree counts are candidate workloads, not scheduled fits; optimizer and likelihood-comparison protocols remain to be specified. No approximate matching, biological subsampling, effect selection, fitted model or inference. Complete source checks and recipe readback are required before reuse.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2),flush=True)


if __name__=='__main__':main()
