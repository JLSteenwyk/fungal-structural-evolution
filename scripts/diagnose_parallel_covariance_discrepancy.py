#!/usr/bin/env python3
"""Read-only source-bound replay of candidate large-cohort numerical failures."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
from scipy import sparse

from ancestral_chain_attempt import sha
from covariance_basis_context import ComponentKernelProducts
from covariance_basis_independent import IndependentKernelProducts
from full_covariance_qualification_sources import jsonl, cohort_rows, design_matrix, folded_operators
from full_expanded_model_design_sources import AXES, digest
from full_expanded_model_input_sources import NUISANCE
from full_entity_operator_sources import table
from readback_full_covariance_qualification import numeric
from reference_measurement_union_sources import bind, verify


def compare(saved, reference, names, n):
    try:
        conservative, error = numeric(saved, reference, names, n)
        result = dict(status='original_frozen_numeric_comparison_passed', conservative=conservative,
                      maximum_projected_absolute_error=error)
    except (AssertionError, ValueError, ArithmeticError) as error:
        result = dict(status='original_frozen_numeric_comparison_failed', error_type=type(error).__name__,
                      error_message=str(error))
    result['gram_mismatches'] = []
    for field, key in [('raw_gram', 'raw'), ('projected_gram', 'projected')]:
        a = np.asarray(saved[field]); b = reference[key]
        for i, j in zip(*np.where(abs(a-b) > 2e-8 + 3e-9*abs(b))):
            result['gram_mismatches'].append(dict(field=field, row=int(i), column=int(j),
                row_kernel=names[i], column_kernel=names[j], saved=float(a[i,j]), independent=float(b[i,j]),
                absolute_difference=float(abs(a[i,j]-b[i,j])),
                original_envelope=float(saved['raw_roundoff_envelope' if field=='raw_gram' else 'projected_roundoff_envelope'][i][j])))
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True); p.add_argument('--receipt', type=Path, required=True)
    args=p.parse_args(); args.output.mkdir(exist_ok=False); assert not args.receipt.exists()
    pp=Path('metadata/parallel_covariance_readback_plan_20261003_v1.json'); plan=json.loads(pp.read_text())
    qpath=Path(plan['source_plan']); q=json.loads(qpath.read_text()); root=Path(q['output'])
    rp=root/'receipt.json'; assert sha(rp)==plan['pins'][str(rp)]
    receipt=json.loads(rp.read_text()); assert receipt['plan_sha256']==sha(qpath)
    assert receipt['scientific_eligibility'] is False
    pins=dict(receipt['source_hashes']); bindings={str(pp):sha(pp),str(qpath):sha(qpath),str(rp):sha(rp)}
    ap=root/'design_covariance_audits.jsonl.gz'; bind(bindings,ap,receipt['artifacts'][ap.name])
    cp=Path('data/software_audits/parallel-covariance-source-discrepancy-20261003-v1/candidate-original-audits.json')
    candidates=json.loads(cp.read_text()); assert len(candidates)==10
    wanted={r['audit_id']:r for r in candidates}; recovered={}
    for record in jsonl(ap):
        if record['audit_id'] in wanted:
            assert record==wanted[record['audit_id']]; recovered[record['audit_id']]=record
            if len(recovered)==len(wanted):break
    assert len(recovered)==len(wanted); bind(bindings,cp)
    configs={}
    for kind in ['design','operator','inputs']:
        fp=Path(q[kind+'_plan']); bind(bindings,fp,pins[str(fp)])
        configs[kind]=json.loads(fp.read_text())
    dr=Path(configs['design']['output']); br=Path(configs['operator']['output']); ir=Path(configs['inputs']['output'])
    def original(path):
        bind(bindings,path,pins[str(path)]); return path
    cohort_manifest=json.loads(original(dr/'cohort_manifest.json').read_text())
    cid=candidates[0]['cohort_id']; cohort=next(c for c in cohort_manifest if c['cohort_id']==cid)
    assert all(r['cohort_id']==cid and r['tree']=='mafft_guide' for r in candidates)
    original(dr/cohort['path']); original(dr/'unique_designs.jsonl')
    recipes={r['design_id']:r for r in jsonl(dr/'unique_designs.jsonl') if r['design_id'] in {r['design_id'] for r in candidates}}
    labels=np.load(original(br/'block_labels.npy'),allow_pickle=False)
    ids=json.loads(original(br/'case_ids.json').read_text())
    operators={}
    for part in json.loads(original(br/'operator_manifest.json').read_text()):
        fp=original(br/part['path']); assert sha(fp)==part['sha256']
        operators[part['mode'],part['kind']]=sparse.load_npz(fp).tocsr()
    manifest=json.loads(original(ir/'partition_manifest.json').read_text()); arrays={}
    pairs={(cohort['mask'],d['order_contrast']) for d in recipes.values()}
    columns=['case_id','input_id',*sorted(set(sum(AXES.values(),[])+NUISANCE))]
    for part in manifest:
        key=(part['mask'],part['order_contrast'])
        if key not in pairs:continue
        fp=original(ir/part['path']); assert sha(fp)==part['sha256']
        values=pq.read_table(fp,columns=columns); assert values['case_id'].to_pylist()==ids
        arrays[key]={k:np.asarray(values[k].to_pylist(),dtype='S64' if k in ['case_id','input_id'] else float) for k in columns}
    vp=original(Path(configs['operator']['covariance_plan'])); vr=Path(json.loads(vp.read_text())['output'])
    indices=table(original(vr/'case_covariance_index.tsv.gz')); assert [r['case_id'] for r in indices]==ids
    patterns=np.asarray([int(r['species_pattern_row']) for r in indices])
    with np.load(original(vr/'mafft_guide.npz'),allow_pickle=False) as a:factor=a['factor'][patterns]
    source=dict(root=dr,ids=ids,labels=labels,operators=operators,arrays=arrays,factors={'mafft_guide':factor})
    rows=cohort_rows(source,cohort); design=design_matrix(source,cohort,rows,recipes[candidates[0]['design_id']])
    assert all(np.array_equal(design,design_matrix(source,cohort,rows,recipes[r['design_id']])) for r in candidates)
    replay=[]
    for mode in ['signed','unsigned']:
        bank=folded_operators(source,rows,mode); names=['residual',*bank,'species']; saved=next(r for r in candidates if r['loading_mode']==mode)
        independent=IndependentKernelProducts(labels[rows],bank).tree(factor[rows])
        reference=independent.project(design); serial=compare(saved['numerical_audit'],reference,names,len(rows))
        fresh=ComponentKernelProducts(labels[rows],bank,np.ones(len(rows))).tree(factor[rows]).audit(design)
        fresh_vs_saved={field:float(np.max(abs(np.asarray(fresh[field])-np.asarray(saved['numerical_audit'][field]))))
                        for field in ['raw_gram','projected_gram','raw_roundoff_envelope','projected_roundoff_envelope']}
        # No mutable producer state is shared; only frozen latent context is read.
        with ThreadPoolExecutor(max_workers=8) as pool:
            threaded=list(pool.map(lambda _:compare(saved['numerical_audit'],independent.project(design),names,len(rows)),range(8)))
        detail=dict(loading_mode=mode,candidate_audit_id=saved['audit_id'],serial=serial,threaded=threaded,
                    fresh_component_vs_original_maximum_differences=fresh_vs_saved,
                    fresh_component_full_audit_equals_saved=fresh==saved['numerical_audit'],
                    serial_matches_each_thread=all(v==serial for v in threaded),source_audit=saved,
                    independent_reference={k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in reference.items()},
                    fresh_component_audit=fresh)
        fp=args.output/(mode+'-replay.json');fp.write_text(json.dumps(detail,indent=2,allow_nan=False)+'\n')
        replay.append({k:v for k,v in detail.items() if k not in ['source_audit','independent_reference','fresh_component_audit','threaded']})
    verify(bindings)
    own=[Path(__file__),Path('scripts/covariance_basis_context.py'),Path('scripts/covariance_basis_independent.py'),
         Path('scripts/readback_full_covariance_qualification.py'),Path('scripts/covariance_basis_audit.py'),
         Path('scripts/full_covariance_qualification_sources.py')]
    bindings.update({str(p):sha(p) for p in own})
    result=dict(status='completed_source_bound_read_only_covariance_discrepancy_diagnostic',checked_utc=datetime.now(timezone.utc).isoformat(),
        cohort_id=cid,records=len(rows),candidate_original_audits=10,unique_fixed_design_matrices=1,
        serial_and_threaded_replays=replay,source_hashes=bindings,artifacts={str(p):sha(p) for p in args.output.iterdir()},
        numeric_tolerances_unchanged=True,original_outputs_unchanged=True,native_inference_runs_started=0,
        scientific_eligibility=False,scope='Diagnostic candidate selection from preserved original traceback, complete '
        'original audit artifact hash and ten exact exported record matches. Selected original plans/cohort/recipes/'
        'CSR/PQ/factor files rehashed against original producer bindings; broader source graph is not freshly replayed. '
        'One distinct fixed design, both modes, one/eight-thread latent calculations and fresh frozen component '
        'calculations. Original numeric assertions remain unchanged. This is not full arithmetic closure, an '
        'accepted numerical repair, production fitting, a biological pilot or a restart of failed native jobs.')
    with args.receipt.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']},indent=2))


if __name__=='__main__':main()
