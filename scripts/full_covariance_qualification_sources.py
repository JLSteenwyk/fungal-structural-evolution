"""Closed complete-grid sources and identities for uniform covariance qualification."""
import gzip
import json
from pathlib import Path
import numpy as np
import pyarrow.parquet as pq
from scipy import sparse
from background_measurement_union_sources import closed_source
from full_entity_operator_sources import KINDS,MODES
from full_expanded_model_design_sources import array_digest,digest,AXES
from full_expanded_model_input_sources import ORDERS,MASKS,NUISANCE,schema
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha

SUMMARY_FIELDS=['logical_cases','model_setting_rows','unique_cohorts','unique_designs',
    'audit_rows','setting_audit_links','audit_status_counts','link_status_counts','trees','loading_modes']
LINK_EXTRA=['loading_mode','tree','audit_id','covariance_disposition','combined_disposition']


def jsonl(path):
    with (gzip.open(path,'rt') if str(path).endswith('.gz') else Path(path).open()) as f:
        for line in f:yield json.loads(line)


def audit_id(contract,design_id,mode,tree):
    return digest(['expanded-uniform-covariance-qualification-v1',contract,design_id,mode,tree])


def load(plan,path):
    bindings=dict(plan['pins']);bind(bindings,path)
    design=closed_source(plan['design_completion'],'complete_verified_full_expanded_model_designs',
        'complete_verified_full_expanded_model_designs_archive',2,bindings)
    bank=closed_source(plan['operator_completion'],'complete_verified_full_entity_operator_bank',
        'complete_verified_full_entity_operator_bank_archive',2,bindings)
    inputs=closed_source(plan['inputs_completion'],'complete_verified_full_expanded_model_inputs',
        'complete_verified_full_expanded_model_inputs_archive',2,bindings)
    configs={}
    for name in ['design','operator','inputs']:
        config_path=plan[name+'_plan'];bind(bindings,config_path)
        configs[name]=json.loads(Path(config_path).read_text())
    root=Path(configs['design']['output']);operator_root=Path(configs['operator']['output']);input_root=Path(configs['inputs']['output'])
    for completed,source_root in [(design,root),(bank,operator_root),(inputs,input_root)]:
        assert completed['producer_receipt']==str(source_root/'receipt.json')
        assert completed['producer_receipt_sha256']==bindings[completed['producer_receipt']]
        assert completed['scientific_eligibility'] is False
    n=design['logical_cases'];assert n==bank['logical_cases']==inputs['logical_cases']==plan['expected']['logical_cases']
    assert design['model_setting_rows']==plan['expected']['model_setting_rows']
    assert design['trees']==bank['trees']==plan['trees']==configs['design']['trees']==configs['operator']['trees']
    ids=json.loads((operator_root/'case_ids.json').read_text());assert len(ids)==len(set(ids))==n
    labels=np.load(operator_root/'block_labels.npy');assert labels.shape==(n,)
    operators={}
    manifest=json.loads((operator_root/'operator_manifest.json').read_text())
    assert len(manifest)==13 and {(r['mode'],r['kind']) for r in manifest}=={(m,k) for m in MODES for k in KINDS}|{('contrast','family_intercept')}
    for r in manifest:
        fp=operator_root/r['path'];assert r['sha256']==bindings[str(fp)]
        z=sparse.load_npz(fp).tocsr();assert list(z.shape)==r['shape'] and z.nnz==r['nnz']
        assert z.shape[0]==n and z.has_canonical_format and np.isfinite(z.data).all()
        operators[r['mode'],r['kind']]=z
    cohorts=json.loads((root/'cohort_manifest.json').read_text())
    assert len(cohorts)==design['unique_cohorts'] and len({r['cohort_id'] for r in cohorts})==len(cohorts)
    parts=json.loads((input_root/'partition_manifest.json').read_text())
    assert [(r['mask'],r['order_contrast']) for r in parts]==[(m,o) for m in MASKS for o in ORDERS]
    arrays={}
    columns=['case_id','input_id',*set(sum(AXES.values(),[])+NUISANCE)]
    for part in parts:
        fp=input_root/part['path'];assert part['sha256']==bindings[str(fp)]
        pf=pq.ParquetFile(fp);assert pf.schema_arrow==schema() and pf.metadata.num_rows==n
        table=pf.read(columns=columns);assert table['case_id'].to_pylist()==ids
        arrays[part['mask'],part['order_contrast']]={name:np.asarray(table[name].to_pylist(),
            dtype='S64' if name in ['case_id','input_id'] else float) for name in columns}
    covariance_plan=configs['operator']['covariance_plan'];assert covariance_plan in bindings
    covariance_root=Path(json.loads(Path(covariance_plan).read_text())['output'])
    from full_entity_operator_sources import table
    case_cov=table(covariance_root/'case_covariance_index.tsv.gz')
    assert [c['case_id'] for c in case_cov]==ids
    pattern_rows=np.asarray([int(c['species_pattern_row']) for c in case_cov]);del case_cov
    factors={}
    for tree in plan['trees']:
        with np.load(covariance_root/(tree+'.npz')) as a:factors[tree]=a['factor'][pattern_rows]
    contract=digest(dict(design_completion=sha(plan['design_completion']),
        operator_completion=sha(plan['operator_completion']),inputs_completion=sha(plan['inputs_completion']),
        schema='expanded-uniform-covariance-qualification-v1',residual='uniform_one',
        loading_modes=MODES,folds=['target_into_uniform_residual','signed_family_zero','unsigned_family_into_contrast_intercept']))
    verify(bindings)
    return dict(root=root,ids=ids,labels=labels,operators=operators,cohorts=cohorts,
        arrays=arrays,factors=factors,design_completion=design,contract=contract),bindings


def cohort_rows(source,record):
    fp=source['root']/record['path'];assert sha(fp)==record['sha256']
    with np.load(fp) as a:assert a.files==['case_rows'];rows=a['case_rows']
    assert rows.dtype==np.dtype('int64') and rows.ndim==1 and len(rows)==record['records']
    assert len(set(rows.tolist()))==len(rows) and np.all((rows>=0)&(rows<len(source['ids'])))
    ids=[source['ids'][i] for i in rows];assert ids==sorted(ids)
    assert array_digest(rows,'<i8')==record['case_rows_sha256']
    assert array_digest(ids,'S64')==record['ordered_case_ids_sha256']
    return rows


def design_matrix(source,cohort,rows,record):
    assert record['cohort_id']==cohort['cohort_id'] and record['mask']==cohort['mask']
    values=source['arrays'][record['mask'],record['order_contrast']]
    columns=['intercept',*AXES[record['sequence_axis']][:record['degree']],*NUISANCE]
    assert record['predictor_columns']==columns and record['records']==len(rows)
    matrix=np.column_stack([np.ones(len(rows)),*[values[k][rows] for k in columns[1:]]])
    assert array_digest(matrix,'<f8')==record['raw_design_sha256']
    assert array_digest(values['input_id'][rows],'S64')==record['ordered_input_ids_sha256']
    return matrix[:,record['active_column_indices']]


def folded_operators(source,rows,mode):
    original={k:source['operators'][mode,k][rows].tocsr() for k in KINDS}
    target=original['target_node'];assert target.nnz==len(rows) and np.all(target.data==1)
    assert np.all(np.diff(target.indptr)==1) and len(np.unique(target.indices))==len(rows)
    family=source['operators']['contrast','family_intercept'][rows].tocsr()
    assert not source['operators']['signed','family'][rows].nnz
    assert not (source['operators']['unsigned','family'][rows]-2*family).nnz
    return {**{k:original[k] for k in ['background_node','model_pair','gene','model']},'family_intercept':family}


def covariance_status(audit):
    statuses=[audit[k]['disposition'] for k in ['raw_diagnostics','reml_diagnostics']]
    return ('qualified_uniform_working_covariance_basis' if all(s=='numerically_independent_covariance_bases' for s in statuses)
        else 'covariance_basis_requires_review')
