"""Additional complete synthetic grid with identifiable nonuniform D bases.

Keep the original five-cohort synthetic grid separately. This additional grid
changes synthetic reuse/operator patterns consistently and supplies an uneven
two-component cohort. No real source, case, QC gate or setting is changed.
"""
import csv
import gzip
import json
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
from scipy import sparse

from ancestral_chain_attempt import sha
from check_full_inverse_reuse_weights import closure, write
from covariance_exact_folds_v2 import certificate, variance_map
from full_expanded_model_design_sources import AXES, NUISANCE, OUTCOMES, array_digest, cohort_id, design_id, digest, fit_id
from full_weighted_covariance_sources_v2 import MODES, STATUSES
from inverse_reuse_weight_controls import calculate
from nonuniform_covariance_cone import record as cone_record
from prepare_full_expanded_model_designs import rank_audit
from run_full_inverse_reuse_weights import record as weight_record


def revise_synthetic_source(path, pair_exception=False):
    plan = json.loads(path.read_text())
    roots = {k:Path(json.loads(Path(plan[k + '_plan']).read_text())['output']) for k in STATUSES}
    op = roots['operator']; labels = np.load(op / 'block_labels.npy'); ids = json.loads((op / 'case_ids.json').read_text())
    manifest = json.loads((op / 'operator_manifest.json').read_text())
    operators = {(r['mode'],r['kind']):sparse.load_npz(op / r['path']).tocsr() for r in manifest}
    # Unequal background-node reuse within each component:2,3,7 rows.
    _, component = np.unique(labels, return_inverse=True); position = np.arange(24) % 12
    group = np.where(position < 2, 0, np.where(position < 5, 1, 2)) + 3 * component
    for mode in MODES:
        background = sparse.csr_matrix((np.full(24, -1. if mode == 'signed' else 1.), (np.arange(24),group)), shape=(24,6))
        target = operators[mode,'target_node']
        operators[mode,'background_node'] = background
        pair_target = target
        if pair_exception:
            columns = target.indices.copy(); columns[[1,13]] = columns[[0,12]]
            pair_target = sparse.csr_matrix((np.ones(24),(np.arange(24),columns)),shape=target.shape)
        operators[mode,'model_pair'] = sparse.hstack([pair_target, background], format='csr')
        operators[mode,'gene'] = .5 * sparse.hstack([target,target,background,background], format='csr')
        operators[mode,'model'] = .5 * sparse.hstack([pair_target,pair_target,background,background], format='csr')
    for r in manifest:
        z = operators[r['mode'],r['kind']]; sparse.save_npz(op / r['path'], z)
        r.update(shape=list(z.shape), nnz=z.nnz, sha256=sha(op / r['path']))
    write(op / 'operator_manifest.json', manifest)
    dr = roots['design']; original_cohorts = json.loads((dr / 'cohort_manifest.json').read_text())
    old_designs = [json.loads(l) for l in (dr / 'unique_designs.jsonl').read_text().splitlines()]
    with gzip.open(dr / 'model_settings.tsv.gz','rt') as f:
        reader = csv.DictReader(f,delimiter='\t'); fields = reader.fieldnames; settings = list(reader)
    cohorts = []; old_to_new = {}; selected_by_id = {}
    for c in original_cohorts:
        with np.load(dr / c['path']) as z:selected = z['case_rows']
        old_id = c['cohort_id']; c = dict(c)
        if len(selected) == 24:
            selected = selected[selected != 0]
            ordered = array_digest([ids[i] for i in selected], 'S64')
            cid = cohort_id(c['guide'], c['mask'], ordered); p = dr / 'cohorts' / (cid + '.npz')
            np.savez_compressed(p, case_rows=selected)
            c.update(cohort_id=cid, ordered_case_ids_sha256=ordered, path=str(p.relative_to(dr)),
                records=len(selected), case_rows_sha256=array_digest(selected,'<i8'), sha256=sha(p))
        old_to_new[old_id] = c['cohort_id']; selected_by_id[c['cohort_id']] = selected; cohorts.append(c)
    cohorts.sort(key=lambda c:c['cohort_id']); byc = {c['cohort_id']:c for c in cohorts}
    parts = json.loads((roots['inputs'] / 'partition_manifest.json').read_text())
    arrays = {(r['mask'],r['order_contrast']):pq.read_table(roots['inputs'] / r['path']).to_pydict() for r in parts}
    contract = digest(['synthetic-input-design-contract']); new_designs = []; new_fits = []; did_by_key = {}; fid_by_key = {}
    for c in cohorts:
        selected = selected_by_id[c['cohort_id']]
        for old in old_designs:
            if old_to_new[old['cohort_id']] != c['cohort_id']:continue
            order,axis,degree = old['order_contrast'],old['sequence_axis'],old['degree']
            values = arrays[c['mask'],order]; columns = ['intercept',*AXES[axis][:degree],*NUISANCE]
            x = np.column_stack([np.ones(len(selected)),*[np.asarray(values[k])[selected] for k in columns[1:]]])
            did = design_id(c['cohort_id'],order,axis,degree,contract)
            d = dict(design_id=did,cohort_id=c['cohort_id'],mask=c['mask'],order_contrast=order,sequence_axis=axis,degree=degree,
                predictor_columns=columns,raw_design_sha256=array_digest(x,'<f8'),
                ordered_input_ids_sha256=array_digest(np.asarray(values['input_id'],dtype='S64')[selected],'S64'),**rank_audit(x,degree))
            new_designs.append(d); did_by_key[c['cohort_id'],order,axis,degree] = did
            for outcome in OUTCOMES:
                y = np.asarray(values[outcome])[selected]; yh = array_digest(y,'<f8'); fid = fit_id(did,outcome,yh)
                status = d['disposition'] if d['disposition'] != 'full_rank_design' else 'constant_response_requires_review' if np.ptp(y) == 0 else 'ready_for_working_covariance_fit'
                new_fits.append(dict(fit_input_id=fid,design_id=did,cohort_id=c['cohort_id'],outcome=outcome,
                    disposition=status,records=len(selected),trees=plan['trees'],response_sha256=yh,
                    response_min=float(y.min()),response_max=float(y.max())))
                fid_by_key[did,outcome] = fid,status
    for s in settings:
        cid = old_to_new[s['cohort_id']]; did = did_by_key[cid,s['order_contrast'],s['sequence_axis'],int(s['degree'])]
        fid,status = fid_by_key[did,s['outcome']]
        s.update(cohort_id=cid,design_id=did,fit_input_id=fid,records=str(byc[cid]['records']),disposition=status)
    write(dr / 'cohort_manifest.json', cohorts)
    for name,records in [('unique_designs.jsonl',new_designs),('unique_fit_inputs.jsonl',new_fits)]:
        (dr / name).write_text(''.join(json.dumps(r)+'\n' for r in records))
    with gzip.open(dr / 'model_settings.tsv.gz','wt') as f:
        w = csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(settings)
    wr = roots['weights']; wc = json.loads((wr / 'receipt.json').read_text())['source_contract']
    cc = json.loads((roots['cone'] / 'receipt.json').read_text())['source_contract']
    controls = []; parents = []; cones = []
    for c in cohorts:
        selected = selected_by_id[c['cohort_id']]
        for mode in MODES:
            ops = {k:operators[mode,k][selected].tocsr() for k in ['target_node','background_node','model_pair','gene','model','family']}
            proof = certificate(ops,operators['contrast','family_intercept'][selected].tocsr(),mode,selected)
            parent = dict(cohort_id=c['cohort_id'],cohort_rows_sha256=c['case_rows_sha256'],ordered_case_ids_sha256=c['ordered_case_ids_sha256'],
                certificate=proof,variance_map=variance_map(proof['relations'],mode),raw_reml_basis_qualification_complete=False,scientific_eligibility=False)
            parents.append(parent);cones.append(cone_record(parent,cc))
        counts,weights,diagonal,groups = calculate([group[selected],group[selected],labels[selected]])
        ap = wr / 'arrays' / (c['cohort_id']+'.npz');np.savez_compressed(ap,case_rows=selected,reuse_counts=counts,weights=weights,reciprocal_diagonals=diagonal)
        jp = wr / 'cohorts' / (c['cohort_id']+'.json');write(jp,weight_record(c,wc,counts,weights,diagonal,groups,selected,ap))
        controls.append(dict(cohort_id=c['cohort_id'],records=len(selected),record_path=str(jp.relative_to(wr)),record_sha256=sha(jp),
            array_path=str(ap.relative_to(wr)),array_sha256=sha(ap)))
    write(wr / 'control_manifest.json',controls)
    for kind,name,records in [('exact','cohort_certificates.jsonl',parents),('cone','cohort_cones.jsonl',cones)]:
        (roots[kind] / name).write_text(''.join(json.dumps(r)+'\n' for r in records))
    plan['expected']['case_row_occurrences'] = sum(c['records'] for c in cohorts)
    for kind in STATUSES:
        cp = Path(plan[kind+'_completion']); old = json.loads(cp.read_text()); archive = json.loads(Path(old['full_hash_archive']).read_text())
        summary = archive['summary']
        if kind == 'weights':summary['case_row_occurrences'] = plan['expected']['case_row_occurrences']
        paths = set(map(Path,archive['source_hashes'])) | {p for p in roots[kind].rglob('*') if p.is_file()
            and p != cp and p.name not in ['completion_archive.json','completion.json'] and not p.name.endswith('.archive.json')}
        closure(cp,STATUSES[kind],sorted(paths),summary)
    write(path,plan)
