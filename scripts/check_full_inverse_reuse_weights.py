#!/usr/bin/env python3
"""Complete declared export/SQL readback and rehashed malformed-control checks."""
import argparse
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path

import numpy as np

from ancestral_chain_attempt import sha
from full_expanded_model_design_sources import array_digest, cohort_id, digest
from full_inverse_reuse_weight_sources import SOURCES, SHARED, load
from independent_inverse_reuse_weights import database, reconstruct, check_policy_records
from inverse_reuse_weight_controls import ARRAYS, POLICIES, calculate, policy_records
from reference_measurement_union_sources import verify
from run_full_inverse_reuse_weights import run


def write(p, value):
    p.parent.mkdir(parents=True, exist_ok=True); p.write_text(json.dumps(value, sort_keys=True) + '\n')


def rejected(action):
    try: action()
    except (AssertionError, ValueError, KeyError, FileExistsError): return
    raise AssertionError('Malformed inverse-reuse controls accepted')


def closure(cp, status, paths, summary):
    ap = cp.with_suffix('.archive.json')
    write(ap, dict(status=status + '_archive', services=[{'synthetic': True}] * 2,
                   source_hashes={str(p): sha(p) for p in paths}, summary=summary))
    write(cp, dict(status=status, **summary, full_hash_archive=str(ap),
                   full_hash_archive_sha256=sha(ap), bound_source_hashes=len(paths),
                   exact_process_journals_checked=2, scientific_eligibility=False))


def setup(root):
    pp = {k: root / (k + '.plan.json') for k in SOURCES}
    cp = {k: root / (k + '.completed.json') for k in SOURCES}
    roots = {k: root / k for k in SOURCES}
    for p in roots.values(): p.mkdir(parents=True)
    n = 24; ids = [digest(['synthetic-logical-case', i]) for i in range(n)]
    assert ids != sorted(ids)
    guide = ['mafft'] * 12 + ['profile'] * 12
    keys = {'background_node': [digest(['node', g, (i % 12) // 3]) for i, g in enumerate(guide)],
            'background_pair': [digest(['versioned-pair', g, (i % 12) // 4]) for i, g in enumerate(guide)],
            'family_component': [digest(['connected-component', g, (i % 12) // 5]) for i, g in enumerate(guide)]}
    common = [dict(case_id=ids[i], physical_case_id=digest(['physical', i]),
        target_id=digest(['target', i]), background_id=keys['background_node'][i], guide=guide[i],
        target_family='OGsynthetic', background_family='OGsynthetic', selection_records=str(i + 1),
        target_same_model='0', background_same_model='0') for i in range(n)]
    def table(p, extra, values):
        fields = SHARED + [extra]
        with gzip.open(p, 'wt') as f:
            f.write('\t'.join(fields) + '\n')
            for i, r in enumerate(common): f.write('\t'.join([r[k] for k in SHARED] + [values[i]]) + '\n')
    case = roots['cases'] / 'case_index.tsv.gz'; table(case, 'background_pair_key', keys['background_pair'])
    cov = roots['covariance'] / 'case_covariance_index.tsv.gz'; table(cov, 'family_component', keys['family_component'])
    ip = roots['operator'] / 'case_ids.json'; write(ip, ids)
    inputp = root / 'inputs.plan.json'
    write(inputp, dict(cases_plan=str(pp['cases']), covariance_plan=str(pp['covariance'])))
    configs = {k: dict(output=str(roots[k])) for k in SOURCES}
    configs['operator'].update(covariance_plan=str(pp['covariance']), covariance_completion=str(cp['covariance']))
    configs['covariance'].update(case_plan=str(pp['cases']), case_completion=str(cp['cases']))
    configs['design']['inputs_plan'] = str(inputp)
    for k, v in configs.items(): write(pp[k], v)
    cohorts = []
    for g, mask, selected in [('mafft','full',list(range(12))), ('mafft','plddt70',list(range(6))),
        ('mafft','full',list(range(4))), ('profile','full',list(range(12,24))),
        ('profile','plddt70',[13,15,16,19,20,22])]:
        selected = np.asarray(sorted(selected, key=lambda i: ids[i]), dtype=np.int64)
        ih = array_digest([ids[i] for i in selected], 'S64'); cid = cohort_id(g, mask, ih)
        relative = 'cohorts/' + cid + '.npz'; p = roots['design'] / relative; p.parent.mkdir(exist_ok=True)
        np.savez_compressed(p, case_rows=selected)
        cohorts.append(dict(cohort_id=cid, guide=g, mask=mask, ordered_case_ids_sha256=ih,
            membership_occurrences=4, path=relative, sha256=sha(p), records=len(selected),
            case_rows_sha256=array_digest(selected, '<i8')))
    cohorts.sort(key=lambda c: c['cohort_id']); mp = roots['design'] / 'cohort_manifest.json'; write(mp, cohorts)
    cones = roots['cone'] / 'cohort_cones.jsonl'
    with cones.open('w') as f:
        for c in cohorts:
            for mode in ['signed', 'unsigned']:
                f.write(json.dumps(dict(cohort_id=c['cohort_id'], loading_mode=mode, records=c['records'],
                    cohort_rows_sha256=c['case_rows_sha256'], ordered_case_ids_sha256=c['ordered_case_ids_sha256'],
                    scientific_eligibility=False, variance_map=dict(residual_diagonal_prepared=False,
                    raw_reml_basis_qualification_complete=False))) + '\n')
    occurrences = sum(c['records'] for c in cohorts)
    source_paths = {'cases':[case], 'covariance':[cov], 'operator':[ip],
        'design':[mp, inputp, *[roots['design']/c['path'] for c in cohorts]], 'cone':[cones]}
    summaries = {k: dict(logical_cases=n) for k in SOURCES}
    summaries['design'].update(unique_cohorts=5, cohort_member_occurrences=occurrences)
    summaries['cone'].update(cohorts=5, certificates=10, case_row_occurrences=2 * occurrences)
    for k in SOURCES: closure(cp[k], SOURCES[k], [pp[k], *source_paths[k]], summaries[k])
    modules = ['inverse_reuse_weight_controls', 'independent_inverse_reuse_weights',
               'full_inverse_reuse_weight_sources', 'run_full_inverse_reuse_weights',
               'check_full_inverse_reuse_weights', 'full_exact_covariance_sources',
               'full_expanded_model_design_sources', 'reference_measurement_union_sources']
    plan = {k + '_plan':str(p) for k,p in pp.items()}
    plan.update({k + '_completion':str(p) for k,p in cp.items()})
    plan.update(expected=dict(logical_cases=n, cohorts=5, case_row_occurrences=occurrences),
        output=str(root/'export'), resources=dict(minimum_free_disk_gib=0),
        pins={'scripts/'+m+'.py':sha('scripts/'+m+'.py') for m in modules},
        scope='Complete declared synthetic 24-case/5-cohort source/export/readback. Parent closures and journals are synthetic. No real weighted qualification or biological acceptance.')
    planp = root / 'plan.json'; write(planp, plan)
    return planp


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True); parser.add_argument('--receipt', type=Path, required=True)
    a = parser.parse_args(); root = a.output.resolve(); root.mkdir(exist_ok=False); assert not a.receipt.exists()
    # Independent random/unequal/equal group counts across all controls.
    rng = np.random.default_rng(20261003); mathematical = 0
    for n in [1, 2, 3, 7, 12, 31, 100, 257]:
        for kind in range(5):
            vectors = [np.asarray([str(v) for v in (np.zeros(n, dtype=int) if kind == 0 else
                np.arange(n) if kind == 1 else np.arange(n) // 3 if kind == 2 else
                rng.integers(0, max(1, n // 4), size=n))], dtype='S64') for _ in range(3)]
            src = dict(ids=list(range(n)), keys=dict(zip(POLICIES[1:], vectors)))
            db = database(src); independent = reconstruct(db, np.arange(n)); db.close()
            direct = calculate(vectors)
            for x,y in zip(direct[:3], independent[:3]): assert np.array_equal(x,y)
            assert direct[3] == independent[3]
            records = policy_records(n, *direct)
            check_policy_records(records,n,*independent); mathematical += 1
    rejected(lambda: calculate([[],[],[]])); rejected(lambda: calculate([[1],[1,2],[1]]))
    rejected(lambda: calculate([[[1]],[[1]],[[1]]]))
    pp = setup(root / 'fixture'); plan = json.loads(pp.read_text()); out = Path(plan['output'])
    produced = run(pp); checked = run(pp, True)
    assert (produced['logical_cases'], produced['cohorts'], produced['case_row_occurrences']) == (24,5,40)
    for p in POLICIES[1:]:
        assert produced['policy_census'][p]['exact_uniform_one_cohorts'] > 0
        assert produced['policy_census'][p]['nonuniform_cohorts'] > 0
    rejected(lambda: run(pp)); rejected(lambda: run(pp, True))
    manifestp = out / 'control_manifest.json'; receiptp = out / 'receipt.json'; readp = out / 'readback.json'
    manifest = json.loads(manifestp.read_text()); jp = out / manifest[0]['record_path']; ap = out / manifest[0]['array_path']
    originals = {p:p.read_bytes() for p in [manifestp, receiptp, readp, jp, ap]}; readp.unlink()
    output_cases = ['omit_cohort','duplicate_cohort','case_order','reuse_count','weight','diagonal',
        'normalization','uniform_class','basis_disposition','selection_multiplicity','confidence_precision',
        'claim_ess','promote_science','wrong_membership','extra_array','array_dtype','foreign_contract']
    for name in output_cases:
        m = json.loads(originals[manifestp]); rec = json.loads(originals[jp])
        with np.load(ap,allow_pickle=False) as z: arrays={k:z[k].copy() for k in z.files}
        if name=='omit_cohort':m.pop()
        elif name=='duplicate_cohort':m[-1]=m[0]
        elif name=='case_order':arrays['case_rows']=arrays['case_rows'][::-1]
        elif name=='reuse_count':arrays['reuse_counts'][0,0]+=1
        elif name=='weight':arrays['weights'][1,0]*=1.01
        elif name=='diagonal':arrays['reciprocal_diagonals'][2,0]*=1.01
        elif name=='normalization':rec['policy_records'][1]['exact_weight_sum']+=1
        elif name=='uniform_class':rec['policy_records'][2]['diagonal_is_exact_uniform_one']=not rec['policy_records'][2]['diagonal_is_exact_uniform_one']
        elif name=='basis_disposition':rec['policy_records'][1]['future_basis_disposition']='numerically_accepted'
        elif name=='selection_multiplicity':rec['recipe']['selection_record_multiplicity_used']=True
        elif name=='confidence_precision':rec['calibrated_measurement_precision']=True
        elif name=='claim_ess':rec['effective_sample_size_estimated']=True
        elif name=='promote_science':rec['scientific_eligibility']=True
        elif name=='wrong_membership':rec['cohort_rows_sha256']='foreign'
        elif name=='extra_array':arrays['unexpected']=np.ones(1)
        elif name=='array_dtype':arrays['reuse_counts']=arrays['reuse_counts'].astype(float)
        else:rec['source_contract']='foreign'
        np.savez_compressed(ap,**arrays)
        rec['array_file_sha256']=sha(ap)
        rec['array_sha256']={k:array_digest(v,'<i8' if k in ARRAYS[:2] else '<f8') for k,v in arrays.items()}
        write(jp,rec)
        m[0]['record_sha256']=sha(jp);m[0]['array_sha256']=sha(ap);write(manifestp,m)
        r=json.loads(originals[receiptp])
        for p in [manifestp,jp,ap]:r['artifacts'][str(p.relative_to(out))]=sha(p)
        write(receiptp,r);rejected(lambda:run(pp,True))
        for p in [manifestp,receiptp,jp,ap]:p.write_bytes(originals[p])
    # Deliberately rehash mutated parent artifacts and archives: semantic guards must reject them.
    casescp=Path(plan['cases_completion']); closurevalue=json.loads(casescp.read_text())
    archive=Path(closurevalue['full_hash_archive']);casefile=Path(json.loads(Path(plan['cases_plan']).read_text())['output'])/'case_index.tsv.gz'
    opcp=Path(plan['operator_completion']);opvalue=json.loads(opcp.read_text());opa=Path(opvalue['full_hash_archive'])
    idp=Path(json.loads(Path(plan['operator_plan']).read_text())['output'])/'case_ids.json'
    original_sources={p:p.read_bytes() for p in [casescp,archive,casefile,opcp,opa,idp]}
    source_cases=['foreign_completion','promote_source','missing_case','foreign_pair','broken_join','reordered_ids']
    for name in source_cases:
        c=json.loads(original_sources[casescp]);arch=json.loads(original_sources[archive])
        if name=='foreign_completion':c['status']='foreign'
        elif name=='promote_source':c['scientific_eligibility']=True
        elif name=='reordered_ids':
            write(idp,list(reversed(json.loads(original_sources[idp]))))
            oa=json.loads(original_sources[opa]);oa['source_hashes'][str(idp)]=sha(idp);write(opa,oa)
            oc=json.loads(original_sources[opcp]);oc['full_hash_archive_sha256']=sha(opa);write(opcp,oc)
        else:
            lines=gzip.decompress(original_sources[casefile]).decode().splitlines();fields=lines[0].split('\t')
            if name=='missing_case':lines.pop()
            else:
                cells=lines[1].split('\t'); cells[fields.index('background_pair_key' if name=='foreign_pair' else 'target_id')]='foreign';lines[1]='\t'.join(cells)
            casefile.write_bytes(gzip.compress(('\n'.join(lines)+'\n').encode(),mtime=0))
            arch['source_hashes'][str(casefile)]=sha(casefile);write(archive,arch);c['full_hash_archive_sha256']=sha(archive)
        write(casescp,c);rejected(lambda:load(plan,pp))
        for p,b in original_sources.items():p.write_bytes(b)
    replayed=run(pp,True);assert readp.read_bytes()==originals[readp]
    for p,b in original_sources.items():assert p.read_bytes()==b
    verify(produced['source_hashes']);verify(replayed['source_hashes'])
    hashes=dict(replayed['source_hashes']);hashes.update({str(p):sha(p) for p in [readp,pp,Path(__file__)]})
    result=dict(status='passed_complete_declared_inverse_reuse_controls_sql_fraction_contracts_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(),mathematical_grouping_cases=mathematical,
        invalid_primitive_inputs_rejected=3,synthetic_logical_cases=24,synthetic_cohorts=5,
        synthetic_case_row_occurrences=40,synthetic_case_control_occurrences=160,
        rehashed_output_cases_rejected=output_cases,rehashed_source_cases_rejected=source_cases,
        completed_producer_reader_restarts_rejected=True,byte_exact_positive_restoration=True,
        source_and_journal_fixtures_synthetic=True,source_hashes=hashes,scientific_eligibility=False,
        scope=plan['scope'])
    with a.receipt.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'}),flush=True)


if __name__=='__main__':main()
