"""Exercise real inputs, restart integrity and retained numerical failures."""
import copy
import json
from pathlib import Path
import tempfile
import pandas as pd
from ancestral_chain_attempt import sha,write_json
import run_selected_matched_residuals as runner


def main():
    root=Path('results/model_validation/selected-matched-interval-grid-inputs-20260928-v1')
    prep=json.loads((root/'receipt.json').read_text())
    frame=pd.read_parquet('results/structural_comparisons/refined-working-model-grid-export-20260928-v2/unique_fits.parquet')
    ids=[frame.sort_values('records').iloc[0].fit_input_id,
         frame[frame.selection=='refined_candidate'].sort_values('records').iloc[0].fit_input_id]
    tasks={}
    for line in (root/'tasks.jsonl').open():
        t=json.loads(line)
        if t['fit_input_id'] in ids:tasks[t['fit_input_id']]=t
    runner.initialize(prep['factors'])
    results=[]
    with tempfile.TemporaryDirectory() as temp:
        for identifier in dict.fromkeys(ids):
            rows={r['tree']:r for r in frame[frame.fit_input_id==identifier].to_dict('records')}
            task=tasks[identifier]
            result=runner.execute(task,rows,'results/model_validation/matched-simulation-input-cache-20260928-v1',temp,'fixture-plan')
            path=Path(result['path']); before=path.stat().st_mtime_ns
            repeated=runner.execute(task,rows,'results/model_validation/matched-simulation-input-cache-20260928-v1',temp,'fixture-plan')
            assert result==repeated and path.stat().st_mtime_ns==before
            payload=json.loads(path.read_text())
            assert all(f['status']=='descriptive_marginal_residual_diagnostics' for f in payload['fits'].values())
            results.append(dict(fit_input_id=identifier,records=int(next(iter(rows.values()))['records']),
                selections={tree:r['selection'] for tree,r in rows.items()},result_sha256=result['sha256']))
            bad=copy.deepcopy(payload['binding']);bad['plan_sha256']='changed'
            try:runner.saved_result(path,bad)
            except AssertionError:pass
            else:raise AssertionError('Changed plan accepted')
            original=path.read_text();path.write_text(original+' ')
            try:runner.saved_result(path,payload['binding'])
            except AssertionError:pass
            else:raise AssertionError('Changed payload accepted')
            path.write_text(original)
        # Inject an inconsistent coefficient while retaining source binding: must remain a review row.
        target_tree=next(iter(rows)); bad_rows=copy.deepcopy(rows)
        bad_rows[target_tree]['selected_intercept']+=1
        failure=runner.execute(task,bad_rows,'results/model_validation/matched-simulation-input-cache-20260928-v1',str(Path(temp)/'failure'),'fixture-plan')
        failed=json.loads(Path(failure['path']).read_text())
        assert failed['fits'][target_tree]['status']=='residual_computation_requires_review'
        assert len(failed['fits'])==5
        # A missing digest is never silently reused or overwritten.
        path=Path(failure['path']);path.with_suffix('.sha256').unlink()
        try:runner.saved_result(path,failed['binding'])
        except AssertionError:pass
        else:raise AssertionError('Incomplete checkpoint accepted')
    result=dict(status='passed_real_input_resume_and_failure_accounting_checks',cases=results,
        changed_plan_rejected=True,changed_payload_rejected=True,incomplete_checkpoint_rejected=True,
        numerical_failure_retained=True,pins={str(p):sha(p) for p in [Path(__file__),Path(runner.__file__)]},
        scope='Two real input groups including refined estimates; full residual output audit and simulation-based reference distributions remain required.')
    write_json(Path('metadata/selected_matched_residual_runner_checks_20260928.json'),result)
    print(json.dumps(result))

if __name__=='__main__':main()
