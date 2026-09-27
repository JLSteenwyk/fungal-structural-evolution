#!/usr/bin/env python3
"""Inventory every audited profile solution as an unconstrained MG94 start."""
import csv,json,hashlib
from pathlib import Path
from audit_local_branch_parameter_profiles import declarations,unchanged_model_text


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    profiles=Path('results/cds/local-mg94-branch-parameter-profiles-20260927-v1')
    fits=Path('results/cds/local-mg94-diagnostics-20260927-v1')
    audit=Path('results/cds/local-mg94-branch-profile-audit-20260927-v1')
    pr=json.loads((profiles/'receipt.json').read_text());ar=json.loads((audit/'receipt.json').read_text())
    assert ar['status']=='passed_full_branch_parameter_profile_artifact_and_numeric_audit'
    assert ar['source_receipt_sha256']==sha(profiles/'receipt.json')
    for name,h in ar['artifacts'].items():assert sha(audit/name)==h
    review=list(csv.DictReader((audit/'case_summary.tsv').open(),delimiter='\t'))
    assert len(review)==ar['cases']==pr['cases']==1632
    source_cases={x['case_id']:x['receipt_sha256'] for x in json.loads((fits/'receipt.json').read_text())['case_receipts']}
    assert set(source_cases)=={r['case_id'] for r in review}
    out=Path('results/cds/local-mg94-unconstrained-starts-20260927-v1');out.mkdir(exist_ok=False)
    manifest=[];seed_count=0;identical=0;case_summaries=[]
    with (out/'parameter_seeds.jsonl').open('w') as seedfile:
        for item in pr['case_receipts']:
            case=item['case_id'];folder=profiles/case
            assert sha(folder/'receipt.json')==item['receipt_sha256']
            cr=json.loads((folder/'receipt.json').read_text())
            assert sha(fits/case/'receipt.json')==source_cases[case]==cr['source_fit_receipt_sha256']
            baseline=json.loads((fits/case/'receipt.json').read_text())
            assert sha(fits/case/'fit.bf')==baseline['artifacts']['fit.bf']
            original=(fits/case/'fit.bf').read_text();original_values,fixed=declarations(original)
            assert not fixed
            bylabel={('unconstrained' if r['fit_kind']=='unconstrained_reoptimization' else 'factor_'+str(r['branch_parameter_multiplier']).replace('.','p')):r for r in cr['rows']}
            assert len(bylabel)==len(cr['proofs'])==8
            distinct=set()
            for proof in cr['proofs']:
                label=proof['fit_kind'];path=folder/label/'profile_fit.bf'
                assert sha(path)==proof['artifacts']['profile_fit.bf']
                text=path.read_text();values,fixed=declarations(text)
                assert set(values)==set(original_values)
                assert unchanged_model_text(text,values)==unchanged_model_text(original,original_values)
                row=bylabel[label]
                target=[n for n in values if n.endswith('.tree_0.'+row['target_node']+'.t')];assert len(target)==1
                assert fixed==(set() if label=='unconstrained' else {target[0]})
                serialized=json.dumps(values,sort_keys=True,separators=(',',':'))
                seed_hash=hashlib.sha256(serialized.encode()).hexdigest();distinct.add(seed_hash)
                seedfile.write(json.dumps(dict(case_id=case,start_label=label,parameter_sha256=seed_hash,parameters=values),sort_keys=True)+'\n')
                manifest.append(dict(case_id=case,start_label=label,original_free_model=str(fits/case/'fit.bf'),original_free_model_sha256=sha(fits/case/'fit.bf'),
                                     profile_seed_model=str(path),profile_seed_model_sha256=sha(path),parameter_sha256=seed_hash,parameters=len(values),
                                     starting_log_likelihood=row['log_likelihood'],target_parameter_name=target[0],starting_target_value=values[target[0]],
                                     constraints_to_release=len(fixed)))
                seed_count+=len(values)
            identical+=8-len(distinct)
            case_summaries.append(dict(case_id=case,starts=8,distinct_parameter_seeds=len(distinct)))
    assert len(manifest)==13056 and len(case_summaries)==1632
    for name,rows in [('start_manifest.tsv',manifest),('case_start_counts.tsv',case_summaries)]:
        with (out/name).open('w') as handle:
            writer=csv.DictWriter(handle,list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
    flagged=[r for r in review if float(r['best_grid_exceeds_unconstrained_by'])>1e-5]
    with (out/'grid_better_than_unconstrained.tsv').open('w') as handle:
        writer=csv.DictWriter(handle,list(review[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(flagged)
    assert len(flagged)==ar['cases_grid_beats_unconstrained_by_gt_1e_5']==10
    # Full serialized-seed readback against independently retained source values and keys.
    expected={(r['case_id'],r['start_label']):r for r in manifest};seen=set()
    for line in (out/'parameter_seeds.jsonl').open():
        row=json.loads(line);key=row['case_id'],row['start_label'];assert key not in seen;seen.add(key)
        entry=expected[key];values,_=declarations(Path(entry['profile_seed_model']).read_text())
        assert row['parameters']==values
        assert hashlib.sha256(json.dumps(row['parameters'],sort_keys=True,separators=(',',':')).encode()).hexdigest()==entry['parameter_sha256']
    assert seen==set(expected)
    result=dict(status='complete_full_local_mg94_unconstrained_start_inventory',cases=1632,starts=13056,
                saved_parameter_values=seed_count,repeated_seed_rows_retained=identical,
                cases_grid_better_than_unconstrained=len(flagged),maximum_grid_likelihood_advantage=max(float(r['best_grid_exceeds_unconstrained_by']) for r in review),
                source_profile_receipt_sha256=sha(profiles/'receipt.json'),source_audit_sha256=sha(audit/'receipt.json'),
                source_fit_receipt_sha256=sha(fits/'receipt.json'),script_sha256=sha(__file__),
                artifacts={p.name:sha(p) for p in out.iterdir()},
                scope='All eight audited seeds for every case; all parameter values and unchanged model content checked. Future fits must load original unconstrained model then assign these values, verify starting likelihood, and optimize every free parameter. No optimization performed here; ten concerns are not the only cases scheduled. No intervals, saturation certificate or selection conclusion.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
