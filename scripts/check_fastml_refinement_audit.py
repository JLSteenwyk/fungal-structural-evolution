#!/usr/bin/env python3
"""Exercise full-readback checks on frozen completed records and altered controls."""
import copy
from collections import defaultdict
import json
from pathlib import Path
from ancestral_chain_attempt import sha,write_json
from audit_full_fastml_refinement import check_record,summarize_group


def main():
    plan_path=Path('metadata/fastml_optimizer_refinement_plan_20260928_v3.json');plan=json.loads(plan_path.read_text())
    jobs={j['id']:j for j in json.loads((Path(plan['preparation'])/'jobs.json').read_text())}
    checked={};groups=defaultdict(list);example=None
    for path in sorted(Path(plan['output']).glob('*/readback.json')):
        r=json.loads(path.read_text());check_record(jobs[r['id']],r);checked[str(path)]=sha(path);groups[r['job']['job_id']].append(r)
        if r['status']=='refined_native_settings_and_replay_checked_requires_optimization_review':example=r
    assert example is not None
    for kind in ['changed_parameter','changed_likelihood','wrong_probability_count']:
        bad=copy.deepcopy(example)
        if kind=='changed_parameter':bad['fitted_parameters']['alpha']+=1
        elif kind=='changed_likelihood':bad['native_rate_likelihood_difference']+=1
        else:bad['probability_rows']+=1
        try:check_record(jobs[bad['id']],bad)
        except AssertionError:pass
        else:raise AssertionError('Accepted '+kind)
    complete=0
    for key,rows in groups.items():
        if len(rows)!=5:continue
        baseline=[]
        for variant in ['precision-only','precision-cache-refresh']:
            path=Path('results/ancestral/full-fastml-native-rate-replay-20260928-v1')/(key+'__'+variant+'.json')
            if path.exists():baseline.append(json.loads(path.read_text())['native_rate_replay_log_likelihood'])
        summarize_group(key,rows,baseline);complete+=1
        try:summarize_group(key,rows[:-1],baseline)
        except AssertionError:pass
        else:raise AssertionError('Accepted missing start')
    out=Path('metadata/fastml_full_refinement_audit_checks_20260928.json');assert not out.exists()
    write_json(out,dict(status='completed_record_checks_and_tamper_rejections_passed',records=len(checked),complete_five_start_groups=complete,
        rejected_mutations=['changed_parameter','changed_likelihood','wrong_probability_count','missing_start'],
        pins={str(p):sha(p) for p in [plan_path,Path(__file__),Path('scripts/audit_full_fastml_refinement.py')]},source_output_hashes=checked,
        scope='Frozen completed subset, not full production qualification. Full780 audit waits for producer terminal success.'))
    print('Checked',len(checked),'records and',complete,'complete groups; altered fields and missing starts rejected.')


if __name__=='__main__':main()
