#!/usr/bin/env python3
"""Census every original scalar log without equating finite output with mixing."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha

TRACKED = ['ASRV.Gamma:alpha', 'RS07:rate', 'RS07:meanLength', '|A|', '#indels', '|indels|', 'prior', 'likelihood', 'posterior']
FIELDS = ['chain_id', 'model_input_identity', 'prior', 'family', 'integrity_disposition',
    'variable', 'logged_iterations', 'last_iteration', 'finite_values', 'nan_values',
    'positive_infinity_values', 'negative_infinity_values', 'first_value', 'last_value',
    'finite_minimum', 'finite_maximum', 'first_nonfinite_iteration', 'source_log', 'source_log_sha256']


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--completion',type=Path,required=True);p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--table',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    assert not a.table.exists() and not a.output.exists()
    c=json.loads(a.completion.read_text()); assert c['status']=='complete_verified_baliphy_initial_horizon_accounting'
    archive=Path(c['full_hash_archive']);assert sha(archive)==c['full_hash_archive_sha256']
    proof=json.loads(archive.read_text());assert len(proof['services'])==c['exact_process_journals_checked']==2
    known={str(Path(q).resolve()):d for q,d in proof['source_hashes'].items()}
    bindings={str(a.completion):sha(a.completion),str(archive):sha(archive),str(Path(__file__)):sha(__file__)}
    def check(path,digest=None):
        path=Path(path); expected=known[str(path.resolve())]
        if digest is not None:assert expected==digest
        assert sha(path)==expected,path;bindings[str(path)]=expected
        return path
    plan=json.loads(check(a.plan).read_text());jobs=json.loads(check(plan['jobs']).read_text())
    native=json.loads(check(Path(plan['output'])/'receipt.json').read_text())
    index={r['chain_id']:r for r in native['chains']}; assert len(jobs)==len(index)==1620
    rows=[]; chain_summaries=[]; nonfinite=Counter(); dispositions=Counter(); values_seen=0
    for job in jobs:
        ch=job['chain'];cid=ch['chain_id'];status=index[cid]['status'];dispositions[status]+=1
        rp=check(index[cid]['receipt'],index[cid]['receipt_sha256']);r=json.loads(rp.read_text())
        candidates=[name for name in r['artifacts'] if Path(name).name=='C1.log'];assert len(candidates)==1
        lp=check(rp.parent/candidates[0],r['artifacts'][candidates[0]])
        if status=='all_saved_alignments_and_candidate_nodes_checked':
            audit=json.loads(check(index[cid]['sample_audit'],index[cid]['sample_audit_sha256']).read_text())
            assert Path(audit['scalar_log']).resolve()==lp.resolve() and audit['scalar_log_sha256']==r['artifacts'][candidates[0]]
        else:assert status=='failed' and r['exit_code']!=0
        stats={k:dict(finite=0,nan=0,posinf=0,neginf=0,minimum=None,maximum=None,first=None,last=None,first_nonfinite=None) for k in TRACKED}
        chain_bad=Counter(); iterations=[]; headers=None
        with lp.open() as handle:
            reader=csv.DictReader(handle,delimiter='\t');headers=reader.fieldnames
            assert len(headers)==len(set(headers)) and set(TRACKED+['iter'])<=set(headers)
            for x in reader:
                assert None not in x and None not in x.values()
                iteration=int(x['iter']);iterations.append(iteration)
                for name, raw in x.items():
                    number=float(raw);values_seen+=1
                    if not math.isfinite(number):chain_bad[name]+=1;nonfinite[name]+=1
                    if name not in stats:continue
                    s=stats[name]
                    if s['first'] is None:s['first']=raw
                    s['last']=raw
                    if math.isfinite(number):
                        s['finite']+=1;s['minimum']=number if s['minimum'] is None else min(s['minimum'],number)
                        s['maximum']=number if s['maximum'] is None else max(s['maximum'],number)
                    else:
                        key='nan' if math.isnan(number) else ('posinf' if number>0 else 'neginf');s[key]+=1
                        if s['first_nonfinite'] is None:s['first_nonfinite']=iteration
        assert iterations==list(range(len(iterations)))
        if status=='all_saved_alignments_and_candidate_nodes_checked':assert iterations==list(range(1001))
        assert iterations
        chain_summaries.append(dict(chain_id=cid,status=status,logged_iterations=len(iterations),last_iteration=iterations[-1],nonfinite_by_variable=dict(chain_bad)))
        for name,s in stats.items():
            rows.append(dict(chain_id=cid,model_input_identity=job['config']['model_input_identity'],prior=ch['prior_label'],family=ch['family'],integrity_disposition=status,
                variable=name,logged_iterations=len(iterations),last_iteration=iterations[-1],finite_values=s['finite'],nan_values=s['nan'],positive_infinity_values=s['posinf'],negative_infinity_values=s['neginf'],
                first_value=s['first'],last_value=s['last'],finite_minimum='' if s['minimum'] is None else s['minimum'],finite_maximum='' if s['maximum'] is None else s['maximum'],
                first_nonfinite_iteration='' if s['first_nonfinite'] is None else s['first_nonfinite'],source_log=str(lp),source_log_sha256=r['artifacts'][candidates[0]]))
    assert dispositions=={'all_saved_alignments_and_candidate_nodes_checked':1617,'failed':3} and len(rows)==14580
    # Rehash exactly the immutable logs/receipts/plans used by this scalar census.
    # Saved alignment/property arrays were fully closed by the input accounting;
    # they are not silently claimed to have been re-parsed here.
    for path,digest in bindings.items():assert sha(path)==digest,path
    with a.table.open('x') as handle:
        w=csv.DictWriter(handle,FIELDS,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    with a.table.open() as handle: saved=list(csv.DictReader(handle,delimiter='\t'))
    assert saved==[{key:str(row[key]) for key in FIELDS} for row in rows]
    result=dict(status='complete_full_original_baliphy_native_scalar_range_census',checked_utc=datetime.now(timezone.utc).isoformat(),
        source_completion=str(a.completion),source_completion_sha256=sha(a.completion),original_chains=1620,tracked_variables=TRACKED,range_rows=14580,
        numeric_log_values_checked=values_seen,logged_iterations=sum(x['logged_iterations'] for x in chain_summaries),
        chains_with_nonfinite_scalar_values=sum(bool(x['nonfinite_by_variable']) for x in chain_summaries),
        integrity_checked_chains_with_nonfinite_scalar_values=sum(bool(x['nonfinite_by_variable']) for x in chain_summaries if x['status']=='all_saved_alignments_and_candidate_nodes_checked'),
        nonfinite_by_variable=dict(nonfinite),chain_dispositions=dict(dispositions),chain_summaries=chain_summaries,
        source_hashes=bindings,artifacts={str(a.table):sha(a.table)},scientific_eligibility=False,
        scope='All1620immutable original C1.logs parsed; every43-field numeric value checked for nonfiniteness, nine predeclared variables exported as raw observed ranges and allthreefailed histories retained. Every1001iteration grid checked for integrity-passing chains; failed logs retain their observed partial grid. No threshold-based filtering, parameter repair, inferred memory requirement, causation or posterior acceptance. Extreme ranges flag follow-up work; they are not probability intervals. Savedalignments were source-closed separately, not re-parsed bythis census. New memory-recovery attempts are separate and are not included whilelive.')
    with a.output.open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','chain_summaries','scope']},indent=2),flush=True)


if __name__=='__main__':main()
