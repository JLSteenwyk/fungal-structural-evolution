#!/usr/bin/env python3
"""Check every merged FCS frame field against its selected source row."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()


def rows(path):
    with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    plan=json.loads(a.plan.read_text()); paths={k:Path(v) for k,v in plan['paths'].items()}
    pins={str(a.plan):sha(a.plan),**plan['pins']}
    result_path=paths['output']/'receipt.json';r=json.loads(result_path.read_text())
    pins[str(result_path)]=sha(result_path)
    pins.update({str(paths['output']/n):d for n,d in r['artifacts'].items()})
    def verify():
        for name,d in pins.items():
            if sha(name)!=d:raise ValueError('Changed source: '+name)
    verify()
    disposition=rows(paths['inputs']/'marker_disposition.tsv')
    ready={x['marker']:x for x in disposition if x['baseline_status']=='ready_for_inference'}
    changed={m for m,x in ready.items() if x['sensitivity_status']=='ready_for_inference'}
    unchanged={m for m,x in ready.items() if x['sensitivity_status']=='unchanged_reuse_baseline_fit'}
    if changed&unchanged or changed|unchanged!=set(ready):raise ValueError('Unexpected disposition')
    baseline=rows(paths['baseline']/'site_rate_exposure.tsv')
    update=rows(paths['changed']/'site_rate_exposure.tsv')
    key=lambda x:(x['marker'],x['paired_column_1based'])
    by_key={key(x):x for x in update}
    if len(by_key)!=len(update) or set(by_key)!={key(x) for x in baseline if x['marker'] in changed}:
        raise ValueError('Changed marker coordinate grid differs')
    expected=[by_key[key(x)] if x['marker'] in changed else x for x in baseline]
    actual=rows(paths['output']/'site_rate_exposure.tsv')
    if expected!=actual:raise ValueError('Merged site fields differ')
    expected_diagnostics=[x for x in rows(paths['baseline']/'fit_diagnostics.tsv') if x['marker'] in unchanged]
    expected_diagnostics += [x for x in rows(paths['changed']/'fit_diagnostics.tsv') if x['marker'] in changed]
    if expected_diagnostics!=rows(paths['output']/'fit_diagnostics.tsv'):
        raise ValueError('Merged fit diagnostics differ')
    provenance=[dict(marker=m,frame_source='changed' if m in changed else 'baseline',
                     source_receipt_sha256=sha(paths['changed' if m in changed else 'baseline']/'receipt.json')) for m in sorted(ready)]
    if provenance!=rows(paths['output']/'marker_provenance.tsv'):raise ValueError('Marker provenance differs')
    counts=Counter();flags=Counter()
    for x in actual:
        counts[x['marker']]+=int(x['observed_taxa'])
        if x['marker_review_status']!='no_current_copy_review_flag':
            flags[x['marker'],x['marker_review_status']]+=1
    access=json.loads((paths['accessibility']/'receipt.json').read_text())
    if (dict(counts)!=access['marker_row_counts'] or sum(counts.values())!=access['rows']
            or r['sites']!=len(actual) or r['markers']!=len(ready)
            or r['changed_markers']!=len(changed) or r['unchanged_markers']!=len(unchanged)):
        raise ValueError('Merged dimensions differ')
    verify()
    receipt=dict(status='passed_full_fcs_frame_source_row_readback',markers=len(ready),sites=len(actual),
                 changed_markers=len(changed),unchanged_markers=len(unchanged),
                 taxon_site_observations=sum(counts.values()),diagnostic_rows=len(expected_diagnostics),
                 flagged_markers=[dict(marker=m,status=s,sites=n) for (m,s),n in sorted(flags.items())],
                 frame_receipt_sha256=sha(result_path),plan_sha256=sha(a.plan),script_sha256=sha(__file__),
                 scope='Every merged site field, fit diagnostic and marker source reference matched to selected baseline or changed-marker rows. Observation totals checked against full retained accessibility receipt. Does not independently repeat source fits, accessibility calculation or FCS annotation.')
    a.output.write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
