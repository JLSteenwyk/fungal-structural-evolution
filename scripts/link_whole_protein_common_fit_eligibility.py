#!/usr/bin/env python3
"""Link all whole-protein triads to duplication events without selecting reference ties."""
import csv
import json
from pathlib import Path
import pandas as pd
from screen_duplication_domain_alignment_coverage import sha


def main():
    completion = Path('metadata/whole_protein_common_fits_completed_20260927.json')
    c = json.loads(completion.read_text())
    assert c['status'] == 'complete_verified_whole_protein_common_fits_and_coverage_summary'
    for path, digest in c['source_hashes'].items():
        assert sha(path) == digest
    pp = Path('metadata/whole_protein_common_fits_plan_20260927.json')
    plan = json.loads(pp.read_text()); root = Path(plan['output'])
    receipt = json.loads((root/'receipt.json').read_text())
    table = root/'common_residue_fits.tsv'
    assert sha(table) == receipt['artifacts'][table.name]
    mp = json.loads(Path(plan['mapping_plan']).read_text()); mr = Path(mp['output'])
    mapping_receipt = json.loads((mr/'receipt.json').read_text())
    assert sha(mr/'receipt.json') == receipt['mapping_receipt_sha256']
    ep = mr/'event_reference_triads.tsv'
    assert sha(ep) == mapping_receipt['artifacts'][ep.name]
    screens = [s['id'] for s in plan['screens']]
    frame = pd.read_csv(table, sep='\t', usecols=['triad_id','distinct_models']+[s+'_pass' for s in screens])
    groups = frame.groupby('triad_id', sort=True)
    assert len(groups) == 17619 and (groups.size() == 32).all()
    assert (groups.distinct_models.nunique() == 1).all()
    triads = groups.distinct_models.first().to_frame()
    for screen in screens:
        triads[screen+'_passing_configurations'] = groups[screen+'_pass'].sum()
        triads[screen+'_all_configurations_pass'] = (triads[screen+'_passing_configurations'] == 32).astype(int)
    # Check every grouped count against the full serialized input, independently.
    scalar = {}
    with table.open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            values = scalar.setdefault(row['triad_id'], [0]*len(screens))
            for i, screen in enumerate(screens):
                values[i] += int(row[screen+'_pass'])
    for key, values in scalar.items():
        assert values == [int(triads.loc[key, s+'_passing_configurations']) for s in screens]
    events = pd.read_csv(ep, sep='\t', dtype=str, keep_default_na=False)
    assert len(events) == 36944 and not events.duplicated().any()
    linked = events.merge(triads.reset_index(), how='left', on='triad_id', validate='many_to_one', sort=False)
    assert len(linked) == len(events) and linked.notna().all().all()
    assert linked[events.columns].equals(events)
    summaries = []
    event_key = ['guide','family','gene_node','gene_a','gene_b']
    for guide, sub in linked.groupby('guide', sort=True):
        for screen in screens:
            eligible = sub[(sub.distinct_models == 3) & (sub[screen+'_all_configurations_pass'] == 1)]
            summaries.append(dict(guide=guide, screen=screen, event_reference_links=len(sub),
                                  unique_events=len(sub[event_key].drop_duplicates()),
                                  three_model_robust_reference_links=len(eligible),
                                  events_with_at_least_one_robust_reference=len(eligible[event_key].drop_duplicates())))
    out = Path('results/structural_comparisons/whole-protein-common-fit-event-links-20260927-v1')
    out.mkdir(exist_ok=False)
    artifacts = {}
    for name, data in [('triad_eligibility.tsv',triads.reset_index()),('event_reference_eligibility.tsv',linked),('eligibility_summary.tsv',pd.DataFrame(summaries))]:
        path = out/name; data.to_csv(path, sep='\t', index=False)
        observed = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
        assert observed.equals(data.reset_index(drop=True).astype(str))
        artifacts[name] = sha(path)
    assert sha(table) == receipt['artifacts'][table.name] and sha(ep) == mapping_receipt['artifacts'][ep.name]
    result = dict(status='complete_whole_protein_common_fit_event_eligibility',
                  source_hashes={str(p):sha(p) for p in [completion,pp,root/'receipt.json',mr/'receipt.json',ep]},
                  script_sha256=sha(__file__),triads=len(triads),event_reference_links=len(linked),
                  summaries=summaries,artifacts=artifacts,
                  scope='Descriptive eligibility for all 32 mask/order/mapping configurations and six screens. All reference ties and shared-model cases retained in linked table; robust-reference counts require three distinct models. Events with one passing reference need not be robust to all references. No hypothesis test or evolutionary asymmetry claim.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
