#!/usr/bin/env python3
"""Require every reference/policy/boundary and fit alternative for event/domain robustness."""
import argparse,csv,json,sqlite3
from collections import defaultdict,Counter
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha
from summarize_common_core_robustness import direction

EVENT=['guide','family','gene_node','gene_a','gene_b']
SEQ=['taxon_id','chosen_reference_gene','duplicate_pair_sequence_distance','sequence_signed_difference','sequence_normalized_contrast','sequence_direction','sequence_covariate_status']


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();p=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        if sha(args.plan)!=ph:raise ValueError('Changed plan')
        for path,h in p['pins'].items():
            if sha(path)!=h:raise ValueError('Changed source: '+path)
    verify();rob=Path(p['robustness']);rr=json.loads((rob/'receipt.json').read_text());audit=json.loads(Path(p['robustness_readback']).read_text())
    if audit['status']!='passed_full_common_core_robustness_readback' or audit['producer_receipt_sha256']!=sha(rob/'receipt.json'):raise ValueError('Unbound robustness audit')
    if sha(rob/'robustness.tsv')!=rr['artifacts']['robustness.tsv']:raise ValueError('Changed robustness table')
    integration=Path(p['integration']);ir=json.loads((integration/'receipt.json').read_text());ia=json.loads(Path(p['integration_readback']).read_text())
    if ia['status']!='passed_full_domain_triad_integration_readback' or ia['producer_receipt_sha256']!=sha(integration/'receipt.json') or sha(integration/'domain_triads.sqlite')!=ir['artifacts']['domain_triads.sqlite']:raise ValueError('Unbound integration')
    maps=Path(p['mappings']);mr=json.loads((maps/'receipt.json').read_text())
    if sha(maps/'event_domain_links.tsv')!=mr['artifacts']['event_domain_links.tsv']:raise ValueError('Changed event links')
    db=sqlite3.connect(f'file:{integration}/domain_triads.sqlite?mode=ro',uri=True);db.row_factory=sqlite3.Row
    units=defaultdict(set);sequences={}
    for x in db.execute('SELECT * FROM triads'):
        event=tuple(x[f] for f in EVENT)
        for boundary in ['alignment','envelope']:units[event].add((x['reference_gene'],x['policy'],boundary))
    for x in db.execute('SELECT * FROM sequence_covariates'):
        event=tuple(x[f] for f in EVENT)
        if event in sequences:raise ValueError('Repeated sequence covariate')
        sequences[event]={f:x[f] for f in SEQ}
    db.close()
    if set(units)!=set(sequences):raise ValueError('Incomplete event covariates')
    summaries={}
    for x in csv.DictReader((rob/'robustness.tsv').open(),delimiter='\t'):
        if x['scope']=='all_alternatives':
            key=x['triad_key'],x['screen']
            if key in summaries:raise ValueError('Repeated interval summary')
            summaries[key]=x
    groups=defaultdict(dict)
    for x in csv.DictReader((maps/'event_domain_links.tsv').open(),delimiter='\t'):
        event=tuple(x[f] for f in EVENT);unit=x['reference_gene'],x['policy'],x['boundary'];key=event+(x['pfam_accession'],)
        if unit not in units[event] or unit in groups[key]:raise ValueError('Repeated/unexpected event/domain alternative')
        groups[key][unit]=x['triad_key']
    out=Path(p['output']);out.mkdir(parents=True,exist_ok=False);counts=Counter();rows=0
    fields=EVENT+['pfam_accession','screen','reference_choices','expected_reference_policy_boundary_units','observed_units','missing_units','complete_units','complete','eligible_contrast_min','eligible_contrast_max']+SEQ+[m['id'] for m in p['margins']]
    with (out/'event_domain_robustness.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for key,links in sorted(groups.items()):
            event=key[:-1];expected=len(units[event]);missing=expected-len(links)
            for screen in p['screens']:
                source=[summaries[tk,screen] for tk in links.values()]
                complete_units=sum(x['complete']=='1' for x in source);complete=complete_units==expected
                lows=[float(x['eligible_contrast_min']) for x in source if x['eligible_contrast_min']!=''];highs=[float(x['eligible_contrast_max']) for x in source if x['eligible_contrast_max']!='']
                lo=min(lows) if lows else '';hi=max(highs) if highs else ''
                r=dict(zip(EVENT,event));r.update(pfam_accession=key[-1],screen=screen,reference_choices=len({x[0] for x in units[event]}),expected_reference_policy_boundary_units=expected,observed_units=len(links),missing_units=missing,complete_units=complete_units,complete=int(complete),eligible_contrast_min=lo,eligible_contrast_max=hi,**sequences[event])
                for margin in p['margins']:
                    label=direction([lo,hi],2,margin['angstrom']) if complete else 'incomplete'
                    r[margin['id']]=label;counts['|'.join([event[0],screen,margin['id'],label])]+=1
                w.writerow(r);rows+=1
    event_domains=Counter(key[:-1] for key in groups)
    with (out/'event_availability.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=EVENT+['reference_choices','domain_candidates','availability']+SEQ,delimiter='\t',lineterminator='\n');w.writeheader()
        for event,u in sorted(units.items()):
            w.writerow(dict(zip(EVENT,event))|dict(reference_choices=len({x[0] for x in u}),domain_candidates=event_domains[event],availability='domain_candidates_present' if event_domains[event] else 'no_common_domain_candidate',**sequences[event]))
    verify()
    result=dict(status='complete_domain_event_robustness_pending_readback',plan_sha256=ph,event_universe=len(units),event_domain_combinations=len(groups),event_domain_screen_rows=rows,events_with_candidates=sum(v>0 for v in event_domains.values()),classification_counts=dict(counts),artifacts={f.name:sha(f) for f in out.iterdir() if f.is_file()},scope='Every provisionally referenced guide/event retained. Each Pfam separate; all tied references, four policies, both boundaries and all 32 fit alternatives required for complete direction label. Missing candidate units explicit. Sequence covariates remain whole-protein and gene-oriented; no ancestral direction, independence or significance claim.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='classification_counts'},indent=2))


if __name__=='__main__':main()
