#!/usr/bin/env python3
"""Summarize all alternative common-core contrasts without selecting favorable orders."""
import argparse,csv,itertools,json,math
from collections import Counter
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha


def direction(values,expected,margin):
    if len(values)!=expected:return 'incomplete'
    lo,hi=min(values),max(values)
    if lo>margin:return 'a_farther_from_reference'
    if hi < -margin:return 'b_farther_from_reference'
    if lo < -margin and hi>margin:return 'changes_direction_across_alternatives'
    return 'touches_or_within_margin_band'


def summarize(rows,screen,expected,margins):
    values=[float(r['rmsd_ar_minus_br']) for r in rows if r[screen+'_pass']=='1']
    if any(not math.isfinite(v) for v in values):raise ValueError('Nonfinite eligible contrast')
    result=dict(expected_alternatives=expected,eligible_alternatives=len(values),complete=int(len(values)==expected),
                eligible_contrast_min=min(values) if values else '',eligible_contrast_max=max(values) if values else '',
                eligible_contrast_range=max(values)-min(values) if values else '')
    for item in margins:result[item['id']]=direction(values,expected,item['angstrom'])
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();p=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        if sha(args.plan)!=ph:raise ValueError('Changed plan')
        for f,h in p['pins'].items():
            if sha(f)!=h:raise ValueError('Changed source: '+f)
    verify();root=Path(p['fits']);r=json.loads((root/'receipt.json').read_text());a=json.loads(Path(p['fit_readback']).read_text())
    if a['status']!='passed_full_common_residue_domain_fit_readback' or a['producer_receipt_sha256']!=sha(root/'receipt.json'):raise ValueError('Unbound numeric audit')
    if sha(root/'common_residue_fits.tsv')!=r['artifacts']['common_residue_fits.tsv']:raise ValueError('Changed fit table')
    out=Path(p['output']);out.mkdir(parents=True,exist_ok=False);seen=set();source_rows=0;output_rows=0;counts=Counter()
    keys=['triad_key','scope','mask','mapping_definition','screen','expected_alternatives','eligible_alternatives','complete','eligible_contrast_min','eligible_contrast_max','eligible_contrast_range']+[m['id'] for m in p['margins']]
    expected={(m,d,*o) for m in ['full','plddt70'] for d in ['reference_common','cycle_consistent'] for o in itertools.product(['0','1'],repeat=3)}
    with (root/'common_residue_fits.tsv').open() as f,(out/'robustness.tsv').open('w') as target:
        writer=csv.DictWriter(target,fieldnames=keys,delimiter='\t',lineterminator='\n');writer.writeheader()
        for tk,it in itertools.groupby(csv.DictReader(f,delimiter='\t'),key=lambda x:x['triad_key']):
            rows=list(it);source_rows+=len(rows)
            if tk in seen:raise ValueError('Noncontiguous/repeated triad')
            seen.add(tk)
            actual={(x['mask'],x['mapping_definition'],x['order_ab'],x['order_ar'],x['order_br']) for x in rows}
            if actual!=expected or len(rows)!=32:raise ValueError('Incomplete alternative grid')
            groups=[('all_alternatives','all','all',rows,32)]
            for mask in ['full','plddt70']:
                for definition in ['reference_common','cycle_consistent']:
                    groups.append(('input_orders',mask,definition,[x for x in rows if x['mask']==mask and x['mapping_definition']==definition],8))
            for scope,mask,definition,group,n in groups:
                for screen in p['screens']:
                    stats=summarize(group,screen,n,p['margins']);record=dict(triad_key=tk,scope=scope,mask=mask,mapping_definition=definition,screen=screen,**stats)
                    writer.writerow(record);output_rows+=1
                    for margin in p['margins']:counts['|'.join([scope,mask,definition,screen,margin['id'],stats[margin['id']]])]+=1
    if source_rows!=a['fit_rows_checked'] or output_rows!=len(seen)*30:raise ValueError('Incomplete summary cohort')
    verify()
    result=dict(status='complete_common_core_robustness_pending_readback',plan_sha256=ph,fit_receipt_sha256=sha(root/'receipt.json'),fit_readback_sha256=sha(p['fit_readback']),triads=len(seen),source_fit_rows=source_rows,summary_rows=output_rows,classification_counts=dict(counts),artifacts={'robustness.tsv':sha(out/'robustness.tsv')},scope='Exhaustive ranges across eight orders within each mask/mapping definition and all 32 alternatives jointly. All alternatives must pass a screen for a direction label; incomplete sets remain explicit. Zero/0.01/0.1 Angstrom margins are descriptive sensitivities, not error calibration or biological cutoffs. Domains, annotation boundaries/policies, references and guides still require linked sensitivity analysis; no independent observations or biological asymmetry claim.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='classification_counts'},indent=2))


if __name__=='__main__':main()
