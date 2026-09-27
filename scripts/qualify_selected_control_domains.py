#!/usr/bin/env python3
"""Retain every matched domain configuration across masks, boundaries and six screens."""
import argparse,csv,gzip,json
from pathlib import Path
from collections import Counter
from screen_duplication_domain_alignment_coverage import sha


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify();inventory=Path(plan['inventory']);qualified=Path(plan['qualification'])
    receipts={}
    for name,root in [('inventory',inventory),('qualification',qualified)]:
        r=json.loads((root/'receipt.json').read_text());audit=json.loads(Path(plan[name+'_audit']).read_text());assert audit['status']==plan[name+'_audit_status'] and audit['source_receipt_sha256']==sha(root/'receipt.json')
        for f,h in r['artifacts'].items():assert sha(root/f)==h
        receipts[name]=sha(root/'receipt.json')
    flags={}
    for cohort in ['target','background']:
        flags[cohort]={}
        for row in csv.DictReader((qualified/(cohort+'.tsv')).open(),delimiter='\t'):
            key=row['pair_key'],row['mask'];assert key not in flags[cohort]
            flags[cohort][key]={s:row[s+'_pass']=='1' for s in plan['screens']}
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);seen=set();count=0;counts=Counter();eligible_counts=Counter()
    fields=['domain_config_id','boundary','mask_cohort','screen','configuration_status','shared_domain_count','eligible_domain_count','eligible_pfams','ineligible_pfams','eligibility_status']
    with gzip.open(out/'configuration_eligibility.tsv.gz','wt') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for line in (inventory/'domain_configurations.jsonl').open():
            c=json.loads(line);cid=c['domain_config_id'];assert cid not in seen;seen.add(cid)
            for boundary in ['alignment','envelope']:
                matches=[m for m in c['matches'] if m['boundary']==boundary];accessions={m['pfam_accession'] for m in matches}
                assert len(accessions)==len(matches) and accessions==set(c['shared_accessions'])
                for mask_cohort in ['full','plddt70','both']:
                    masks=['full','plddt70'] if mask_cohort=='both' else [mask_cohort]
                    for screen in plan['screens']:
                        eligible=[]
                        for m in matches:
                            passes=all(flags[side][m[side+'_domain_pair'],mask][screen] for side in ['target','background'] for mask in masks)
                            if passes and c['status']=='shared_domain_comparisons':eligible.append(m['pfam_accession'])
                        eligible=sorted(eligible);ineligible=sorted(accessions-set(eligible))
                        status=c['status'] if c['status']!='shared_domain_comparisons' else 'all_shared_domains_pass' if len(eligible)==len(accessions) else 'some_shared_domains_pass' if eligible else 'no_shared_domains_pass'
                        w.writerow(dict(domain_config_id=cid,boundary=boundary,mask_cohort=mask_cohort,screen=screen,configuration_status=c['status'],shared_domain_count=len(accessions),eligible_domain_count=len(eligible),eligible_pfams=json.dumps(eligible,separators=(',',':')),ineligible_pfams=json.dumps(ineligible,separators=(',',':')),eligibility_status=status))
                        counts[boundary+'|'+mask_cohort+'|'+screen+'|'+status]+=1
                        eligible_counts[boundary+'|'+mask_cohort+'|'+screen]+=len(eligible);count+=1
            if len(seen)%20000==0:print('Qualified configurations',len(seen),flush=True)
    assert len(seen)==json.loads((inventory/'receipt.json').read_text())['distinct_configurations']
    assert count==len(seen)*2*3*len(plan['screens']);verify()
    result=dict(status='complete_selected_control_domain_qualification_pending_readback',plan_sha256=ph,configurations=len(seen),rows=count,status_counts=dict(counts),eligible_domain_occurrences=dict(eligible_counts),source_receipts=receipts,artifacts={'configuration_eligibility.tsv.gz':sha(out/'configuration_eligibility.tsv.gz')},scope='Every distinct matched protein-pair/policy configuration, both domain boundaries, six thresholds and full/pLDDT70/common-mask cohorts. Shared accessions must pass for both target and background and both alignment orders. Identical models and empty cohorts retained explicitly. Counts are repeated configuration/domain occurrences, not independent events or duplication effects.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['status_counts','eligible_domain_occurrences']},indent=2))

if __name__=='__main__':main()
