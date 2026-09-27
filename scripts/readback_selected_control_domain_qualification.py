#!/usr/bin/env python3
"""Check the full domain-qualification grid using independently intersected pass sets."""
import argparse,csv,gzip,hashlib,json,time,subprocess,itertools
from pathlib import Path
from collections import Counter,defaultdict
import psutil


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--launch',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    plan=json.loads(a.plan.read_text());ph=sha(a.plan);dep=json.loads(a.launch.read_text());lh=sha(a.launch);assert dep['plan_sha256']==ph
    while True:
        try:
            p=psutil.Process(dep['pid'])
            if p.create_time()!=dep['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==dep['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',dep['unit'],'-p','ActiveState','-p','ExecMainStatus'],text=True).splitlines());assert state=={'ActiveState':'inactive','ExecMainStatus':'0'}
    assert sha(a.plan)==ph and sha(a.launch)==lh
    for p,h in plan['pins'].items():assert sha(p)==h,p
    root=Path(plan['output']);r=json.loads((root/'receipt.json').read_text());assert r['status']=='complete_selected_control_domain_qualification_pending_readback' and r['plan_sha256']==ph
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    for name in ['inventory','qualification']:assert sha(Path(plan[name])/'receipt.json')==r['source_receipts'][name]
    passing=defaultdict(set);universe={}
    for side in ['target','background']:
        seen=set()
        for row in csv.DictReader((Path(plan['qualification'])/(side+'.tsv')).open(),delimiter='\t'):
            key=row['pair_key'],row['mask'];assert key not in seen;seen.add(key)
            for screen in plan['screens']:
                if row[screen+'_pass']=='1':passing[side,row['mask'],screen].add(row['pair_key'])
        universe[side]=seen
    configs={};expected={}
    for line in (Path(plan['inventory'])/'domain_configurations.jsonl').open():
        c=json.loads(line);cid=c['domain_config_id'];assert cid not in configs;configs[cid]=c
        for boundary in ['alignment','envelope']:
            match=[m for m in c['matches'] if m['boundary']==boundary]
            assert {m['pfam_accession'] for m in match}==set(c['shared_accessions'])
            for screen in plan['screens']:
                eligible={}
                for mask in ['full','plddt70']:
                    sets=[]
                    for side in ['target','background']:
                        assert all((m[side+'_domain_pair'],mask) in universe[side] for m in match)
                        sets.append({m['pfam_accession'] for m in match if m[side+'_domain_pair'] in passing[side,mask,screen]})
                    eligible[mask]=sets[0].intersection(sets[1]) if c['status']=='shared_domain_comparisons' else set()
                eligible['both']=eligible['full'].intersection(eligible['plddt70'])
                for mask,values in eligible.items():expected[cid,boundary,mask,screen]=tuple(sorted(values))
    seen=set();counts=Counter();occurrences=Counter()
    with gzip.open(root/'configuration_eligibility.tsv.gz','rt') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            key=tuple(row[k] for k in ['domain_config_id','boundary','mask_cohort','screen']);assert key in expected and key not in seen;seen.add(key)
            cid,boundary,mask,screen=key;c=configs[cid];eligible=list(expected[key]);missing=sorted(set(c['shared_accessions'])-set(eligible))
            assert json.loads(row['eligible_pfams'])==eligible and json.loads(row['ineligible_pfams'])==missing
            assert int(row['eligible_domain_count'])==len(eligible) and int(row['shared_domain_count'])==len(c['shared_accessions'])
            assert row['configuration_status']==c['status']
            status=c['status'] if c['status']!='shared_domain_comparisons' else 'no_shared_domains_pass' if not eligible else 'some_shared_domains_pass' if missing else 'all_shared_domains_pass'
            assert row['eligibility_status']==status
            counts[boundary+'|'+mask+'|'+screen+'|'+status]+=1;occurrences[boundary+'|'+mask+'|'+screen]+=len(eligible)
    assert seen==set(expected) and len(seen)==r['rows']==len(configs)*36 and len(configs)==r['configurations']
    assert dict(counts)==r['status_counts'] and dict(occurrences)==r['eligible_domain_occurrences']
    for p,h in plan['pins'].items():assert sha(p)==h,p
    result=dict(status='passed_full_selected_control_domain_qualification_readback',rows=len(seen),configurations=len(configs),source_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(__file__),status_counts=dict(counts),eligible_domain_occurrences=dict(occurrences),scope='All configuration/boundary/mask/screen cells, exact eligible/ineligible accession sets, identical-model exclusions and status counts independently reconstructed by target/background pass-set intersections; both-mask cohort uses accession intersection. Not independent events or effect estimates.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['status_counts','eligible_domain_occurrences']},indent=2),flush=True)

if __name__=='__main__':main()
