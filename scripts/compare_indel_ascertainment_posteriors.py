#!/usr/bin/env python3
"""Compare matched ancestral gap probabilities under two ascertainment assumptions."""
import csv,gzip,json,subprocess,time
from pathlib import Path
import numpy as np
import psutil
from prepare_case_ancestral_neighborhoods import sha


def main():
    pp=Path('metadata/indel_ascertainment_comparison_plan_20260927.json');plan=json.loads(pp.read_text())
    for p,h in plan['pins'].items():assert sha(p)==h,p
    identity=json.loads(Path(plan['auditor_launch']).read_text())
    while True:
        state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',identity['unit'],'--property=ActiveState,MainPID,Result,ExecMainStatus'],text=True).splitlines())
        if state['ActiveState']=='inactive':
            assert state['Result']=='success' and state['ExecMainStatus']=='0';break
        assert state['ActiveState']=='active',state
        p=psutil.Process(identity['pid']);assert p.create_time()==identity['created'] and p.cmdline()==identity['cmdline'] and int(state['MainPID'])==p.pid
        print('waiting_for_verified_complete_posterior_audit',flush=True);time.sleep(30)
    auditroot=Path(plan['audit_output']);ar=json.loads((auditroot/'receipt.json').read_text())
    assert ar['status']=='all_selected_fit_gap_probabilities_independently_verified' and ar['input_dispositions']==312
    root=Path(plan['posterior_output']);sr=json.loads((root/'receipt.json').read_text());assert sha(root/'receipt.json')==ar['source_receipt_sha256']
    with Path(plan['mapping']).open() as h:mappings=list(csv.DictReader(h,delimiter='\t'))
    models=json.loads(Path(plan['fit_plan']).read_text())['jobs'];by={j['job_id']:j for j in models}
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);summaries=[];count=0;switches=0;strong=0
    with gzip.open(out/'character_comparisons.tsv.gz','wt') as h:
        writer=csv.writer(h,delimiter='\t',lineterminator='\n');writer.writerow(['encoding','level','character','all_taxa_gap_probability','observed_mask_gap_probability','absolute_difference','map_switch','opposed_states_both_at_least_090'])
        for jid,job in by.items():
            if job['correction']!='all_taxa':continue
            other=jid.removesuffix('-all_taxa')+'-observed_mask';assert other in by
            assert job['tree']==by[other]['tree'] and job['characters']==by[other]['characters']
            receipts=[]
            for name in [jid,other]:
                rp=root/name/'receipt.json';assert sha(rp)==sr['job_receipts'][name+'/receipt.json']
                r=json.loads(rp.read_text());ap=auditroot/(name+'.json');assert sha(ap)==ar['artifacts'][ap.name]
                assert json.loads(ap.read_text())['source_receipt_sha256']==sha(rp)
                for f,digest in r['artifacts'].items():assert sha(root/name/f)==digest
                receipts.append(r)
            encoding=jid.removesuffix('-all_taxa')
            if not job['character_count']:
                summaries.append(dict(encoding=encoding,status='no_coded_characters',node_characters=0,maximum_difference=0,map_switches=0,strong_conflicts=0));continue
            selected=sorted([r for r in mappings if r['job_id']==jid and r['guide']=='profile' and r['status']=='matched_unrooted_vertex'],key=lambda r:int(r['level']))
            assert len(selected)==3
            with np.load(root/jid/'posterior.npz',allow_pickle=False) as a,np.load(root/other/'posterior.npz',allow_pickle=False) as b:
                assert np.array_equal(a['character'],b['character']) and np.array_equal(a['node_names'],b['node_names'])
                names=list(a['node_names']);differences=[];sw=0;st=0
                for row in selected:
                    counterpart=[r for r in mappings if r['job_id']==other and r['guide']=='profile' and r['level']==row['level']];assert len(counterpart)==1
                    assert counterpart[0]['posterior_node']==row['posterior_node'] and counterpart[0]['incident_tip_partition_sha256']==row['incident_tip_partition_sha256']
                    index=names.index(row['posterior_node']);x=a['probability_gap'][index];y=b['probability_gap'][index];d=abs(x-y);flips=(x>=.5)!=(y>=.5);conflicts=((x>=.9)&(y<=.1))|((y>=.9)&(x<=.1))
                    for col in range(len(x)):writer.writerow([encoding,row['level'],col+1,float(x[col]),float(y[col]),float(d[col]),int(flips[col]),int(conflicts[col])])
                    differences.extend(d.tolist());sw+=int(flips.sum());st+=int(conflicts.sum())
                n=len(differences);count+=n;switches+=sw;strong+=st
                summaries.append(dict(encoding=encoding,status='matched_candidate_nodes_compared',node_characters=n,maximum_difference=max(differences),map_switches=sw,strong_conflicts=st))
    assert len(summaries)==156
    with (out/'encoding_summary.tsv').open('w') as h:w=csv.DictWriter(h,list(summaries[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(summaries)
    receipt=dict(status='complete_156_ascertainment_probability_comparisons',encodings=156,node_characters=count,map_switches=switches,strong_conflicts=strong,maximum_probability_difference=max(r['maximum_difference'] for r in summaries),plan_sha256=sha(pp),audit_receipt_sha256=sha(auditroot/'receipt.json'),auditor_terminal_state=state,artifacts={p.name:sha(p) for p in out.iterdir()},scope=plan['scope'])
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':main()
