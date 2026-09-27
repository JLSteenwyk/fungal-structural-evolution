#!/usr/bin/env python3
"""Bind full wood-decay bootstrap audits and preserve their limited inferential scope."""
import json,subprocess
from pathlib import Path
import pandas as pd
from screen_duplication_domain_alignment_coverage import sha


def main():
    states={}
    for name in ['bootstrap','bootstrap-readback','bootstrap-summary','bootstrap-summary-readback']:
        unit='fungal-wood-decay-'+name+'-20260927.service'
        state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',unit,'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
        assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0');states[unit]=state
    base=Path('results/ecology');summary=base/'wood-decay-bootstrap-summary-20260927-v1';rp=summary/'receipt.json';r=json.loads(rp.read_text())
    ap=Path('metadata/wood_decay_bootstrap_summary_readback_20260927.json');a=json.loads(ap.read_text())
    assert a['status']=='passed_full_wood_decay_bootstrap_summary_readback' and a['source_receipt_sha256']==sha(rp)
    for k in ['split_rows','distribution_rows','paired_coding_rows','coding_tree_combinations']:assert a[k]==r[k]
    for path,digest in r['sources'].items():assert sha(path)==digest
    for name,digest in r['artifacts'].items():assert sha(summary/name)==digest
    auditroot=base/'wood-decay-bootstrap-readback-20260927-v1';arp=auditroot/'receipt.json';ar=json.loads(arp.read_text())
    assert ar['status']=='passed_full_bootstrap_wood_decay_edge_network_flow_readback'
    assert ar['producer_receipt_sha256']==sha(base/'wood-decay-bootstrap-edges-20260927-v1/receipt.json')
    assert ar['trees']==r['trees']==2000 and ar['edge_rows']==10490000 and ar['constrained_costs']==41960000
    assert len(ar['proofs'])==2000
    for name,digest in ar['proofs'].items():assert sha(auditroot/name)==digest
    d=pd.read_csv(summary/'tree_metric_distributions.tsv',sep='\t');primary=d[d.coding.eq('uncertain_unknown')]
    assert primary[primary.metric.eq('minimum_changes')].value.eq(2).all()
    assert primary[primary.metric.eq('required_change')].value.eq(0).all()
    splits=pd.read_csv(summary/'split_uncertainty.tsv',sep='\t');required=splits[splits.required_change.gt(0)]
    assert len(required)==4 and set(required.coding)=={'jaapia_1_botryobasidium_0','jaapia_1_botryobasidium_1'}
    assert required.required_change.eq(996).all() and required.present.eq(1000).all()
    assert required.canonical_side.eq('F104355;F202697;F38799;F5364').all()
    result=dict(status='complete_verified_wood_decay_bootstrap_uncertainty',terminal_states=states,source_hashes={str(p):sha(p) for p in [rp,ap,arp]},script_sha256=sha(__file__),trees=2000,coding_tree_combinations=10000,independently_checked_edges=10490000,independently_checked_endpoint_costs=41960000,summary_artifacts=r['artifacts'],primary_required_changes_per_tree=0,primary_minimum_changes_per_tree=2,required_split_under_jaapia_brown=required.to_dict('records'),scope='Primary uncertain labels remain unknown: no edge is required to change in any primary bootstrap optimum set. A required split appears in 996/1000 trees per ensemble only when Jaapia is forced brown. Frequencies are topology/coding sensitivity, not posterior origins, independent ecological transitions or structural effects.')
    p=Path('metadata/wood_decay_bootstrap_completed_20260927.json')
    with p.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='required_split_under_jaapia_brown'},indent=2))


if __name__=='__main__':main()
