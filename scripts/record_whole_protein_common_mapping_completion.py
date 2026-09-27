#!/usr/bin/env python3
"""Bind successful full triad mapping/readback without promoting maps to biological inference."""
import json,subprocess
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha


def main():
    pp=Path('metadata/whole_protein_common_residues_plan_20260927.json');plan=json.loads(pp.read_text());states={}
    for unit in ['fungal-whole-protein-common-residues-20260927.service','fungal-whole-protein-common-residues-readback-20260927.service']:
        state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',unit,'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
        assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0'),state;states[unit]=state
    for path,digest in plan['pins'].items():assert sha(path)==digest,path
    root=Path(plan['output']);rp=root/'receipt.json';r=json.loads(rp.read_text());ap=Path('metadata/whole_protein_common_residues_readback_20260927.json');a=json.loads(ap.read_text())
    assert r['status']=='complete_whole_protein_common_residue_inventory_pending_readback' and a['status']=='passed_full_whole_protein_common_residue_readback'
    assert r['plan_sha256']==a['plan_sha256']==sha(pp) and a['producer_receipt_sha256']==sha(rp)
    for name,digest in r['artifacts'].items():assert sha(root/name)==digest
    fields=['distinct_oriented_model_triads','event_reference_links','mask_order_dispositions','common_reference_residue_occurrences','cycle_consistent_residue_occurrences','counts']
    for field in fields:assert r[field]==a[field],field
    assert r['distinct_oriented_model_triads']==17619 and r['event_reference_links']==36944 and r['mask_order_dispositions']==281904
    assert sum(r['counts'].values())==281904
    result=dict(status='complete_verified_whole_protein_common_residue_mapping',terminal_states=states,source_hashes={str(p):sha(p) for p in [pp,rp,ap]},artifacts=r['artifacts'],**{field:a[field] for field in fields},dispositions_with_nonempty_reference_intersection=a['dispositions_with_nonempty_reference_intersection'],dispositions_with_nonempty_consistent_mapping=a['dispositions_with_nonempty_consistent_mapping'],dispositions_with_any_mapping_disagreement=a['dispositions_with_any_mapping_disagreement'],script_sha256=sha(__file__),scope='Complete original-position mapping and independent full record reconstruction verified, including all missing/numerical exclusions and shared-model identities. Residue occurrences repeat across masks/orders/events. Common-core fitting, coverage/geometry qualification and biological inference remain separate stages.')
    with Path('metadata/whole_protein_common_mapping_completed_20260927.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['counts','artifacts','source_hashes']},indent=2))


if __name__=='__main__':main()
