#!/usr/bin/env python3
"""Close full primary geometry validation and reconcile all degenerate mappings."""
import csv,json,subprocess
from collections import Counter
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha


def main():
    paths=['metadata/primary_diagnostic_geometry_plan_20260927.json','metadata/primary_diagnostic_geometry_readback_plan_20260927.json','metadata/primary_diagnostic_geometry_launch_20260927.json','metadata/primary_diagnostic_geometry_readback_launch_20260927.json','metadata/primary_alignment_short_geometry_readback_20260927.json','metadata/primary_rmsd_diagnostic_and_short_completed_20260927.json']
    states={}
    for lp in paths[2:4]:
        launch=json.loads(Path(lp).read_text());unit=launch['unit']
        state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',unit,'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
        assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0'),state
        states[unit]=state
    pp,ap=map(Path,paths[:2]);plan=json.loads(pp.read_text());auditplan=json.loads(ap.read_text())
    assert auditplan['source_plan']==str(pp)
    for spec in [plan,auditplan]:
        for path,digest in spec['pins'].items():assert sha(path)==digest,path
    root=Path(plan['output']);rp=root/'receipt.json';r=json.loads(rp.read_text())
    proof=Path(auditplan['output']);a=json.loads(proof.read_text())
    assert r['status']=='complete_primary_diagnostic_geometry_pending_independent_readback' and r['plan_sha256']==sha(pp)
    assert a['status']=='passed_full_primary_diagnostic_geometry_readback' and a['plan_sha256']==sha(ap)
    assert a['producer_receipt_sha256']==sha(rp)
    assert a['alignments_checked']==r['alignments']==387646 and a['counts']==r['counts']
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    short=json.loads(Path(paths[4]).read_text());assert short['status']=='passed_all_short_primary_alignment_analytic_rmsd_checks'
    assert short['source_receipt_sha256']==r['diagnostic_receipt_sha256']
    key=lambda x:(x['pair_key'],x['mask'],int(x['order']))
    shortrows={key(x):x for x in short['rows']};degenerate={key(x):x for x in a['degenerate_alignments']}
    assert len(shortrows)==len(short['rows'])==len(degenerate)==len(a['degenerate_alignments'])==327
    assert set(shortrows)==set(degenerate)
    for k in shortrows:
        assert shortrows[k]['aligned_length']==degenerate[k]['aligned_length']
        assert shortrows[k]['rmsd_status']==degenerate[k]['rmsd_status']
    counts=Counter();rmsd=Counter();seen=set();short_seen=set()
    with (root/'alignment_geometry.tsv').open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            k=key(row);assert k not in seen;seen.add(k)
            counts[row['mask']+':'+row['geometry_status']]+=1;rmsd[row['rmsd_status']]+=1
            if row['geometry_status']=='degenerate_at_numeric_tolerance':
                short_seen.add(k);assert int(row['aligned_length'])<3
            else:assert row['geometry_status']=='unique_at_numeric_tolerance' and int(row['aligned_length'])>=3
    assert len(seen)==387646 and short_seen==set(shortrows) and dict(counts)==r['counts']
    assert sum(v for k,v in rmsd.items() if k!='within_printed_rounding')==27
    paths += [str(rp),str(root/'alignment_geometry.tsv'),str(proof)]
    result=dict(status='complete_verified_primary_diagnostic_geometry',alignments=387646,numerically_unique_rotations=387319,degenerate_short_alignments=327,rmsd_classification_counts=dict(rmsd),counts=dict(counts),terminal_states=states,source_hashes={p:sha(p) for p in paths},maximum_scaled_quaternion_curvature_error=a['maximum_scaled_quaternion_curvature_error'],near_zero_quaternion_gaps=a['near_zero_quaternion_gaps'],script_sha256=sha(__file__),scope='Full geometry production and serialized independent readback passed. All 327 degenerate mappings exactly match the separately checked short-alignment census; all 27 RMSD discrepancies remain retained. Unique rotation under numerical tolerance is not prediction accuracy, confidence/coverage qualification or biological inference. Original strict audit failure remains unchanged.')
    with Path('metadata/primary_diagnostic_geometry_completed_20260927.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
