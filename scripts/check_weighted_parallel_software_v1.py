#!/usr/bin/env python3
"""Compare full serial/parallel grids and independently reject private corruptions."""
import argparse
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
import shutil

from ancestral_chain_attempt import sha
from full_weighted_covariance_qualification import run as serial, SUMMARY
from full_weighted_covariance_qualification_parallel_v1 import run as parallel
from reference_measurement_union_sources import bind, verify


def rejected(action):
    try: action()
    except (AssertionError, ValueError, KeyError, FileExistsError, StopIteration):return
    raise AssertionError('Private corruption or completed-stage restart was accepted')


def audits(root,entry):
    with gzip.open(root/entry['audit_file'],'rt') as f:
        return [json.loads(line) for line in f]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fixtures',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--receipt',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(exist_ok=False);assert not a.receipt.exists()
    prepared=json.loads(a.fixtures.read_text());verify(prepared['source_hashes'])
    assert prepared['status']=='prepared_three_complete_private_parallel_software_grids_v1'
    pins=dict(prepared['source_hashes']);bind(pins,a.fixtures);bind(pins,Path(__file__))
    compared=linked=0;grids=[];dimensions=set();policies=set();workers=set()
    for grid in prepared['grids']:
        s=serial(grid['serial']); sr=serial(grid['serial'],True)
        q=parallel(grid['parallel']); qr=parallel(grid['parallel'],True)
        assert all(s[k]==sr[k]==q[k]==qr[k] for k in SUMMARY)
        assert (q['cohorts'],q['designs'],q['settings'],q['numerical_audit_rows'],q['setting_audit_links'])==(5,150,600,6000,24000)
        roots={key:Path(json.loads(Path(grid[key]).read_text())['output']) for key in ['serial','parallel']}
        sm=json.loads((roots['serial']/'cohort_manifest.json').read_text())
        qm=json.loads((roots['parallel']/'cohort_manifest.json').read_text())
        group_workers=set()
        for se,qe in zip(sm,qm):
            assert se['cohort_id']==qe['cohort_id']
            left=audits(roots['serial'],se);right=audits(roots['parallel'],qe)
            identities={}
            for l,r in zip(left,right):
                identities[l['audit_id']]=r['audit_id']
                assert {k:v for k,v in l.items() if k not in ['audit_id','source_contract']}=={k:v for k,v in r.items() if k not in ['audit_id','source_contract']}
                if grid['qualified'] and r['disposition']=='numerically_qualified_exact_four_control_covariance_basis':
                    dimensions.add(len(r['retained_kernel_names']))
                    if not r['residual_diagonal_is_exact_uniform_one']:policies.add(r['control_policy'])
                compared+=1
            assert len(left)==len(right)==1200
            with gzip.open(roots['serial']/se['link_file'],'rt') as sf,gzip.open(roots['parallel']/qe['link_file'],'rt') as qf:
                header=sf.readline().rstrip('\n').split('\t');assert qf.readline().rstrip('\n').split('\t')==header
                index=header.index('audit_id')
                for l,r in zip(sf,qf):
                    ls=l.rstrip('\n').split('\t');rs=r.rstrip('\n').split('\t')
                    ls[index]=identities[ls[index]];assert ls==rs;linked+=1
                assert next(sf,None) is None and next(qf,None) is None
            for cp in (roots['parallel']/'checkpoints').glob('*.producer.json'):
                c=json.loads(cp.read_text());assert c['cached_numeric_inputs_preserved']
                assert c['worker_address_space_limit_bytes']==12*2**30
                group_workers.add((c['worker']['pid'],c['worker']['created']))
        assert len(group_workers)==2;workers.update(group_workers)
        for receipt in [sr,qr]:
            verify(receipt['source_hashes'])
            for path,digest in receipt['source_hashes'].items():bind(pins,path,digest)
        for root in roots.values():
            for path in [root/'receipt.json',root/'readback.json']:bind(pins,path)
        rejected(lambda:parallel(grid['parallel']));rejected(lambda:parallel(grid['parallel'],True))
        grids.append(dict(name=grid['name'],audits=6000,links=24000,distinct_producer_workers=len(group_workers),parallel_readback=str(roots['parallel']/'readback.json')))
    assert (compared,linked)==(18000,72000) and dimensions=={4,5,6}
    assert policies=={'background_node','background_pair','family_component'}
    rejected_cases=[];captures=0
    for item in prepared['negatives']:
        case=item['case'];plan=item['plan'];parallel(plan)
        root=Path(json.loads(Path(plan).read_text())['output'])
        rp=root/'receipt.json';mp=root/'cohort_manifest.json';m=json.loads(mp.read_text())
        entry=next(e for e in m if any(r['numerical_audit'] is not None for r in audits(root,e)))
        ap=root/entry['audit_file'];lp=root/entry['link_file']
        checkpoint=next(cp for cp in (root/'checkpoints').glob('*.producer.json') if json.loads(cp.read_text())['cohort_id']==entry['cohort_id'])
        archive=root/'original-private-serialized-fixture';archive.mkdir()
        for path in [rp,mp,ap,lp,checkpoint]:shutil.copyfile(path,archive/path.name)
        if case in ['foreign_link','missing_link','duplicate_link','changed_original_setting','changed_setting_ordinal']:
            lines=gzip.decompress(lp.read_bytes()).decode().splitlines();fields=lines[0].split('\t')
            if case=='missing_link':lines.pop()
            elif case=='duplicate_link':lines[-1]=lines[1]
            else:
                cells=lines[1].split('\t');field='audit_id' if case=='foreign_link' else 'scenario_id' if case=='changed_original_setting' else 'source_setting_ordinal'
                cells[fields.index(field)]='foreign' if field!='source_setting_ordinal' else str(int(cells[fields.index(field)])+1)
                lines[1]='\t'.join(cells)
            lp.write_bytes(gzip.compress(('\n'.join(lines)+'\n').encode(),mtime=0));entry['link_sha256']=sha(lp)
        else:
            records=audits(root,entry);r=next(r for r in records if r['numerical_audit'] is not None);value=r['numerical_audit']
            if case in ['raw_gram','projected_gram']:value[case][0][0]+=.01
            elif case=='raw_envelope':value['raw_roundoff_envelope'][0][0]*=1.02
            elif case=='projected_envelope':value['projected_roundoff_envelope'][0][0]*=1.02
            elif case=='diagnostic_rank':value['raw_diagnostics']['rank']+=1
            elif case=='retained_names':r['retained_kernel_names']=list(reversed(r['retained_kernel_names']))
            elif case=='uniform_class':r['residual_diagonal_is_exact_uniform_one']=not r['residual_diagonal_is_exact_uniform_one']
            elif case=='control_label':r['control_policy']='foreign'
            elif case=='diagonal_hash':r['diagonal_sha256']='foreign'
            elif case=='source_design_status':r['source_design_disposition']='empty_setting'
            elif case=='promote_science':r['scientific_eligibility']=True
            elif case=='omit_last_audit':records.pop()
            elif case=='duplicate_audit':records[-1]=records[0]
            elif case=='foreign_record_count':r['records']+=1
            else:raise AssertionError(case)
            ap.write_bytes(gzip.compress(''.join(json.dumps(r)+'\n' for r in records).encode(),mtime=0));entry['audit_sha256']=sha(ap)
        mp.write_text(json.dumps(m)+'\n')
        c=json.loads(checkpoint.read_text());c['entry']=entry;checkpoint.write_text(json.dumps(c)+'\n')
        receipt=json.loads(rp.read_text())
        for path in [ap,lp,mp,checkpoint]:receipt['artifacts'][str(path.relative_to(root))]=sha(path)
        rp.write_text(json.dumps(receipt)+'\n')
        rejected(lambda:parallel(plan,True));assert not (root/'readback.json').exists()
        failures=list((root/'failures').glob('*/failure.json'))
        for path in failures:
            assert not json.loads(path.read_text())['scientific_eligibility']
            assert (path.parent/'original_numeric_inputs.npz').is_file();captures+=1
        rejected_cases.append(case)
        print('parallel_private_corruption_rejected',case,flush=True)
    assert len(rejected_cases)==19 and captures>=10
    verify(pins)
    result=dict(status='passed_complete_parallel_four_control_software_qualification_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(),complete_grids=grids,
        serial_parallel_identical_audit_records=compared,serial_parallel_identical_setting_links=linked,
        independent_latent_readbacks=3,qualified_basis_dimensions=sorted(dimensions),
        qualified_nonuniform_policies=sorted(policies),distinct_producer_worker_identities=len(workers),
        rehashed_private_corruptions_rejected=rejected_cases,retained_original_numeric_failure_captures=captures,
        completed_stage_restart_refusals=6,worker_address_space_gib=12,
        cached_numeric_source_arrays_preserved_all_complete_cohorts=True,
        unchanged_comparison_rtol=3e-9,unchanged_comparison_atol=2e-8,
        source_hashes=pins,scientific_eligibility=False,full_real_data_numerical_qualification_complete=False,
        scope='Full three synthetic grids compared to frozen serial component arithmetic and '
              'reconstructed by independent latent arithmetic with two actual fork workers. '
              'All nineteen private corruption cases have distinct producer/output namespaces '
              'and retained original fixtures. No production model fit, biological pilot or repair '
              'of the earlier unexplained native write.')
    with a.receipt.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
