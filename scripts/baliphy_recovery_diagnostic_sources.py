"""Full recovery diagnostic source graph and whole-quartet reuse rules."""
from collections import Counter, defaultdict
import json
from pathlib import Path
from audit_baliphy_memory_recovery_20261002_v2 import build_overlay, CHECKED
from background_measurement_union_sources import closed_source
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def dispatch(overlay):
    rows=overlay['rows'];index={r['chain']['chain_id']:r for r in rows}
    assert len(rows)==len(index)==1620
    groups={q['model_input_identity']:q for q in overlay['quartets']}
    assert len(groups)==len(overlay['quartets'])==405
    result={};seen=set()
    for group,q in sorted(groups.items()):
        members=[index[cid] for cid in q['chain_ids']]
        assert len(members)==4 and not seen.intersection(q['chain_ids'])
        assert all(r['model_input_identity']==group for r in members);seen.update(q['chain_ids'])
        ready=all(r['selected_disposition']['status']==CHECKED for r in members)
        if not ready:
            assert q['status']=='unresolved_failed_chain'
            status='unresolved_failed_native_chain_retained'
        elif all(r['selected_disposition']==r['original_disposition'] for r in members):
            assert q['status']=='complete_integrity_checked_quartet'
            status='reuse_closed_original_unchanged_quartet'
        else:
            assert q['status']=='complete_integrity_checked_quartet'
            assert any(r['selection_origin']=='new_whole_same_seed_integrity_checked_attempt' for r in members)
            status='compute_complete_recovery_quartet'
        result[group]=dict(status=status,chain_ids=q['chain_ids'])
    assert seen==set(index)
    return result


def load(plan,path):
    bindings=dict(plan['pins']);bind(bindings,path)
    recovery=closed_source(plan['recovery_completion'],
        'complete_verified_full_baliphy_memory_recovery_accounting',
        'complete_verified_baliphy_memory_recovery_accounting_archive',2,bindings)
    post=closed_source(plan['postprocessing_completion'],
        'complete_verified_baliphy_postprocessing_accounting',
        'complete_verified_baliphy_postprocessing_accounting_archive',3,bindings)
    source=json.loads(Path(plan['producer_plan']).read_text())
    jobs=json.loads(Path(source['jobs']).read_text())
    original=json.loads((Path(source['output'])/'receipt.json').read_text())
    retry_plan=json.loads(Path(plan['recovery_plan']).read_text())
    retry=json.loads((Path(retry_plan['output'])/'receipt.json').read_text())
    rows,quartets,summary=build_overlay(jobs,original['chains'],retry['recovery_chains'],retry_plan['failed_chain_ids'])
    overlay=json.loads(Path(recovery['overlay']).read_text());bind(bindings,recovery['overlay'],recovery['overlay_sha256'])
    assert overlay['rows']==rows and overlay['quartets']==quartets and overlay['summary']==summary
    assert all(recovery[k]==v for k,v in summary.items()) and source['iterations']==1000
    assert post['state_chains']==1617 and post['categorical_quartets']==post['length_quartets']==402
    original_plans={k:json.loads(Path(p).read_text()) for k,p in plan['original_plans'].items()}
    reports={k:json.loads((Path(p['output'])/'receipt.json').read_text()) for k,p in original_plans.items()}
    routes=dispatch(overlay)
    assert Counter(r['status'] for r in routes.values())=={
        'reuse_closed_original_unchanged_quartet':402,'compute_complete_recovery_quartet':1,
        'unresolved_failed_native_chain_retained':2}
    assert set(reports['states']['chains'])=={r['chain']['chain_id'] for r in rows}
    for name in ['scalar','length','categorical']:
        assert set(reports[name]['groups'])==set(routes)
    verify(bindings)
    return dict(source=source,rows={r['chain']['chain_id']:r for r in rows},
        summary=summary,overlay=overlay,routes=routes,reports=reports,original_plans=original_plans),bindings


def document(path,bindings,digest=None):
    path=Path(path);bind(bindings,path,digest);record=json.loads(path.read_text())
    for field in ['pins','evidence','source_hashes']:
        for p,d in record.get(field,{}).items():bind(bindings,p,d)
    for name,d in record.get('artifacts',{}).items():bind(bindings,path.parent/name,d)
    return record


def summarize(groups,states,source,bindings):
    assert set(groups)==set(source['routes']) and set(states)==set(source['rows'])
    totals={k:Counter() for k in ['scalar_status_counts','length_status_counts',
        'categorical_coordinate_counts','categorical_pattern_counts','coordinate_status_counts',
        'pattern_status_counts','quartets_passing_every_scalar','quartets_passing_every_length_scalar']}
    checked=0;observations=anchors=unanchored=0
    for cid,row in source['rows'].items():
        selected=row['selected_disposition'];info=states[cid]
        if selected['status']!=CHECKED:
            assert info==dict(status='unresolved_failed_native_chain_retained',native_disposition=selected)
            continue
        if selected==row['original_disposition']:
            original=source['reports']['states']['chains'][cid]
            assert info==dict(status='reused_closed_original_state_trace',receipt=original['receipt'],receipt_sha256=original['receipt_sha256'])
        else:assert info['status']=='new_recovery_whole_chain_state_trace'
        receipt=document(info['receipt'],bindings,info['receipt_sha256'])
        assert receipt['status']=='complete_independently_checked_anchored_state_traces' and receipt['chains']==1
        assert len(receipt['summaries'])==1
        summary=receipt['summaries'][0];assert summary['chain']==cid and summary['samples']==101
        assert summary['source_audit_sha256']==selected['sample_audit_sha256']
        assert Path(summary['source_audit']).resolve()==Path(selected['sample_audit']).resolve()
        for p,d in summary['artifacts'].items():bind(bindings,p,d)
        assert summary['state_observations']==summary['candidate_anchor_coordinates']*101
        checked+=1;observations+=summary['state_observations'];anchors+=summary['candidate_anchor_coordinates']
        unanchored+=summary['unanchored_residue_observations']
    complete=0
    for group,route in source['routes'].items():
        info=groups[group];assert info['chain_ids']==route['chain_ids'] and info['origin']==route['status']
        if route['status']=='unresolved_failed_native_chain_retained':
            assert info['status']==route['status'] and set(info)=={'status','origin','chain_ids','scientific_eligibility'}
            assert info['scientific_eligibility'] is False;continue
        assert info['status']=='complete_scalar_length_category_screens_not_posterior_qualification'
        assert info['scientific_eligibility'] is False;complete+=1
        for kind in ['scalar','length']:
            if route['status']=='reuse_closed_original_unchanged_quartet':
                old=source['reports'][kind]['groups'][group]
                assert info[kind]==dict(receipt=old['receipt'],receipt_sha256=old['receipt_sha256'])
            record=document(info[kind]['receipt'],bindings,info[kind]['receipt_sha256'])
            assert set(record['outputs'])=={'250','500'}
            for cutoff,value in record['outputs'].items():
                report=document(value['path'],bindings,value['sha256'])
                assert report['status']=='scalar_diagnostics_complete_not_posterior_qualification'
                assert report['discard_through_iteration']==int(cutoff)
                expected_samples=(750 if cutoff=='250' else 500) if kind=='scalar' else (75 if cutoff=='250' else 50)
                assert all(v['draws_per_chain']==expected_samples for v in report['variables'].values())
                assert report['thresholds']==dict(rhat_strict_upper=1.01,bulk_ess_minimum=400,tail_ess_minimum=400)
                counts=Counter(v['status'] for v in report['variables'].values())
                totals[kind+'_status_counts'].update({cutoff+':'+k:v for k,v in counts.items()})
                totals['quartets_passing_every_'+('scalar' if kind=='scalar' else 'length_scalar')][cutoff]+=all(
                    v['status']=='passes_scalar_screen_only' for v in report['variables'].values())
        category=document(info['categorical']['receipt'],bindings,info['categorical']['receipt_sha256'])
        if route['status']=='reuse_closed_original_unchanged_quartet':
            old=source['reports']['categorical']['groups'][group]
            assert info['categorical']==dict(receipt=old['receipt'],receipt_sha256=old['receipt_sha256'])
        assert category['status']=='verified_quartet_categorical_reports_complete_not_posterior_qualification' and category['group']==group
        manifest=document(category['manifest'],bindings,category['manifest_sha256'])
        report=document(category['report'],bindings,category['report_sha256'])
        assert manifest['status']=='provenance_checked_state_quartet'
        assert {r['chain_id'] for r in manifest['chains']}==set(route['chain_ids']) and len(manifest['chains'])==4
        assert len({r['seed'] for r in manifest['chains']})==4 and all(r['model_input_identity']==group for r in manifest['chains'])
        assert manifest['expected_iterations']==list(range(0,1001,10))
        bind(bindings,manifest['arrays'],manifest['arrays_sha256'])
        coordinates=manifest['coordinates'];assert coordinates['nodes']==sorted(set(coordinates['nodes'])) and len(coordinates['nodes'])==4
        expected=4*sum(r['length'] for r in coordinates['tips'])
        for cid in route['chain_ids']:
            state=document(states[cid]['receipt'],bindings,states[cid]['receipt_sha256'])['summaries'][0]
            cp=[p for p in state['artifacts'] if Path(p).name=='coordinates.json'];assert len(cp)==1
            assert document(cp[0],bindings,state['artifacts'][cp[0]])==coordinates
            assert state['candidate_anchor_coordinates']==expected
        assert set(report['outputs'])=={'250','500'}
        for cutoff,value in report['outputs'].items():
            summary=document(Path(category['report']).parent/('discard-'+cutoff)/'summary.json',bindings)
            assert summary['discard_through']==int(cutoff)
            assert summary['coordinates']==value['coordinates']==expected
            assert summary['patterns']==value['patterns']
            assert summary['retained_samples_per_chain']==value['retained_samples_per_chain']==(75 if cutoff=='250' else 50)
            assert sum(summary['coordinate_status_counts'].values())==expected
            assert sum(summary['pattern_status_counts'].values())==value['patterns']
            totals['categorical_coordinate_counts'][cutoff]+=expected;totals['categorical_pattern_counts'][cutoff]+=value['patterns']
            for key in ['coordinate_status_counts','pattern_status_counts']:
                totals[key].update({cutoff+':'+k:v for k,v in summary[key].items()})
    assert checked==source['summary']['selected_checked_chains'] and complete==source['summary']['complete_quartets']
    return dict(**source['summary'],state_chains=checked,state_observations=observations,
        per_chain_anchor_occurrences=anchors,unanchored_residue_observations=unanchored,
        scalar_quartets=complete,length_quartets=complete,categorical_quartets=complete,
        **{k:dict(v) for k,v in totals.items()})
