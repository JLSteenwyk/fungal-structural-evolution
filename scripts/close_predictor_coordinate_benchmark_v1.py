#!/usr/bin/env python3
"""Close complete physical predictor controls, summary and full tree-path handoff."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind,verify

STAGES=[('overlap_coordinate_benchmark',81365,'fungal-overlap-coordinates-20261004-v1.service'),
        ('overlap_coordinate_readback',61890,'fungal-overlap-coordinates-reader-20261004-v1.service'),
        ('overlap_coordinate_summary',75735,'fungal-overlap-coordinates-summary-20261004-v1.service'),
        ('overlap_coordinate_summary_readback',68866,'fungal-overlap-coordinates-summary-reader-20261004-v1.service'),
        ('predictor_coordinate_phylogeny_link',8451,'fungal-coordinate-phylogeny-links-20261004-v1.service'),
        ('predictor_coordinate_phylogeny_link_readback',5100,'fungal-coordinate-phylogeny-links-reader-20261004-v1.service')]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--receipt',type=Path,required=True)
    a=p.parse_args();assert not a.receipt.exists();pins={};stages=[];cpu=0.
    for prefix,session,unit in STAGES:
        path=lambda middle:Path('metadata',prefix+middle+'_20261004_v1.json')
        ep,rp,tp,pp=[path(v) for v in ['_execution','','_transport','_original_tool_payloads']]
        e,r,t,tool=[json.loads(q.read_text()) for q in [ep,rp,tp,pp]]
        for q in [ep,rp,tp,pp]:bind(pins,q)
        assert e['exit_code']==0 and not e['timed_out'] and e['status']=='exited_zero_with_receipt'
        assert e['receipt_sha256']==sha(rp)==t['validation_sha256']
        assert t['original_tool_session_id']==tool['original_tool_session_id']==session
        assert tool['initial']['session_id']==session and tool['terminal']['exit_code']==0
        assert t['unit']==unit and t['invocation_id']==e['invocation_id']
        assert t['whole_wrapper_initial_and_terminal_payloads_matched']
        assert t['manager_start_records']==t['manager_completion_records']==1
        assert t['exact_wrapper_pid_journal_entries']==2 and r['scientific_eligibility'] is False
        for mapping in [e['source_hashes'],e['artifacts'],r['source_hashes'],t['source_hashes']]:
            for q,h in mapping.items():bind(pins,q,h)
        journal=ep.with_suffix('')/'original-invocation-journal.jsonl'
        lines=[json.loads(l) for l in journal.read_text().splitlines()]
        entries=[j for j in lines if j.get('_PID')==str(e['wrapper']['pid']) and j.get('_CMDLINE')==' '.join(e['wrapper']['cmdline'])]
        assert len(entries)==2
        assert json.loads(entries[0]['MESSAGE'])==dict(original_wrapper=e['wrapper'],invocation_id=e['invocation_id'])
        assert json.loads(entries[1]['MESSAGE'])=={k:v for k,v in e.items() if k not in ['source_hashes','artifacts','command','wrapper','child','scope']}
        cpu+=e['child_cpu_seconds']
        stages.append(dict(stage=prefix,original_tool_session_id=session,original_tool_exit_code=0,
            unit=unit,invocation_id=e['invocation_id'],native_wall_seconds=e['wall_seconds'],native_cpu_seconds=e['child_cpu_seconds'],
            child_peak_rss_bytes=e['child_peak_rss_bytes'],receipt=str(rp),receipt_sha256=sha(rp)))
    def read(name):
        q=Path('metadata',name+'_20261004_v1.json');bind(pins,q);return json.loads(q.read_text())
    coordinate=read('overlap_coordinate_readback');summary=read('overlap_coordinate_summary_readback');links=read('predictor_coordinate_phylogeny_link_readback')
    assert coordinate['model_pairs']==643 and coordinate['confidence_rows']==7716 and coordinate['geometry_rows_compared']==5171
    assert coordinate['coordinate_validation_rejected_pairs']==0
    assert summary['threshold_summary_rows']==12 and summary['descriptive_correlation_rows']==96
    assert links['comparison_cases']==8750 and links['original_internal_branch_slots']==4523750
    assert links['mapped_internal_paths']==31290 and links['paired_distribution_occurrences']==187740
    assert links['single_original_branch_paths']==13235 and links['shared_original_branch_paths']==18055
    bind(pins,Path(__file__));verify(pins)
    result=dict(status='complete_verified_full_matched_predictor_coordinate_benchmark_and_phylogenetic_handoff',
        checked_utc=datetime.now(timezone.utc).isoformat(),model_pairs=643,confidence_settings=12,coordinate_rows=7716,
        coordinate_compared_rows=5171,insufficient_coverage_rows=2545,coordinate_rejected_pairs=0,
        source_linked_markers=78,predictor_overlap_taxa=21,branch_ready_markers=71,
        full_marker_slots=125,full_candidate_tree_views=70,full_project_taxon_entries=526,
        full_marker_view_cases=8750,original_internal_branch_slots=4523750,mapped_internal_paths=31290,
        single_original_branch_paths=13235,shared_original_branch_paths=18055,paired_distribution_occurrences=187740,
        original_closed_stages=stages,aggregate_child_cpu_seconds=cpu,all_source_and_output_bindings=len(pins),source_hashes=pins,
        scientific_eligibility=False,full_weighted_biological_fits=0,gpu=False,all_eight_aims_incomplete=True,
        scope='Six original bounded waits closed against exact complete wrapper payloads, manager start/end and all source/output '
              'bindings. Full raw-coordinate/mask/geometry replay, full scalar summary replay and independent full taxon-set '
              'tree-path reconstruction passed. PNG/PDF visually inspected separately. Predictor controls cover selected overlap, '
              'not all fungi/outgroups. Original merged paths/missing data retained. Original formatter/SVD mutation remains '
              'unrepaired; all24native comparison closure, adequate ancestral posterior, accepted phylogenetic framework, full '
              'atlas, timing-gated full models and inferential calibration remain required.')
    with a.receipt.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','original_closed_stages']},indent=2))


if __name__=='__main__':main()
