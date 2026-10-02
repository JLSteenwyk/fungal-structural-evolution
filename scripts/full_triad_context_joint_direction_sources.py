"""Closed joint directions and unchanged original context lineage for all 27 scenarios."""
import itertools
import json
from pathlib import Path
from full_triad_joint_direction_sources import MASKS,METHOD_SCENARIOS,CORE_SCENARIOS,bundle
from full_triad_context_contrast_sources import GATES,POLICIES,PHYSICAL_DIRECTIONS,STATES,SUMMARY_FIELDS
from reference_measurement_union_sources import bind,verify

SCENARIOS=list(itertools.product(MASKS,METHOD_SCENARIOS,CORE_SCENARIOS))
SCENARIO_IDS=['|'.join(s) for s in SCENARIOS]


def load_sources(plan,path):
    bindings=dict(plan['pins']);bind(bindings,path)
    c,root=bundle(plan['joint_completion'],'complete_verified_full_triad_joint_directions',
        'complete_full_triad_joint_directions_pending_independent_readback',
        'passed_full_triad_joint_directions_raw_fit_sql_readback',plan['joint_plan'],bindings)
    assert c['measured_triads']==27056 and c['joint_direction_groups']==730512
    config=json.loads(Path(plan['joint_plan']).read_text())
    assert config['screens']==plan['screens'] and config['contrast_numerical_tolerance_angstrom']==1e-9
    structural=json.loads(Path(config['structural_robustness_plan']).read_text())
    assert structural['triad_work_plan']==plan['triad_work_plan'] and structural['triad_work_completion']==plan['triad_work_completion']
    wc=json.loads(Path(plan['triad_work_completion']).read_text())
    assert wc['status']=='complete_verified_full_reference_triad_work_design' and len(wc['services'])==2
    for p,d in wc['source_hashes'].items():bind(bindings,p,d)
    assert plan['scenario_axes']==[list(s) for s in SCENARIOS]
    expected=dict(target_contexts=283409,reference_tie_records=214461,duplicate_reference_links=428922,
        measured_triads=27056,measured_contrast_groups=730512,logical_reference_screen_decisions=34742682,
        context_screen_states=91824516,context_policy_decisions=459122580,summary_rows=32400)
    assert plan['expected']==expected
    for k in ['target_contexts','reference_tie_records','duplicate_reference_links']:assert wc['summary'][k]==expected[k]
    work=json.loads(Path(plan['triad_work_plan']).read_text());contexts=Path(work['output'])/'context_triad_design.jsonl.gz'
    assert str(contexts) in bindings
    for p in [plan['triad_work_plan'],plan['triad_work_completion']]:bind(bindings,p)
    verify(bindings);return contexts,root/'joint_directions.jsonl.gz',wc['summary'],bindings
