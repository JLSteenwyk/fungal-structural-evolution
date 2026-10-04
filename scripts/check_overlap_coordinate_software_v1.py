#!/usr/bin/env python3
"""Qualify coordinate arithmetic and full source universe before the complete benchmark."""
import argparse
import csv
from datetime import datetime,timezone
import json
from pathlib import Path

import numpy as np
from scipy.spatial.distance import pdist
from scipy.spatial.transform import Rotation

from overlap_coordinate_benchmark_v1 import measure,load_model,identity
from compare_predictor_alphabets import compare_states
from reference_measurement_union_sources import bind,verify


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--receipt',type=Path,required=True)
    a=p.parse_args();assert not a.receipt.exists();pins={}
    context=Path('results/phylogeny/paired-source-model-context-20260926-v1/source_model_context.tsv')
    features=Path('results/phylogeny/overlap-model-feature-comparison-20260926-v1/model_pair_comparisons.tsv')
    for q in [context,features,Path(__file__),Path('scripts/overlap_coordinate_benchmark_v1.py'),
              Path('scripts/compare_overlap_coordinates_v1.py'),Path('scripts/compare_marker_structures.py'),
              Path('scripts/compare_predictor_alphabets.py'),Path('scripts/extract_domain_coordinates.py')]:bind(pins,q)
    links=list(csv.DictReader(context.open(),delimiter='\t'));pairs={identity(r):r for r in links}
    assert len(links)==673 and len(pairs)==643
    assert all(r['full_sequence_status']=='identical_complete_encoded_sequence' for r in links)
    old={}
    for r in csv.DictReader(features.open(),delimiter='\t'):
        key=tuple(r[k] for k in ['reference_model_id','reference_version','local_model_id','local_version'])
        conf=(int(r['plddt_cutoff']),r['pae_cutoff']);assert (key,conf) not in old;old[key,conf]=r
    assert set(old)=={(k,(c,p)) for k in pairs for c in [0,70,90] for p in ['unfiltered','5','10','15']}
    x=np.array([[0.,0,0],[1,0,0],[0,2,0],[0,0,3],[2,3,4]])
    pos=np.array([1,2,3,4,9]);r=np.array([[0.,-1,0],[1,0,0],[0,0,1]])
    controls=[]
    for label,y in [('identity',x.copy()),('proper_rotation_translation',x@r+np.array([8.,10,-7])),
                    ('reflection',x*np.array([-1.,1,1])),('scale',2*x)]:
        result=measure(x,y,pos)
        assert result['all_distance_pairs']==10
        rot,_=Rotation.align_vectors(x-x.mean(0),y-y.mean(0))
        independent=float(np.sqrt(np.mean(np.sum((rot.apply(y-y.mean(0))-(x-x.mean(0)))**2,axis=1))))
        assert abs(independent-result['ca_superposition_rmsd_angstrom'])<1e-13
        assert abs(float(np.sqrt(np.mean((pdist(x)-pdist(y))**2)))-result['all_distance_rms_change_angstrom'])<1e-13
        if label in ['identity','proper_rotation_translation']:
            assert max(v for k,v in result.items() if k.endswith('_angstrom') and v!='')<1e-13
        elif label=='reflection':
            assert result['ca_superposition_rmsd_angstrom']>.1 and result['all_distance_rms_change_angstrom']==0
        else:assert result['all_distance_rms_change_angstrom']>0
        controls.append(dict(control=label,metrics=result))
    gap=measure(x,2*x,np.array([1,100,200,300,400]))
    assert gap['sequence_local_distance_pairs']==0 and gap['sequence_local_distance_rms_change_angstrom']==''
    refused=0
    for xx,pp in [(x[:2],pos[:2]),(np.full_like(x,np.nan),pos),(x,np.array([1,1,2,3,4])),
                  (x,pos.astype(float)),(x,np.array([0,1,2,3,4]))]:
        try:measure(xx,xx,pp)
        except (AssertionError,ValueError):refused+=1
        else:raise AssertionError('Invalid coordinate/position control accepted')
    assert refused==5
    key=sorted(pairs)[0];loaded=[load_model(pairs[key],prefix,pins) for prefix in ['reference','local']]
    parser_rows=0
    for cutoff in [0,70,90]:
        for pae in [None,5,10,15]:
            state,mask=compare_states(loaded[0][0],loaded[1][0],cutoff,pae)
            assert all(old[key,(cutoff,str(pae) if pae is not None else 'unfiltered')][k]==str(v) for k,v in state.items())
            if state['status']=='compared':measure(loaded[0][1][mask],loaded[1][1][mask],np.flatnonzero(mask)+1)
            parser_rows+=1
    verify(pins)
    result=dict(status='passed_overlap_coordinate_software_controls_v1',checked_utc=datetime.now(timezone.utc).isoformat(),
        controls=controls,original_position_gap_control=gap,invalid_controls_refused=refused,
        full_declared_pairs=643,full_declared_rows=7716,actual_parser_control_pairs=1,actual_parser_control_masks=parser_rows,
        prior_software_failure=dict(stage='initial_unsaved_pure_geometry_control',failure='METRICS schema AssertionError',
            correction='Two sequence_local metric labels now match the produced distance fields',native_full_stage_launched=False),
        arithmetic_control_absolute_tolerance=1e-13,source_hashes=pins,scientific_eligibility=False,
        scope='Software controls and static complete-cohort contracts; one original parser fixture exercises twelve masks. '
              'Machine-rounding tolerance only applies to these arithmetic controls. Full benchmark and full independent '
              'raw-coordinate replay still required. No scientific inference or broader sampler/SVD repair.')
    with a.receipt.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','controls']},indent=2))


if __name__=='__main__':main()
