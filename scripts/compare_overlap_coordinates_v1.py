#!/usr/bin/env python3
"""Benchmark all643 matched model pairs against coordinates at all12 confidence settings."""
import argparse
from collections import Counter
import csv
from datetime import datetime,timezone
import json
from pathlib import Path

import numpy as np

from ancestral_chain_attempt import sha
from overlap_coordinate_benchmark_v1 import METRICS,load_model,measure,identity,pair_id
from compare_predictor_alphabets import compare_states
from reference_measurement_union_sources import bind,verify


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True)
    p.add_argument('--qualification',type=Path,required=True)
    a=p.parse_args();assert not a.output.exists() and not a.receipt.exists()
    context=Path('results/phylogeny/paired-source-model-context-20260926-v1')
    features=Path('results/phylogeny/overlap-model-feature-comparison-20260926-v1')
    pins={}
    qualification=json.loads(a.qualification.read_text())
    assert qualification['status']=='passed_overlap_coordinate_software_controls_v1'
    assert qualification['scientific_eligibility'] is False
    for q,d in qualification['source_hashes'].items():bind(pins,q,d)
    bind(pins,a.qualification)
    for root,proof,status in [(context,'metadata/paired_source_model_context_completed_20260926.json',
        'passed_paired_source_model_context_readback'),(features,'metadata/overlap_model_feature_comparison_completed_20260926.json',
        'passed_full_overlap_model_feature_comparison_readback')]:
        rp=root/'receipt.json';r=json.loads(rp.read_text());ap=Path(proof);audit=json.loads(ap.read_text())
        assert audit['status']==status and audit['producer_receipt_sha256']==sha(rp)
        assert audit['producer_receipt']==str(rp)
        for name,d in r['artifacts'].items():bind(pins,root/name,d)
        for q,d in r.get('source_sha256',{}).items():bind(pins,q,d)
        for q in [rp,ap]:bind(pins,q)
    verify(pins)
    pairs={};links=list(csv.DictReader((context/'source_model_context.tsv').open(),delimiter='\t'))
    assert len(links)==673
    for row in links:
        assert row['full_sequence_status']=='identical_complete_encoded_sequence'
        key=identity(row)
        if key in pairs:
            for field in row:
                if field.endswith(('_model_path','_model_sha256','_encoding_path','_encoding_sha256','_sequence_sha256','_length')):
                    assert row[field]==pairs[key][field]
        pairs[key]=row
    assert len(pairs)==643
    original={}
    for row in csv.DictReader((features/'model_pair_comparisons.tsv').open(),delimiter='\t'):
        key=tuple(row[k] for k in ['reference_model_id','reference_version','local_model_id','local_version'])
        conf=(int(row['plddt_cutoff']),row['pae_cutoff'])
        assert (key,conf) not in original;original[key,conf]=row
    assert len(original)==7716
    a.output.mkdir();table=a.output/'all_model_pair_coordinates.tsv';rows=Counter();n=0
    with table.open('x') as f:
        writer=None
        for index,(key,link) in enumerate(sorted(pairs.items())):
            error='';loaded=[]
            try:
                loaded=[load_model(link,prefix,pins) for prefix in ['reference','local']]
            except (ValueError,KeyError,AssertionError) as problem:
                error=type(problem).__name__+': '+str(problem)
            for cutoff in [0,70,90]:
                for pae in [None,5,10,15]:
                    old=original[key,(cutoff,str(pae) if pae is not None else 'unfiltered')]
                    metrics={name:'' for name in METRICS}
                    if error:status='coordinate_validation_rejected'
                    else:
                        (ea,x),(eb,y)=loaded;state,mask=compare_states(ea,eb,cutoff,pae)
                        assert all(old[k]==str(v) for k,v in state.items()),'Confidence/state cohort differs'
                        if state['status']!='compared':status='insufficient_common_coverage'
                        else:
                            indices=np.flatnonzero(mask)
                            metrics=measure(x[indices],y[indices],indices+1);status='coordinates_compared'
                    out=dict(model_pair_id=pair_id(key),reference_model_id=key[0],reference_version=key[1],
                        local_model_id=key[2],local_version=key[3],sequence_sha256=link['reference_sequence_sha256'],
                        plddt_cutoff=cutoff,pae_cutoff=pae if pae is not None else 'unfiltered',
                        protein_length=int(old['protein_length']),retained_residues=int(old['retained_residues']),
                        feature_coverage_status=old['status'],state_mismatches=int(old['state_mismatches']),
                        partner_changes=int(old['partner_changes']),same_partner_residues=int(old['same_partner_residues']),
                        same_partner_state_mismatches=int(old['same_partner_state_mismatches']),
                        changed_partner_residues=int(old['changed_partner_residues']),
                        changed_partner_state_mismatches=int(old['changed_partner_state_mismatches']),
                        status=status,coordinate_validation_issue=error,**metrics,scientific_eligibility=False)
                    if writer is None:writer=csv.DictWriter(f,fieldnames=list(out),delimiter='\t');writer.writeheader()
                    writer.writerow(out);rows[status]+=1;n+=1
            print('overlap_coordinate_pairs',index+1,'/643',flush=True)
    assert n==7716 and sum(rows.values())==7716
    for q in [Path(__file__),Path('scripts/overlap_coordinate_benchmark_v1.py'),
              Path('scripts/compare_marker_structures.py'),Path('scripts/compare_predictor_alphabets.py'),
              Path('scripts/extract_domain_coordinates.py')]:bind(pins,q)
    verify(pins)
    result=dict(status='complete_full_matched_predictor_coordinate_benchmark_pending_independent_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(),model_pairs=643,marker_taxon_links=673,
        confidence_settings=12,comparison_rows=n,status_counts=dict(rows),table=str(table),table_sha256=sha(table),
        source_hashes=pins,scientific_eligibility=False,new_native_fits=0,gpu=False,all_eight_aims_incomplete=True,
        scope='Every original exact-complete-sequence matched model pair and all original12feature-confidence '
              'settings retained. Raw coordinate sequences, all atoms/CA and confidence values checked; '
              'proper-rotation superposition plus all-pair, spatial-contact union<=15A/separation>=3, '
              'and original-sequence-gap3to10distance differences measured on exactly the original '
              'state-confidence masks. Whole proteins, not only fitted marker columns or conserved '
              'domains; orientation and local conformation may both contribute. Source-file and numeric '
              'memory checks preserve inputs. No additivity, evolutionary branch assignment, experimental '
              'accuracy, causation or independent model-pair sampling inferred. Full independent '
              'coordinate/mask/numeric replay and linkage to original phylogenetic paths remain required.')
    with a.receipt.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
