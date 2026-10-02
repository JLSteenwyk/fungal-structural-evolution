#!/usr/bin/env python3
"""Independently count all raw-MSA-derived masks and coordinate work."""
import argparse
import base64
from collections import Counter
import json
from pathlib import Path
from readback_full_triad_sequence_alignments import reconstruct
from full_triad_sequence_fit_sources import load_sources,iterate_native,PREFLIGHT_FIELDS
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def inspect(plan_path,output):
    plan=json.loads(Path(plan_path).read_text())
    ready,inputs,links,sets,baseline,database,closure,bindings=load_sources(plan,plan_path)
    rp=Path(plan['output'])/'receipt.json';receipt=json.loads(rp.read_text())
    assert receipt['status']=='complete_full_triad_sequence_fit_preflight_pending_independent_readback'
    assert receipt['plan_sha256']==sha(plan_path)
    counters=Counter();occurrences=Counter();needed=set();states=rows=unique=total=0;last=None;seen=set()
    for triad,link,native,digest in iterate_native(ready,links,sets,database):
        if triad['triad_id']!=last:last=triad['triad_id'];seen=set()
        if native['status']=='valid_sequence_alignment':
            raw=base64.b64decode(native['stdout_base64'],validate=True)
            rebuilt=reconstruct(raw,sets[link['sequence_set_id']]['sequences'])
            assert rebuilt['common_full_triples']==native['common_full_triples']
            source=rebuilt['common_full_triples']
        else:source=[]
        for mask in ['full','plddt70']:
            positions=[inputs[(*m,mask)]['original_positions'] for m in triad['models']]
            tuples=[]
            for original in source:
                actual=tuple(original[j] for j in link['role_to_sequence_indices'])
                if all(p in positions[i] for i,p in enumerate(actual)):tuples.append(actual)
            status='native_sequence_alignment_unusable' if native['status']!='valid_sequence_alignment' else 'pending_common_coordinate_geometry'
            if native['status']=='valid_sequence_alignment' and len(tuples)<3:status='fewer_than_three_common_residues'
            if any(inputs[(*m,mask)]['status']=='source_rejected' for m in triad['models']):status='source_mask_rejected'
            key=mask+':'+native['method'];counters[key+':'+status]+=1;occurrences[key]+=len(tuples);rows+=1
            if status=='pending_common_coordinate_geometry':
                marker=mask,tuple(tuples)
                if marker not in seen:seen.add(marker);unique+=1;total+=len(tuples)
                for model in triad['models']:needed.add((*model,'full'))
        states+=1
        if states%10000==0:print('Independent raw sequence-fit preflight',states,'/',closure['native_alignment_states'],flush=True)
    summary=dict(ordered_model_triads=closure['ordered_model_triads'],source_ready_triads=len(ready),
        native_alignment_states=states,fit_rows=rows,input_status_counts=dict(counters),triple_occurrences=dict(occurrences),
        unique_eligible_cores=unique,unique_eligible_residue_occurrences=total,proper_pair_fits_per_implementation=3*unique,
        needed_full_pdb_inputs=len(needed),needed_full_pdb_bytes=sum(Path(inputs[k]['path']).stat().st_size for k in needed))
    assert rows==2*states==2*closure['native_alignment_states']
    assert all(receipt[k]==summary[k] for k in PREFLIGHT_FIELDS)
    bind(bindings,rp);verify(bindings)
    result=dict(status='passed_full_triad_sequence_fit_raw_preflight_readback',**summary,
        plan_sha256=sha(plan_path),producer_receipt_sha256=sha(rp),source_hashes=bindings,
        scientific_eligibility=False,scope=plan['scope'])
    with Path(output).open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'}),flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();inspect(a.plan,a.output)
