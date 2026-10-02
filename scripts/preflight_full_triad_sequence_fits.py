#!/usr/bin/env python3
"""Exhaustively inventory sequence-derived fit work before coordinate fitting."""
import argparse
from collections import Counter
import json
from pathlib import Path
from full_triad_sequence_fit_sources import load_sources,iterate_native,triple_sha
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def run(plan_path):
    plan=json.loads(Path(plan_path).read_text())
    ready,inputs,links,sets,baseline,database,closure,bindings=load_sources(plan,plan_path)
    counts=Counter();occurrences=Counter();needed=set();states=rows=unique=residues=0;previous=None;seen=set()
    for triad,link,native,digest in iterate_native(ready,links,sets,database):
        if triad['triad_id']!=previous:previous=triad['triad_id'];seen=set()
        full=[tuple(t[i] for i in link['role_to_sequence_indices']) for t in (native['common_full_triples'] or [])]
        for mask in ['full','plddt70']:
            allowed=[set(inputs[(*m,mask)]['original_positions']) for m in triad['models']]
            triples=[t for t in full if all(t[i] in allowed[i] for i in range(3))]
            assert all(len({t[i] for t in triples})==len(triples) for i in range(3))
            status=('native_sequence_alignment_unusable' if native['status']!='valid_sequence_alignment'
                    else 'fewer_than_three_common_residues' if len(triples)<3 else 'pending_common_coordinate_geometry')
            if any(inputs[(*m,mask)]['status']=='source_rejected' for m in triad['models']):status='source_mask_rejected'
            key=mask+':'+native['method'];counts[key+':'+status]+=1;occurrences[key]+=len(triples);rows+=1
            if status!='pending_common_coordinate_geometry':continue
            marker=mask,triple_sha(triples)
            if marker not in seen:seen.add(marker);unique+=1;residues+=len(triples)
            needed.update((*m,'full') for m in triad['models'])
        states+=1
        if states%10000==0:print('Full sequence-fit work preflight',states,'/',closure['native_alignment_states'],flush=True)
    assert states==closure['native_alignment_states'] and rows==2*states
    assert all(inputs[k]['status']=='ready' for k in needed)
    summary=dict(ordered_model_triads=closure['ordered_model_triads'],source_ready_triads=len(ready),
        native_alignment_states=states,fit_rows=rows,input_status_counts=dict(counts),triple_occurrences=dict(occurrences),
        unique_eligible_cores=unique,unique_eligible_residue_occurrences=residues,proper_pair_fits_per_implementation=3*unique,
        needed_full_pdb_inputs=len(needed),needed_full_pdb_bytes=sum(Path(inputs[k]['path']).stat().st_size for k in needed))
    verify(bindings);out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    result=dict(status='complete_full_triad_sequence_fit_preflight_pending_independent_readback',**summary,
        plan_sha256=sha(plan_path),source_hashes=bindings,artifacts={},scientific_eligibility=False,scope=plan['scope'])
    with (out/'receipt.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'}),flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);run(p.parse_args().plan)
