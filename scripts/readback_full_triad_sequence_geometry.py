#!/usr/bin/env python3
"""Rebuild all residue matches from raw MSA and fit independent quaternion rotations."""
import argparse
import base64
from collections import Counter
import csv
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
import json
from pathlib import Path
import sqlite3
import zlib
import numpy as np

from duplication_alignment_numeric_readback import load_pdb
from readback_full_triad_sequence_alignments import reconstruct
from readback_whole_protein_common_fits import core_metrics,compare_row
from full_triad_sequence_fit_sources import load_sources,iterate_native,require_preflight,fields,triple_sha,METRICS,SUMMARY_FIELDS
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def run(plan_path,output):
    plan=json.loads(Path(plan_path).read_text());ph=sha(plan_path)
    ready,inputs,links,sets,baseline,database,closure,bindings=load_sources(plan,plan_path)
    preflight=require_preflight(plan,bindings)
    out=Path(plan['output']);rp=out/'receipt.json';receipt=json.loads(rp.read_text());rh=sha(rp)
    assert receipt['status']=='complete_full_triad_sequence_geometry_pending_independent_readback' and receipt['plan_sha256']==ph
    assert receipt['native_completion_sha256']==sha(plan['native_completion'])
    assert receipt['preflight_completion_sha256']==sha(plan['preflight_completion'])
    assert receipt['artifacts']=={name:sha(out/name) for name in ['sequence_common_residue_fits.tsv.gz','fit_checkpoints.sqlite']}
    config=json.loads((out/'configuration.json').read_text());assert config==dict(plan_sha256=ph,source_hashes=bindings)
    fitted=set();counts=Counter();passes=Counter();core_passes=Counter();occurrences=Counter();rows=states=0;maximum=0.
    @lru_cache(maxsize=1024)
    def residues(key):
        source=inputs[key];assert source['status']=='ready';bind(bindings,source['path'],source['sha256']);fitted.add(key)
        seq,coords,conf=load_pdb(source)
        return {p:(coords[i],seq[i],conf[i]) for i,p in enumerate(source['original_positions'])}
    @lru_cache(maxsize=64)
    def recompute(models,triples):
        sources=[residues((*m,'full')) for m in models]
        parts=[[sources[column][t[column]] for t in triples] for column in range(3)]
        return core_metrics([np.array([v[0] for v in part]) for part in parts],
            [''.join(v[1] for v in part) for part in parts],[[v[2] for v in part] for part in parts])
    db=sqlite3.connect('file:'+str(out/'fit_checkpoints.sqlite')+'?mode=ro',uri=True)
    try:
        with gzip.open(out/'sequence_common_residue_fits.tsv.gz','rt') as f:
            reader=csv.DictReader(f,delimiter='\t');assert reader.fieldnames==fields(plan)
            for triad,link,native,digest in iterate_native(ready,links,sets,database):
                models=tuple(tuple(m) for m in triad['models'])
                if native['status']=='valid_sequence_alignment':
                    raw=base64.b64decode(native['stdout_base64'],validate=True)
                    rebuilt=reconstruct(raw,sets[link['sequence_set_id']]['sequences'])
                    assert all(rebuilt[k]==native[k] for k in ['alignment_columns','common_residues','common_full_triples'])
                    full=[tuple(original[column] for column in link['role_to_sequence_indices']) for original in rebuilt['common_full_triples']]
                else:
                    assert native['common_full_triples'] is None;full=[]
                lengths=[inputs[(*model,'full')]['original_length'] for model in models]
                for mask in ['full','plddt70']:
                    positions=[inputs[(*model,mask)]['original_positions'] for model in models]
                    triples=[]
                    for match in full:
                        if all(p in positions[i] for i,p in enumerate(match)):triples.append(match)
                    triples=tuple(triples);n=len(triples);retained=[len(p) for p in positions]
                    high=[inputs[(*m,'plddt70')]['original_positions'] for m in models]
                    wanted=dict(triad_id=triad['triad_id'],sequence_set_id=link['sequence_set_id'],method=native['method'],permutation=native['permutation'],
                        native_payload_sha256=digest,sequence_native_status=native['status'],mask=mask,triples_sha256=triple_sha(triples),
                        source_full_common_residues=len(full),removed_for_joint_mask=len(full)-n,common_residues=n,
                        joint_authoritative_plddt70_fraction=sum(all(p in high[i] for i,p in enumerate(t)) for t in triples)/n if n else '')
                    for i,role in enumerate(['a','b','reference']):
                        wanted.update({f'model_{role}':models[i][0],f'version_{role}':models[i][1],f'length_{role}':lengths[i],
                            f'retained_{role}':retained[i],f'mask_input_status_{role}':inputs[(*models[i],mask)]['status'],
                            f'coverage_{role}':n/lengths[i],f'retained_coverage_{role}':n/retained[i] if retained[i] else ''})
                    wanted.update({k:'' for k in METRICS})
                    if any(inputs[(*m,mask)]['status']=='source_rejected' for m in models):wanted['fit_status']='source_mask_rejected'
                    elif native['status']!='valid_sequence_alignment':wanted['fit_status']='native_sequence_alignment_unusable'
                    elif n<3:wanted['fit_status']='fewer_than_three_common_residues'
                    else:wanted.update(recompute(models,triples))
                    for screen in plan['screens']:
                        sid=screen['id'];why=[];cutoff=Fraction(str(screen['minimum_original_coverage']))
                        if n<screen['minimum_aligned_residues']:why.append('short_common_core')
                        if any(n*cutoff.denominator<length*cutoff.numerator for length in lengths):why.append('low_original_protein_coverage')
                        status=wanted['fit_status']
                        if status=='computed_degenerate_geometry':why.append('degenerate_common_core_geometry')
                        elif status!='computed_unique_at_numeric_tolerance':why.append(status)
                        inherited=baseline[triad['triad_id'],mask][sid]
                        wanted.update({sid+'_core_pass':int(not why),sid+'_core_exclusions':';'.join(why),
                            sid+'_three_pair_pass':int(not inherited),sid+'_three_pair_exclusions':';'.join(inherited),
                            sid+'_pass':int(not why and not inherited),sid+'_exclusions':';'.join(why+inherited)})
                        key=mask+':'+native['method']+':'+sid;core_passes[key]+=not why;passes[key]+=not why and not inherited
                    actual=next(reader,None);assert actual is not None,'Missing full sequence fit row'
                    maximum=max(maximum,compare_row(actual,wanted))
                    saved=db.execute('SELECT payload,digest FROM fits WHERE triad_id=? AND method=? AND permutation=? AND mask=?',
                        (triad['triad_id'],native['method'],native['permutation'],mask)).fetchone();assert saved is not None
                    raw=zlib.decompress(saved[0]);assert hashlib.sha256(raw).hexdigest()==saved[1]
                    checkpoint={k:str(v) for k,v in json.loads(raw).items()};assert actual==checkpoint,'Checkpoint/export disagree'
                    rows+=1;counts[mask+':'+native['method']+':'+wanted['fit_status']]+=1;occurrences[mask+':'+native['method']]+=n
                states+=1
                if states%10000==0:print('Independent raw-MSA/quaternion sequence geometry',states,'/',closure['native_alignment_states'],flush=True)
            assert next(reader,None) is None,'Extra full sequence fit row'
        assert db.execute('SELECT COUNT(*) FROM fits').fetchone()[0]==rows
    finally:db.close()
    assert states==closure['native_alignment_states'] and rows==2*states==preflight['fit_rows']
    assert dict(occurrences)==preflight['triple_occurrences'] and len(fitted)==preflight['needed_full_pdb_inputs']
    bind(bindings,out/'configuration.json')
    summary=dict(ordered_model_triads=closure['ordered_model_triads'],source_ready_triads=len(ready),native_alignment_states=states,fit_rows=rows,
        counts=dict(counts),core_screen_pass_counts=dict(core_passes),screen_pass_counts=dict(passes),triple_occurrences=dict(occurrences),fitted_pdb_inputs=len(fitted))
    assert all(receipt[k]==summary[k] for k in SUMMARY_FIELDS) and receipt['source_hashes']==bindings
    verify(bindings);assert sha(rp)==rh
    for name,digest in receipt['artifacts'].items():assert sha(out/name)==digest
    bind(bindings,rp)
    result=dict(status='passed_full_triad_sequence_raw_msa_quaternion_geometry_readback',**summary,plan_sha256=ph,producer_receipt_sha256=rh,
        maximum_absolute_rmsd_or_contrast_difference=maximum,source_hashes=bindings,scientific_eligibility=False,scope=plan['scope'])
    with Path(output).open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'}),flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.plan,a.output)
