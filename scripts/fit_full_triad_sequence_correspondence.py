#!/usr/bin/env python3
"""Fit the complete sequence-derived residue grid with checked checkpoints."""
import argparse
from collections import Counter
import csv
import fcntl
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import sqlite3
import zlib

from duplication_alignment_numeric_readback import load_pdb
from fit_whole_protein_common_residues import fit_triplet
from full_triad_sequence_fit_sources import load_sources,iterate_native,require_preflight,fields,triple_sha,METRICS
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def run(plan_path):
    plan=json.loads(Path(plan_path).read_text());ph=sha(plan_path)
    assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
    ready,inputs,links,sets,baseline,database,closure,bindings=load_sources(plan,plan_path)
    preflight=require_preflight(plan,bindings)
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=True)
    assert not (out/'receipt.json').exists(),'Completed fit output must remain immutable'
    lock=(out/'run.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    configuration=out/'configuration.json'
    config=dict(plan_sha256=ph,source_hashes=bindings)
    if configuration.exists():assert json.loads(configuration.read_text())==config,'Changed checkpoint inputs'
    else:
        with configuration.open('x') as f:f.write(json.dumps(config,indent=2)+'\n')
    dbpath=out/'fit_checkpoints.sqlite';db=sqlite3.connect(dbpath)
    db.execute('CREATE TABLE IF NOT EXISTS fits (triad_id TEXT, method TEXT, permutation TEXT, mask TEXT, payload BLOB, digest TEXT, PRIMARY KEY(triad_id,method,permutation,mask))');db.commit()
    fitted=set();counts=Counter();passes=Counter();core_passes=Counter();occurrences=Counter();rows=states=reused=0
    @lru_cache(maxsize=1024)
    def pdb(key):
        source=inputs[key];assert source['status']=='ready'
        bind(bindings,source['path'],source['sha256']);fitted.add(key)
        letters,coords,confidence=load_pdb(source)
        return letters,coords,confidence,{p:i for i,p in enumerate(source['original_positions'])}
    @lru_cache(maxsize=64)
    def calculate(models,triples):
        coords=[];letters=[];confidence=[]
        for column,model in enumerate(models):
            seq,xyz,conf,index=pdb((*model,'full'));indices=[index[t[column]] for t in triples]
            coords.append(xyz[indices]);letters.append(''.join(seq[i] for i in indices));confidence.append(conf[indices])
        return fit_triplet(coords,letters,confidence)
    try:
        with gzip.open(out/'sequence_common_residue_fits.tsv.gz','wt',compresslevel=1) as f:
            writer=csv.DictWriter(f,fieldnames=fields(plan),delimiter='\t',lineterminator='\n');writer.writeheader()
            for triad,link,native,digest in iterate_native(ready,links,sets,database):
                models=tuple(tuple(m) for m in triad['models'])
                full=tuple(tuple(t[i] for i in link['role_to_sequence_indices']) for t in (native['common_full_triples'] or []))
                if native['status']!='valid_sequence_alignment':assert not full
                lengths=[inputs[(*m,'full')]['original_length'] for m in models]
                for mask in ['full','plddt70']:
                    allowed=[set(inputs[(*m,mask)]['original_positions']) for m in models]
                    triples=tuple(t for t in full if all(t[i] in allowed[i] for i in range(3)));n=len(triples)
                    assert all(len({t[i] for t in triples})==n for i in range(3))
                    retained=[inputs[(*m,mask)]['retained_residues'] for m in models]
                    high=[set(inputs[(*m,'plddt70')]['original_positions']) for m in models]
                    row=dict(triad_id=triad['triad_id'],sequence_set_id=link['sequence_set_id'],method=native['method'],permutation=native['permutation'],
                        native_payload_sha256=digest,sequence_native_status=native['status'],mask=mask,triples_sha256=triple_sha(triples),
                        source_full_common_residues=len(full),removed_for_joint_mask=len(full)-n,common_residues=n,
                        joint_authoritative_plddt70_fraction=sum(all(t[i] in high[i] for i in range(3)) for t in triples)/n if n else '')
                    for i,role in enumerate(['a','b','reference']):
                        row.update({f'model_{role}':models[i][0],f'version_{role}':models[i][1],f'length_{role}':lengths[i],
                            f'retained_{role}':retained[i],f'mask_input_status_{role}':inputs[(*models[i],mask)]['status'],
                            f'coverage_{role}':n/lengths[i],f'retained_coverage_{role}':n/retained[i] if retained[i] else ''})
                    row.update({k:'' for k in METRICS})
                    if any(inputs[(*m,mask)]['status']=='source_rejected' for m in models):row['fit_status']='source_mask_rejected'
                    elif native['status']!='valid_sequence_alignment':row['fit_status']='native_sequence_alignment_unusable'
                    elif n<3:row['fit_status']='fewer_than_three_common_residues'
                    else:row.update(calculate(models,triples))
                    for screen in plan['screens']:
                        sid=screen['id'];why=[];cutoff=Fraction(str(screen['minimum_original_coverage']))
                        if n<screen['minimum_aligned_residues']:why.append('short_common_core')
                        if any(Fraction(n,length)<cutoff for length in lengths):why.append('low_original_protein_coverage')
                        if row['fit_status']=='computed_degenerate_geometry':why.append('degenerate_common_core_geometry')
                        elif row['fit_status']!='computed_unique_at_numeric_tolerance':why.append(row['fit_status'])
                        inherited=baseline[triad['triad_id'],mask][sid]
                        row.update({sid+'_core_pass':int(not why),sid+'_core_exclusions':';'.join(why),
                            sid+'_three_pair_pass':int(not inherited),sid+'_three_pair_exclusions':';'.join(inherited),
                            sid+'_pass':int(not why and not inherited),sid+'_exclusions':';'.join(why+inherited)})
                        key=mask+':'+native['method']+':'+sid;passes[key]+=row[sid+'_pass'];core_passes[key]+=row[sid+'_core_pass']
                    key=(triad['triad_id'],native['method'],native['permutation'],mask)
                    raw=json.dumps(row,separators=(',',':')).encode();rh=hashlib.sha256(raw).hexdigest()
                    previous=db.execute('SELECT payload,digest FROM fits WHERE triad_id=? AND method=? AND permutation=? AND mask=?',key).fetchone()
                    if previous is None:db.execute('INSERT INTO fits VALUES(?,?,?,?,?,?)',(*key,zlib.compress(raw),rh))
                    else:
                        stored=zlib.decompress(previous[0]);assert hashlib.sha256(stored).hexdigest()==previous[1]
                        assert stored==raw,'Changed or false fit checkpoint';reused+=1
                    writer.writerow(row);rows+=1;counts[mask+':'+native['method']+':'+row['fit_status']]+=1;occurrences[mask+':'+native['method']]+=n
                states+=1
                if states%12==0:db.commit()
                if states%128==0:
                    db.commit()
                    state=dict(status='running_full_triad_sequence_geometry',completed=rows,total=plan['expected']['fit_rows'],counts=dict(counts),checked_reused_fit_rows=reused,plan_sha256=ph)
                    temp=out/'state.tmp';temp.write_text(json.dumps(state,indent=2)+'\n');temp.replace(out/'state.json')
                    assert dbpath.stat().st_size+f.tell()<=plan['resources']['maximum_output_gib']*2**30,'Output budget exceeded'
                    assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
                if states%10000==0:print('Full sequence-derived geometry',states,'/',closure['native_alignment_states'],flush=True)
        db.commit();assert db.execute('SELECT COUNT(*) FROM fits').fetchone()[0]==rows
        assert states==closure['native_alignment_states'] and rows==2*states==preflight['fit_rows']
        assert dict(occurrences)==preflight['triple_occurrences'] and len(fitted)==preflight['needed_full_pdb_inputs']
    finally:db.close();lock.close()
    bind(bindings,configuration);verify(bindings)
    summary=dict(ordered_model_triads=closure['ordered_model_triads'],source_ready_triads=len(ready),native_alignment_states=states,fit_rows=rows,
        counts=dict(counts),core_screen_pass_counts=dict(core_passes),screen_pass_counts=dict(passes),triple_occurrences=dict(occurrences),fitted_pdb_inputs=len(fitted))
    result=dict(status='complete_full_triad_sequence_geometry_pending_independent_readback',**summary,plan_sha256=ph,
        native_completion_sha256=sha(plan['native_completion']),preflight_completion_sha256=sha(plan['preflight_completion']),
        checked_reused_fit_rows=reused,source_hashes=bindings,artifacts={name:sha(out/name) for name in ['sequence_common_residue_fits.tsv.gz','fit_checkpoints.sqlite']},
        scientific_eligibility=False,scope=plan['scope'])
    with (out/'receipt.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'}),flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);run(p.parse_args().plan)
