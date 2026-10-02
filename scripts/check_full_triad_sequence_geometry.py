#!/usr/bin/env python3
"""Complete synthetic sequence/order/mask grid, raw readbacks and false exports."""
import argparse
import base64
import contextlib
import csv
import gzip
import hashlib
import io
import itertools
import json
from pathlib import Path
import sqlite3
import tempfile
import zlib
import numpy as np

import fit_full_triad_sequence_correspondence as producer
import readback_full_triad_sequence_geometry as reader
import preflight_full_triad_sequence_fits as preflight
import readback_full_triad_sequence_fit_preflight as preflight_reader
from triad_sequence_native import parse_alignment
from run_ortholog_pair_guide_comparison import sha


def run(output):
    with tempfile.TemporaryDirectory(prefix='sequence-geometry-fixtures-') as temp:
        base=Path(temp);inputs={};ready=[];links={};sets={};baseline={};native_db=base/'native.sqlite'
        screens=[dict(id=f'n{n}_c{int(c*100)}',minimum_aligned_residues=n,minimum_original_coverage=c) for n in [30,50] for c in [.5,.7,.9]]
        db=sqlite3.connect(native_db)
        db.execute('CREATE TABLE alignments(sequence_set_id TEXT,method TEXT,permutation TEXT,payload BLOB,payload_sha256 TEXT,PRIMARY KEY(sequence_set_id,method,permutation))')
        for case in range(5):
            models=[[f'case{case}_{name}',1] for name in ['z_a','a_b','m_reference']]
            tid=hashlib.sha256(json.dumps(models).encode()).hexdigest();ready.append(dict(triad_id=tid,models=models))
            for column,model in enumerate(models):
                length=[60,80,65][column];positions=list(range(1,length+1));t=np.array(positions)
                coords=np.column_stack((t*.8,np.sin(t)*2,np.cos(t*.7)))
                if case==2:coords=np.column_stack((t*.8,np.zeros(length),np.zeros(length)))
                if column==1:
                    coords=coords@np.array([[0.,-1,0],[1.,0,0],[0,0,1.]])+[4,3,2];coords[0,0]+=.4
                if column==2 and case==1:coords[:,0]*=-1
                lines=[];seq=[]
                for serial,(x,y,z) in enumerate(coords,1):
                    aa='GLY' if column==1 and serial%4==0 else 'ALA';seq.append('G' if aa=='GLY' else 'A')
                    confidence=69.996 if serial==2 else 85. if serial%2==0 else 65.
                    lines.append(f'ATOM  {serial:5d}  CA  {aa} A{serial:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{confidence:6.2f}          C  \n')
                path=base/(model[0]+'.pdb');path.write_text(''.join(lines));sequence=''.join(seq)
                inputs[(*model,'full')]=dict(path=str(path),sha256=sha(path),model_id=model[0],version=1,mask='full',
                    sequence=sequence,original_positions=positions,original_length=length,retained_residues=length,status='ready')
                kept=list(range(4,length+1,2))
                inputs[(*model,'plddt70')]=dict(original_positions=kept,original_length=length,retained_residues=len(kept),
                    sequence=''.join(sequence[i-1] for i in kept),status='source_rejected' if case==4 else 'ready')
            sorted_models=sorted(models);sid=hashlib.sha256(json.dumps(sorted_models).encode()).hexdigest()
            sequences=[inputs[(*m,'full')]['sequence'] for m in sorted_models]
            sets[sid]=dict(sequence_set_id=sid,models=sorted_models,sequences=sequences)
            links[tid]=dict(sequence_set_id=sid,role_to_sequence_indices=[sorted_models.index(m) for m in models])
            for mask in ['full','plddt70']:
                baseline[tid,mask]={s['id']:['br:order1_not_numerically_usable'] if case==0 and s['id']=='n30_c70' else [] for s in screens}
            if case==3:
                # Two common columns; every remaining letter occurs in exactly one sequence's block.
                width=sum(len(s)-2 for s in sequences)
                aligned=[];offset=0
                for seq in sequences:
                    tail=seq[2:];aligned.append(seq[:2]+'-'*offset+tail+'-'*(width-offset-len(tail)));offset+=len(tail)
            else:aligned=[seq+'-'*(max(map(len,sequences))-len(seq)) for seq in sequences]
            for method in ['mafft_auto','famsa_default']:
                for order in itertools.permutations(range(3)):
                    permutation=''.join(map(str,order));raw=''.join('>m%d\n%s\n'%(i,aligned[i]) for i in order).encode()
                    invalid=case==1 and method=='famsa_default' and permutation=='210'
                    metrics={k:None for k in ['alignment_columns','common_residues','common_full_triples']} if invalid else parse_alignment(raw,sequences)
                    record=dict(sequence_set_id=sid,models=sorted_models,method=method,permutation=permutation,
                        stdout_base64=base64.b64encode(raw).decode(),status='native_nonzero_exit' if invalid else 'valid_sequence_alignment',**metrics)
                    payload=json.dumps(record,separators=(',',':')).encode()
                    db.execute('INSERT INTO alignments VALUES(?,?,?,?,?)',(sid,method,permutation,zlib.compress(payload),hashlib.sha256(payload).hexdigest()))
        # Reuse an unordered physical source under the opposite duplicate-role projection.
        original=ready[0];models=[original['models'][1],original['models'][0],original['models'][2]]
        tid='swapped_duplicate_roles';ready.append(dict(triad_id=tid,models=models));link=links[original['triad_id']]
        links[tid]=dict(sequence_set_id=link['sequence_set_id'],role_to_sequence_indices=[link['role_to_sequence_indices'][1],link['role_to_sequence_indices'][0],link['role_to_sequence_indices'][2]])
        for mask in ['full','plddt70']:baseline[tid,mask]=baseline[original['triad_id'],mask]
        db.commit();db.close()
        completion=base/'native_completion.json';completion.write_text('{"fixture_only":true}\n')
        pc=base/'preflight_completion.json';pc.write_text('{"fixture_only":true}\n')
        preplan=base/'preflight_plan.json';fitplan=base/'fit_plan.json'
        common=dict(native_completion=str(completion),preflight_completion=str(pc),screens=screens,
            expected=dict(fit_rows=len(ready)*24),resources=dict(minimum_free_disk_gib=0,maximum_output_gib=1),scope='Synthetic software contracts; source closure I/O stubbed, real raw MSA, SQLite, PDB, SVD/quaternion and screens exercised. Not production proof or a taxon pilot.')
        preplan.write_text(json.dumps(dict(common,output=str(base/'preflight'))));fitplan.write_text(json.dumps(dict(common,output=str(base/'geometry'))))
        closure=dict(ordered_model_triads=len(ready)+1,source_ready_triads=len(ready),native_alignment_states=12*len(ready))
        def sources(config,path):
            bindings={str(path):sha(path),str(native_db):sha(native_db),str(completion):sha(completion)}
            return ready,inputs,links,sets,baseline,native_db,closure,bindings
        for module in [producer,reader,preflight,preflight_reader]:module.load_sources=sources
        with contextlib.redirect_stdout(io.StringIO()):
            pr=preflight.run(preplan);pa=preflight_reader.inspect(preplan,base/'preflight_readback.json')
        assert pr['fit_rows']==pa['fit_rows']==144
        def gate(config,bindings):bindings[str(pc)]=sha(pc);return pr
        producer.require_preflight=reader.require_preflight=gate
        # Interrupt only this owned synthetic iterator after one committed triad.
        actual_iterator=producer.iterate_native
        def interrupted(*args):
            for i,item in enumerate(actual_iterator(*args)):
                if i==16:raise RuntimeError('owned_fixture_interruption')
                yield item
        producer.iterate_native=interrupted
        try:
            with contextlib.redirect_stdout(io.StringIO()):producer.run(fitplan)
        except RuntimeError as error:assert str(error)=='owned_fixture_interruption'
        else:raise AssertionError('Fixture interruption did not occur')
        producer.iterate_native=actual_iterator
        with contextlib.redirect_stdout(io.StringIO()):
            made=producer.run(fitplan);checked=reader.run(fitplan,base/'geometry_readback.json')
        assert made['fit_rows']==checked['fit_rows']==144
        assert made['checked_reused_fit_rows']==24
        table=base/'geometry'/'sequence_common_residue_fits.tsv.gz';rp=table.parent/'receipt.json'
        with gzip.open(table,'rt') as f:records=list(csv.DictReader(f,delimiter='\t'))
        find=lambda **kw:next(i for i,r in enumerate(records) if all(r[k]==v for k,v in kw.items()))
        normal=find(model_a='case0_z_a',mask='full',method='mafft_auto',permutation='012')
        mirrored=find(model_a='case1_z_a',mask='full',method='mafft_auto',permutation='012')
        degenerate=find(fit_status='computed_degenerate_geometry');short=find(fit_status='fewer_than_three_common_residues')
        excluded=find(fit_status='source_mask_rejected');unusable=find(fit_status='native_sequence_alignment_unusable')
        masked=find(model_a='case0_z_a',mask='plddt70',method='mafft_auto',permutation='012')
        swapped=find(model_a='case0_a_b',mask='full',method='mafft_auto',permutation='012')
        assert records[normal]['coverage_b']=='0.75' and records[normal]['n30_c70_core_pass']=='1' and records[normal]['n30_c70_pass']=='0'
        assert records[masked]['common_residues']=='29' and records[masked]['removed_for_joint_mask']=='31'
        assert float(records[normal]['joint_authoritative_plddt70_fraction'])==29/60
        assert float(records[normal]['joint_plddt70_fraction'])==30/60
        assert float(records[mirrored]['rmsd_ar'])>.1
        assert abs(float(records[normal]['rmsd_ar_minus_br']))>.001
        assert np.isclose(float(records[normal]['rmsd_ar_minus_br']),-float(records[swapped]['rmsd_ar_minus_br']))
        assert all(records[i]['rmsd_ab']=='' for i in [short,excluded,unusable])
        pristine_table=table.read_bytes();pristine_receipt=rp.read_bytes();columns=list(records[0])
        changes=[('changed_rmsd',normal,'rmsd_ab','9'),('changed_signed_contrast',normal,'rmsd_ar_minus_br','1'),
            ('wrong_original_denominator',normal,'coverage_b','1'),('wrong_mask_positions',masked,'common_residues','30'),
            ('changed_residue_identity',normal,'sequence_identity_ab','0'),('changed_confidence',normal,'mean_plddt_a','0'),
            ('changed_triple_hash',normal,'triples_sha256','0'*64),('swapped_model_role',normal,'model_a',records[normal]['model_b']),
            ('promoted_degenerate',degenerate,'fit_status','computed_unique_at_numeric_tolerance'),
            ('promoted_short',short,'fit_status','computed_unique_at_numeric_tolerance'),('cleared_source_rejection',excluded,'mask_input_status_a','ready'),
            ('cleared_native_failure',unusable,'sequence_native_status','valid_sequence_alignment'),
            ('promoted_inherited_screen',normal,'n30_c70_three_pair_pass','1'),('changed_native_input_order',normal,'permutation','021'),
            ('removed_state',None,None,None),('duplicated_state',None,None,None)]
        rejected=[]
        for name,index,field,value in changes:
            altered=[dict(r) for r in records]
            if name=='removed_state':altered.pop()
            elif name=='duplicated_state':altered.append(dict(altered[-1]))
            else:altered[index][field]=value
            with gzip.open(table,'wt') as f:
                writer=csv.DictWriter(f,fieldnames=columns,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(altered)
            r=json.loads(pristine_receipt);r['artifacts']['sequence_common_residue_fits.tsv.gz']=sha(table);rp.write_text(json.dumps(r))
            try:
                with contextlib.redirect_stdout(io.StringIO()):reader.run(fitplan,base/(name+'.json'))
            except (AssertionError,ValueError):rejected.append(name)
            else:raise AssertionError('Accepted false geometry export '+name)
        table.write_bytes(pristine_table);rp.write_bytes(pristine_receipt)
        pre_receipt=base/'preflight'/'receipt.json';original_pre=pre_receipt.read_bytes();rejected_pre=[]
        for field in ['fit_rows','unique_eligible_cores','needed_full_pdb_bytes']:
            r=json.loads(original_pre);r[field]+=1;pre_receipt.write_text(json.dumps(r))
            try:
                with contextlib.redirect_stdout(io.StringIO()):preflight_reader.inspect(preplan,base/('false_pre_'+field+'.json'))
            except (AssertionError,ValueError):rejected_pre.append(field)
            else:raise AssertionError('Accepted false preflight '+field)
        result=dict(status='passed_full_triad_sequence_geometry_software_contracts',source_model_sets=len(sets),ordered_role_triads=len(ready),native_alignment_states=72,fit_rows=144,
            raw_msa_and_pdb_algorithms_exercised=True,original_length_coverage_and_nonconsecutive_masks=True,
            rounded70_excluded_by_authoritative_source_mask=True,duplicate_role_sign_reversal=True,proper_rotation_reflection_case=True,
            interrupted_checkpoint_restart_passed=True,maximum_absolute_rmsd_or_contrast_difference=checked['maximum_absolute_rmsd_or_contrast_difference'],
            checked_reused_fit_rows=made['checked_reused_fit_rows'],
            rejected_rehashed_false_geometry_exports=rejected,rejected_false_preflight_summaries=rejected_pre,scope=common['scope'])
        with Path(output).open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result,indent=2));return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);run(p.parse_args().output)
