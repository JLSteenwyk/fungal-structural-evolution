#!/usr/bin/env python3
"""Complete order-pair software grid with actual raw MSA/PDB/fit source algorithms."""
import argparse
import contextlib
import csv
from fractions import Fraction
import gzip
import io
import itertools
import json
from pathlib import Path
import numpy as np

import check_full_triad_sequence_geometry as geometry_fixture
import compare_full_triad_correspondence as producer
import readback_full_triad_correspondence_comparison as reader
from duplication_alignment_numeric_readback import load_pdb
from fit_whole_protein_common_residues import fit_triplet
from full_triad_fit_sources import fields as structural_fields,METRICS
from full_triad_correspondence_comparison_sources import blocks,NUMERIC_FIELDS
from full_triad_sequence_fit_sources import iterate_native,triple_sha
from run_ortholog_pair_guide_comparison import sha


def cases(geometry_plan,base):
    gp=json.loads(Path(geometry_plan).read_text())
    ready,inputs,links,sets,baseline,native_db,closure,_=geometry_fixture.reader.load_sources(gp,geometry_plan)
    maps_root=base/'structural_maps';maps_root.mkdir();fits=base/'structural_fits.tsv.gz';mapping=[];rows=[]
    full={}
    for triad,link,native,_ in iterate_native(ready,links,sets,native_db):
        if native['method']=='mafft_auto' and native['permutation']=='012':
            full[triad['triad_id']]=[tuple(t[i] for i in link['role_to_sequence_indices']) for t in native['common_full_triples']]
    coordinates={k:load_pdb(v) for k,v in inputs.items() if k[-1]=='full'}
    for triad in ready:
        models=triad['models'];tid=triad['triad_id'];lengths=[inputs[(*m,'full')]['original_length'] for m in models]
        for mask in ['full','plddt70']:
            allowed=[set(inputs[(*m,mask)]['original_positions']) for m in models]
            all_triples=[t for t in full[tid] if all(t[i] in allowed[i] for i in range(3))]
            for index,order in enumerate(itertools.product([0,1],repeat=3)):
                reference=all_triples[index%3:];cycle=reference[:2] if index==3 else reference[:max(2,len(reference)-index%4)]
                problems=['fixture_source_mask_rejection'] if mask=='plddt70' and inputs[(*models[0],mask)]['status']=='source_rejected' else []
                if models[0][0]=='case1_z_a' and mask=='plddt70' and order[0]:problems.append('fixture_mapping_exclusion')
                inherited=[{s['id']:list(baseline[tid,mask][s['id']]) if edge==2 else [] for s in gp['screens']} for edge in range(3)]
                mapping.append(dict(triad_id=tid,models=models,role_order=['a','b','reference'],mask=mask,orders=list(order),
                    original_lengths=lengths,retained_input_lengths=[inputs[(*m,mask)]['retained_residues'] for m in models],
                    edge_provenance=[dict(mapping_exclusions=problems),dict(mapping_exclusions=[]),dict(mapping_exclusions=[])],
                    edge_pair_screen_exclusions=inherited,edge_pair_screen_pass=[{sid:not why for sid,why in e.items()} for e in inherited],
                    all_three_pair_screen_pass={s['id']:not baseline[tid,mask][s['id']] for s in gp['screens']},
                    reference_common_triples=reference,cycle_consistent_triples=cycle,common_reference_count=len(reference),cycle_consistent_count=len(cycle),
                    common_core_fit_input_status={d:'source_excluded' if problems else 'fewer_than_three_common_residues' if len(t)<3 else 'pending_common_coordinate_geometry' for d,t in [('reference_common',reference),('cycle_consistent',cycle)]}))
                for definition,triples in [('reference_common',reference),('cycle_consistent',cycle)]:
                    n=len(triples);row=dict(triad_id=tid,mask=mask,order_ab=order[0],order_ar=order[1],order_br=order[2],mapping_definition=definition,
                        triples_sha256=triple_sha(triples),common_residues=n,source_exclusions=';'.join(problems),mapping_disagreement_count=len(reference)-len(cycle))
                    for i,role in enumerate(['a','b','reference']):
                        retained=inputs[(*models[i],mask)]['retained_residues']
                        row.update({f'model_{role}':models[i][0],f'version_{role}':models[i][1],f'length_{role}':lengths[i],f'retained_{role}':retained,
                            f'coverage_{role}':n/lengths[i],f'retained_coverage_{role}':n/retained if retained else ''})
                    row.update({k:'' for k in METRICS})
                    if problems:row['fit_status']='source_excluded'
                    elif n<3:row['fit_status']='fewer_than_three_common_residues'
                    else:
                        coords=[];letters=[];confidence=[]
                        for i,model in enumerate(models):
                            seq,xyz,conf=coordinates[(*model,'full')];ix=[t[i]-1 for t in triples]
                            coords.append(xyz[ix]);letters.append(''.join(seq[j] for j in ix));confidence.append(conf[ix])
                        row.update(fit_triplet(coords,letters,confidence))
                    for screen in gp['screens']:
                        sid=screen['id'];why=list(problems);cutoff=Fraction(str(screen['minimum_original_coverage']))
                        if n<screen['minimum_aligned_residues']:why.append('short_common_core')
                        if any(Fraction(n,l)<cutoff for l in lengths):why.append('low_original_protein_coverage')
                        if row['fit_status']=='computed_degenerate_geometry':why.append('degenerate_common_core_geometry')
                        elif row['fit_status']=='fewer_than_three_common_residues':why.append('fewer_than_three_common_residues')
                        inherited_reasons=baseline[tid,mask][sid]
                        row.update({sid+'_core_pass':int(not why),sid+'_core_exclusions':';'.join(why),sid+'_three_pair_pass':int(not inherited_reasons),
                            sid+'_three_pair_exclusions':';'.join(inherited_reasons),sid+'_pass':int(not why and not inherited_reasons),sid+'_exclusions':';'.join(why+inherited_reasons)})
                    rows.append(row)
    with gzip.open(maps_root/'common_residue_maps.jsonl.gz','wt') as f:
        for r in mapping:f.write(json.dumps(r)+'\n')
    with gzip.open(fits,'wt') as f:
        w=csv.DictWriter(f,fieldnames=structural_fields(gp),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    source=dict(ready=ready,inputs=inputs,links=links,sets=sets,native_db=native_db,mapping_root=maps_root,sequence_source=Path(gp['output'])/'sequence_common_residue_fits.tsv.gz',
        structural_source=fits,sequence_plan=gp,structural_plan=gp,ordered_model_triads=7)
    out=base/'comparisons';plan_path=base/'comparison_plan.json'
    plan=dict(output=str(out),screens=gp['screens'],resources=dict(minimum_free_disk_gib=0,maximum_output_gib=1),
        expected=dict(ordered_model_triads=7,source_ready_triads=6,sequence_fit_rows=144,structural_fit_rows=192,sequence_groups=24,comparison_groups=48,pair_states=2304),
        scope='Synthetic full order-grid software check. Closed-source I/O stubbed; real original raw MSA reconstruction, authoritative masks, PDB parsing/SVD fits, source field/grid checks, independent SQL/overlap algorithms and checkpoint/export validation exercised. No production source qualification or biological pilot.')
    plan_path.write_text(json.dumps(plan))
    def sources(config,path):
        paths=[path,native_db,source['sequence_source'],fits,maps_root/'common_residue_maps.jsonl.gz',*[v['path'] for k,v in inputs.items() if k[-1]=='full']]
        return source,{str(p):sha(p) for p in paths}
    producer.load_sources=reader.load_sources=sources
    actual_blocks=producer.blocks
    def interrupt(*args,**kw):
        for i,item in enumerate(actual_blocks(*args,**kw)):
            if i==2:raise RuntimeError('owned_comparison_fixture_interruption')
            yield item
    producer.blocks=interrupt
    try:producer.run(plan_path)
    except RuntimeError as e:assert str(e)=='owned_comparison_fixture_interruption'
    else:raise AssertionError('Missing fixture interruption')
    producer.blocks=actual_blocks
    produced=producer.run(plan_path);checked=reader.run(plan_path,base/'comparison_readback.json')
    assert produced['pair_states']==checked['pair_states']==2304 and produced['checked_reused_triad_blocks']==2
    # All orders of every method/mask/core survive; exceptions and disagreement remain measurable.
    pair_file=out/'correspondence_pair_states.tsv.gz'
    with gzip.open(pair_file,'rt') as f:pairs=list(csv.DictReader(f,delimiter='\t'))
    assert any(r['identical_correspondence']=='0' and int(r['intersection_triples'])>0 for r in pairs)
    assert any(r['strict_numeric_contrast_sign_agreement']=='unavailable' for r in pairs)
    assert any(r['n30_c70_sequence_pass']=='0' for r in pairs)
    receipt=out/'receipt.json';pristine_receipt=receipt.read_bytes();pristine_pair=pair_file.read_bytes();fields=list(pairs[0])
    good=next(i for i,r in enumerate(pairs) if r['both_unique_numeric']=='1' and r['sequence_minus_structural_rmsd_ab']!='')
    changes=[('changed_numeric_difference',good,'sequence_minus_structural_rmsd_ab','9'),('changed_signed_difference',good,'sequence_minus_structural_rmsd_ar_minus_br','1'),
        ('wrong_overlap',good,'intersection_triples','999'),('wrong_jaccard',good,'triple_jaccard','0'),('wrong_sequence_order',good,'sequence_permutation','210'),
        ('wrong_structural_order',good,'structural_order','111'),('wrong_mask',good,'mask','plddt70'),('wrong_triple_hash',good,'sequence_triples_sha256','0'*64),
        ('promoted_screen',good,'n30_c70_both_pass','1'),('removed_pair',None,None,None),('duplicated_pair',None,None,None)]
    rejected=[]
    for name,index,field,value in changes:
        altered=[dict(r) for r in pairs]
        if name=='removed_pair':altered.pop()
        elif name=='duplicated_pair':altered.append(dict(altered[-1]))
        else:
            if altered[index][field]==value:value='0' if value=='1' else 'changed'
            altered[index][field]=value
        with gzip.open(pair_file,'wt') as f:
            w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(altered)
        r=json.loads(pristine_receipt);r['artifacts'][pair_file.name]=sha(pair_file);receipt.write_text(json.dumps(r))
        try:reader.run(plan_path,base/(name+'.json'))
        except (AssertionError,ValueError,StopIteration):rejected.append(name)
        else:raise AssertionError('Accepted false pair export '+name)
    pair_file.write_bytes(pristine_pair);receipt.write_bytes(pristine_receipt)
    group_files=['sequence_order_robustness.jsonl.gz','correspondence_comparison_groups.jsonl.gz'];rejected_groups=[]
    for name in group_files:
        path=out/name;pristine=path.read_bytes()
        with gzip.open(path,'rt') as f:records=[json.loads(s) for s in f]
        for mode in ['range','bitmap','removed_group']:
            altered=json.loads(json.dumps(records))
            if mode=='removed_group':altered.pop()
            elif mode=='range':
                key='metric_ranges' if name.startswith('sequence') else 'metric_difference_ranges'
                field='rmsd_ab' if key=='metric_ranges' else 'sequence_minus_structural_rmsd_ab'
                altered[0][key][field]['mean']+=1
            else:
                screen=altered[0]['screens']['n30_c70'];key='all_orders_pass' if name.startswith('sequence') else 'all_pairs_pass';screen[key]=not screen[key]
            with gzip.open(path,'wt') as f:
                for row in altered:f.write(json.dumps(row)+'\n')
            r=json.loads(pristine_receipt);r['artifacts'][name]=sha(path);receipt.write_text(json.dumps(r))
            try:reader.run(plan_path,base/(name+'_'+mode+'.json'))
            except (AssertionError,ValueError,StopIteration):rejected_groups.append(name+':'+mode)
            else:raise AssertionError('Accepted false group export '+name+':'+mode)
        path.write_bytes(pristine);receipt.write_bytes(pristine_receipt)
    return dict(status='passed_full_triad_correspondence_comparison_software_contracts',source_ready_triads=6,sequence_fit_rows=144,structural_fit_rows=192,
        sequence_groups=24,comparison_groups=48,pair_states=2304,checked_reused_triad_blocks=2,
        maximum_absolute_numeric_difference=checked['maximum_absolute_numeric_difference'],rejected_rehashed_pair_exports=rejected,rejected_rehashed_group_exports=rejected_groups,
        raw_msa_pdb_svd_sources_and_independent_sql_overlap_exercised=True,scope=plan['scope'])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    original=geometry_fixture.reader.run;result=None
    def hook(plan,output):
        nonlocal result
        checked=original(plan,output)
        if Path(output).name=='geometry_readback.json':result=cases(plan,Path(plan).parent)
        return checked
    geometry_fixture.reader.run=hook
    with contextlib.redirect_stdout(io.StringIO()):geometry_fixture.run(args.output.with_suffix('.geometry_fixture.json'))
    assert result is not None
    with args.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
