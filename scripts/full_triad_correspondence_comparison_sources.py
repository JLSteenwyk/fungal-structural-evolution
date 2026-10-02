"""Complete closed geometry grids and residue sets for correspondence comparison."""
import base64
import csv
import gzip
import itertools
import json
from pathlib import Path

from full_triad_sequence_sources import load_original
from full_triad_sequence_fit_sources import iterate_native,triple_sha,fields as sequence_fields
from full_triad_fit_sources import iterate_maps,fields as structural_fields
from full_triad_robustness_sources import NUMERIC_FIELDS,DEFINITIONS
from readback_full_triad_sequence_alignments import reconstruct
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha

METHODS=['famsa_default','mafft_auto']
PERMUTATIONS=[''.join(map(str,p)) for p in itertools.permutations(range(3))]
SEQUENCE_NUMERIC_FIELDS=NUMERIC_FIELDS+['joint_authoritative_plddt70_fraction']
SUMMARY_FIELDS=['ordered_model_triads','source_ready_triads','sequence_fit_rows','structural_fit_rows',
    'sequence_groups','comparison_groups','pair_states','pair_status_counts','sequence_direction_counts',
    'joint_direction_counts','both_screen_pass_counts','all_pair_screen_pass_group_counts','maximum_absolute_metric_difference']
PAIR_FIELDS=['triad_id','mask','sequence_method','sequence_permutation','mapping_definition','structural_order',
    'sequence_fit_status','structural_fit_status','sequence_triples_sha256','structural_triples_sha256',
    'sequence_common_residues','structural_common_residues','intersection_triples','union_triples','triple_jaccard',
    'identical_correspondence','both_unique_numeric','strict_numeric_contrast_sign_agreement']


def pair_fields(plan):
    return PAIR_FIELDS+['sequence_minus_structural_'+f for f in NUMERIC_FIELDS]+[
        s['id']+suffix for s in plan['screens'] for suffix in ['_sequence_pass','_structural_pass','_both_pass']]


def closed(path,status,source_plan,bindings,archive_status=None):
    c=json.loads(Path(path).read_text());assert c['status']==status and c['exact_process_journals_checked']==2
    bind(bindings,path);bind(bindings,c['full_hash_archive'],c['full_hash_archive_sha256'])
    archive=json.loads(Path(c['full_hash_archive']).read_text())
    assert archive['status']==(archive_status or status+'_archive') and len(archive['services'])==2
    assert len(archive['source_hashes'])==c['bound_source_hashes']
    for q,d in archive['source_hashes'].items():bind(bindings,q,d)
    rp=Path(c['producer_receipt']);ap=Path(c['independent_readback'])
    bind(bindings,rp,c['producer_receipt_sha256']);bind(bindings,ap,c['independent_readback_sha256'])
    r,a=[json.loads(p.read_text()) for p in [rp,ap]]
    expected_statuses={
        'complete_verified_full_triad_sequence_correspondence_geometry':(
            'complete_full_triad_sequence_geometry_pending_independent_readback','passed_full_triad_sequence_raw_msa_quaternion_geometry_readback'),
        'complete_verified_full_triad_same_residue_geometry':(
            'complete_full_triad_same_residue_fits_pending_independent_readback','passed_full_triad_same_residue_quaternion_readback')}
    assert (r['status'],a['status'])==expected_statuses[status]
    assert r['plan_sha256']==a['plan_sha256']==sha(source_plan)
    assert a['producer_receipt_sha256']==sha(rp) and r['scientific_eligibility'] is a['scientific_eligibility'] is False
    for k,v in archive['summary'].items():
        if k in r and k in a:assert r[k]==a[k]==c[k]==v
    for name,d in r['artifacts'].items():bind(bindings,rp.parent/name,d)
    return c,r,rp.parent


def load_sources(plan,plan_path):
    cp=Path(plan['catalog_source_plan']);catalog_plan=json.loads(cp.read_text())
    catalog,ready,inputs,mapping,bindings=load_original(catalog_plan,cp)
    bind(bindings,plan_path)
    for q,d in plan['pins'].items():bind(bindings,q,d)
    sc,sr,sequence_root=closed(plan['sequence_geometry_completion'],'complete_verified_full_triad_sequence_correspondence_geometry',plan['sequence_geometry_plan'],bindings)
    tc,tr,structural_root=closed(plan['structural_geometry_completion'],'complete_verified_full_triad_same_residue_geometry',plan['structural_geometry_plan'],bindings,
        'complete_verified_full_triad_same_residue_fit_archive')
    assert sc['source_ready_triads']==sr['source_ready_triads']==len(ready)==27056
    assert sc['ordered_model_triads']==len(catalog)==31235 and sc['fit_rows']==sr['fit_rows']==649344
    assert tc['correspondence_work_triads']==tr['correspondence_work_triads']==len(ready) and tc['fit_rows']==tr['fit_rows']==865792
    assert plan['expected']==dict(ordered_model_triads=31235,source_ready_triads=27056,sequence_fit_rows=649344,structural_fit_rows=865792,
        sequence_groups=108224,comparison_groups=216448,pair_states=10389504)
    sequence_plan=json.loads(Path(plan['sequence_geometry_plan']).read_text());structural_plan=json.loads(Path(plan['structural_geometry_plan']).read_text())
    assert plan['structural_geometry_plan']==catalog_plan['geometry_source_plan']
    assert sequence_plan['catalog_source_plan']==plan['catalog_source_plan']
    assert sequence_plan['screens']==structural_plan['screens']==plan['screens']
    assert sr['native_completion_sha256']==sha(sequence_plan['native_completion'])
    native_c=json.loads(Path(sequence_plan['native_completion']).read_text())
    assert native_c['status']=='complete_verified_full_triad_sequence_alignment_dispositions' and native_c['exact_process_journals_checked']==2
    catalog_root=Path(catalog_plan['output'])
    links={r['triad_id']:r for r in map(json.loads,(catalog_root/'triad_sequence_links.jsonl').open())}
    sets={r['sequence_set_id']:r for r in map(json.loads,(catalog_root/'sequence_sets.jsonl').open())}
    assert len(links)==len(catalog) and len(sets)==27056
    for t in ready:
        link=links[t['triad_id']];record=sets[link['sequence_set_id']]
        assert [record['models'][i] for i in link['role_to_sequence_indices']]==t['models']
        assert record['sequences']==[inputs[(*m,'full')]['sequence'] for m in record['models']]
    mapping_root=Path(json.loads(Path(structural_plan['mapping_plan']).read_text())['output'])
    native_db=Path(native_c['producer_receipt']).parent/'native_alignments.sqlite'
    assert all(str(p) in bindings for p in [native_db,mapping_root/'common_residue_maps.jsonl.gz',
        catalog_root/'triad_sequence_links.jsonl',catalog_root/'sequence_sets.jsonl',
        sequence_root/'sequence_common_residue_fits.tsv.gz',structural_root/'common_residue_fits.tsv.gz'])
    verify(bindings)
    return dict(ready=ready,inputs=inputs,links=links,sets=sets,native_db=native_db,mapping_root=mapping_root,
        sequence_source=sequence_root/'sequence_common_residue_fits.tsv.gz',structural_source=structural_root/'common_residue_fits.tsv.gz',
        sequence_plan=sequence_plan,structural_plan=structural_plan,ordered_model_triads=len(catalog)),bindings


def blocks(source,raw_msa=False):
    ready=source['ready'];inputs=source['inputs']
    native=iter(iterate_native(ready,source['links'],source['sets'],source['native_db']))
    maps=iter(iterate_maps(source['mapping_root'],ready,inputs))
    with gzip.open(source['sequence_source'],'rt') as sf,gzip.open(source['structural_source'],'rt') as tf:
        seq=csv.DictReader(sf,delimiter='\t');struct=csv.DictReader(tf,delimiter='\t')
        assert seq.fieldnames==sequence_fields(source['sequence_plan']) and struct.fieldnames==structural_fields(source['structural_plan'])
        for triad in ready:
            sequence={};structural={};sequence_triples={};structural_triples={}
            for method in METHODS:
                for permutation in PERMUTATIONS:
                    t,link,n,digest=next(native);assert t==triad and (n['method'],n['permutation'])==(method,permutation)
                    if n['status']=='valid_sequence_alignment':
                        full=n['common_full_triples']
                        if raw_msa:
                            rebuilt=reconstruct(base64.b64decode(n['stdout_base64'],validate=True),source['sets'][link['sequence_set_id']]['sequences'])
                            assert rebuilt['common_full_triples']==full;full=rebuilt['common_full_triples']
                    else:assert n['common_full_triples'] is None;full=[]
                    for mask in ['full','plddt70']:
                        row=next(seq,None);assert row is not None
                        assert (row['triad_id'],row['method'],row['permutation'],row['mask'])==(triad['triad_id'],method,permutation,mask)
                        assert row['native_payload_sha256']==digest and row['sequence_set_id']==link['sequence_set_id'] and row['sequence_native_status']==n['status']
                        positions=[set(inputs[(*m,mask)]['original_positions']) for m in triad['models']]
                        projected=[tuple(t[i] for i in link['role_to_sequence_indices']) for t in full]
                        triples=tuple(t for t in projected if all(t[i] in positions[i] for i in range(3)))
                        assert row['triples_sha256']==triple_sha(triples) and int(row['common_residues'])==len(triples)
                        key=mask,method,permutation;sequence[key]=row;sequence_triples[key]=triples
            for mask in ['full','plddt70']:
                for index,order in enumerate(itertools.product([0,1],repeat=3)):
                    mapping=next(maps);assert (mapping['triad_id'],mapping['mask'],mapping['orders'])==(triad['triad_id'],mask,list(order))
                    for definition,field in [('reference_common','reference_common_triples'),('cycle_consistent','cycle_consistent_triples')]:
                        row=next(struct,None);assert row is not None
                        assert (row['triad_id'],row['mask'],row['mapping_definition'])==(triad['triad_id'],mask,definition)
                        assert tuple(int(row['order_'+e]) for e in ['ab','ar','br'])==order
                        triples=tuple(tuple(t) for t in mapping[field]);assert row['triples_sha256']==triple_sha(triples) and int(row['common_residues'])==len(triples)
                        key=mask,definition,index;structural[key]=row;structural_triples[key]=triples
            for row in [*sequence.values(),*structural.values()]:
                assert [[row['model_'+r],int(row['version_'+r])] for r in ['a','b','reference']]==triad['models']
                assert [int(row['length_'+r]) for r in ['a','b','reference']]==[inputs[(*m,'full')]['original_length'] for m in triad['models']]
            yield triad,sequence,structural,sequence_triples,structural_triples
        assert next(seq,None) is next(struct,None) is None
    assert next(native,None) is next(maps,None) is None
