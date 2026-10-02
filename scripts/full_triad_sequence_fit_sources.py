"""Closed sequence/coordinate sources and original three-pair coverage gates."""
import hashlib
import itertools
import json
from pathlib import Path
import sqlite3
import zlib

from full_triad_sequence_sources import load_original
from full_triad_fit_sources import iterate_maps, METRICS
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha

SUMMARY_FIELDS=['ordered_model_triads','source_ready_triads','native_alignment_states','fit_rows',
                'counts','core_screen_pass_counts','screen_pass_counts','triple_occurrences','fitted_pdb_inputs']
PREFLIGHT_FIELDS=['ordered_model_triads','source_ready_triads','native_alignment_states','fit_rows',
                 'input_status_counts','triple_occurrences','unique_eligible_cores',
                 'unique_eligible_residue_occurrences','proper_pair_fits_per_implementation',
                 'needed_full_pdb_inputs','needed_full_pdb_bytes']
BASE_FIELDS=['triad_id','sequence_set_id','method','permutation','native_payload_sha256',
             'sequence_native_status','mask','triples_sha256','source_full_common_residues',
             'removed_for_joint_mask','common_residues','joint_authoritative_plddt70_fraction']


def fields(plan):
    result=BASE_FIELDS.copy()
    for role in ['a','b','reference']:
        result += ['model_'+role,'version_'+role,'length_'+role,'retained_'+role,
                   'mask_input_status_'+role,'coverage_'+role,'retained_coverage_'+role]
    result+=METRICS
    for screen in plan['screens']:
        result += [screen['id']+suffix for suffix in ['_core_pass','_core_exclusions',
            '_three_pair_pass','_three_pair_exclusions','_pass','_exclusions']]
    return result


def triple_sha(triples):
    return hashlib.sha256(json.dumps(triples,separators=(',',':')).encode()).hexdigest()


def load_sources(plan,plan_path):
    cp=Path(plan['catalog_source_plan']);catalog_plan=json.loads(cp.read_text())
    catalog,ready,inputs,mapping,bindings=load_original(catalog_plan,cp)
    bind(bindings,plan_path)
    for path,digest in plan['pins'].items():bind(bindings,path,digest)
    closure=json.loads(Path(plan['native_completion']).read_text())
    assert closure['status']=='complete_verified_full_triad_sequence_alignment_dispositions'
    assert closure['ordered_model_triads']==31235 and closure['source_ready_triads']==27056
    assert closure['native_alignment_states']==324672 and closure['exact_process_journals_checked']==2
    assert plan['expected']==dict(ordered_model_triads=31235,source_ready_triads=27056,
                                 native_alignment_states=324672,fit_rows=649344)
    bind(bindings,plan['native_completion'])
    bind(bindings,closure['full_hash_archive'],closure['full_hash_archive_sha256'])
    archive=json.loads(Path(closure['full_hash_archive']).read_text());assert len(archive['services'])==2
    for path,digest in archive['source_hashes'].items():bind(bindings,path,digest)
    np=Path(plan['native_source_plan']);native_plan=json.loads(np.read_text())
    rp=Path(closure['producer_receipt']);receipt=json.loads(rp.read_text())
    assert receipt['plan_sha256']==sha(np)
    assert sha(rp)==closure['producer_receipt_sha256']
    assert closure['producer_receipt']==str(Path(native_plan['output'])/'receipt.json')
    # Closure includes the catalog's complete source proof and link artifacts.
    catalog_root=Path(catalog_plan['output'])
    links=[json.loads(line) for line in (catalog_root/'triad_sequence_links.jsonl').open()]
    assert len(links)==len(catalog)==31235
    link_index={row['triad_id']:row for row in links};assert len(link_index)==31235
    sets={row['sequence_set_id']:row for row in map(json.loads,(catalog_root/'sequence_sets.jsonl').open())}
    assert len(sets)==closure['unique_sequence_model_sets']==27056
    for triad in ready:
        link=link_index[triad['triad_id']];record=sets[link['sequence_set_id']]
        assert link['sequence_control_disposition']=='scheduled_full_sequence_alignment'
        assert [record['models'][i] for i in link['role_to_sequence_indices']]==triad['models']
        for column,model in enumerate(record['models']):
            assert record['sequences'][column]==inputs[(*model,'full')]['sequence']
    original=json.loads(Path(catalog_plan['geometry_source_plan']).read_text())
    assert original['screens']==plan['screens']
    map_root=Path(json.loads(Path(original['mapping_plan']).read_text())['output'])
    baseline={}
    for row in iterate_maps(map_root,ready,inputs):
        key=row['triad_id'],row['mask']
        reasons={screen['id']:[label+':'+reason for label,edge in zip(['ab','ar','br'],row['edge_pair_screen_exclusions'])
            for reason in edge[screen['id']]] for screen in plan['screens']}
        assert all(row['all_three_pair_screen_pass'][sid]==(not why) for sid,why in reasons.items())
        if key in baseline:assert baseline[key]==reasons,'Both-order pair gates vary across structural order views'
        else:baseline[key]=reasons
    assert len(baseline)==2*len(ready)
    database=rp.parent/'native_alignments.sqlite'
    assert archive['source_hashes'][str(database)]==receipt['artifacts']['native_alignments.sqlite']
    verify(bindings)
    return ready,inputs,link_index,sets,baseline,database,closure,bindings


def iterate_native(ready,links,sets,database):
    db=sqlite3.connect('file:'+str(database)+'?mode=ro',uri=True)
    try:
        assert db.execute('SELECT COUNT(*) FROM alignments').fetchone()[0]==12*len(sets)
        for triad in ready:
            link=links[triad['triad_id']];sid=link['sequence_set_id'];seen=set()
            for method,permutation,payload,digest in db.execute(
                'SELECT method,permutation,payload,payload_sha256 FROM alignments WHERE sequence_set_id=? ORDER BY method,permutation',(sid,)):
                raw=zlib.decompress(payload);assert hashlib.sha256(raw).hexdigest()==digest;record=json.loads(raw)
                assert record['sequence_set_id']==sid and record['method']==method and record['permutation']==permutation
                assert record['models']==sets[sid]['models']
                assert record['status'] in ['valid_sequence_alignment','native_timeout','native_nonzero_exit','invalid_native_alignment']
                seen.add((method,permutation))
                yield triad,link,record,digest
            assert seen=={(method,''.join(map(str,order))) for method in ['mafft_auto','famsa_default'] for order in itertools.permutations(range(3))}
    finally:db.close()


def require_preflight(plan,bindings):
    c=json.loads(Path(plan['preflight_completion']).read_text())
    assert c['status']=='complete_verified_full_triad_sequence_fit_preflight' and c['exact_process_journals_checked']==2
    assert c['source_ready_triads']==27056 and c['native_alignment_states']==324672 and c['fit_rows']==649344
    bind(bindings,plan['preflight_completion']);bind(bindings,c['full_hash_archive'],c['full_hash_archive_sha256'])
    archive=json.loads(Path(c['full_hash_archive']).read_text());assert len(archive['services'])==2
    assert archive['source_hashes'][plan['native_completion']]==sha(plan['native_completion'])
    for path,digest in archive['source_hashes'].items():bind(bindings,path,digest)
    verify(bindings)
    return c
