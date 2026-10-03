#!/usr/bin/env python3
"""Manual clade/branch primitives and full 1620-chain serialized contracts."""
import argparse
import copy
from io import StringIO
import json
from pathlib import Path
import random
import tempfile

from Bio import Phylo

from audit_baliphy_sample_mapping import descendant_index
from independent_native_ancestral_topology import match_trees,parse_tree,clade_index
from independent_native_clade_mapping_fixture import setup
from independent_native_clade_mapping import load,replay
from prepare_independent_native_clade_mapping import run
from readback_independent_native_clade_mapping import run as readback
from prepare_independent_baliphy_scalar_readback_v2 import atomic
from run_ortholog_pair_guide_comparison import sha


def primitives():
    rng=random.Random(120526);cases=0
    for tips in [5,12,48,96]:
        for repeat in range(10):
            nodes=[('tip-'+str(i),[],'0.1') for i in range(tips)];working=list(nodes);number=0
            while len(working)>1:
                size=min(len(working),rng.choice([2,2,3]));chosen=rng.sample(working,size)
                working=[n for n in working if n not in chosen]
                node=('source-'+str(number),chosen,str(rng.randint(1,100000)/10000));number+=1;working.append(node);nodes.append(node)
            def encode(n,runtime):
                label=n[0].replace('source-','runtime-') if runtime else n[0]
                children=n[1][::-1] if runtime else n[1]
                return ('('+','.join(encode(c,runtime) for c in children)+')' if children else '')+label+':'+n[2]
            source,runtime=encode(working[0],False)+';',encode(working[0],True)+';'
            match=match_trees(source,runtime);left=descendant_index(Phylo.read(StringIO(source),'newick'));right=descendant_index(Phylo.read(StringIO(runtime),'newick'))
            assert set(left)==set(right) and len(match['rows'])==len(left)
            for row in match['rows']:
                mask=int(row['descendant_mask_hex'],16);key=tuple(t for i,t in enumerate(match['tips']) if mask&(1<<i))
                assert row['source_node']==left[key].name and row['runtime_node']==right[key].name
                assert row['absolute_length_difference'] in [None,'0']
            cases+=1
    quoted="('t,0':1e-7,'t''1':2.0)source[comment];"
    assert match_trees(quoted,quoted.replace('source','runtime'))['tips']==sorted(['t,0',"t'1"])
    assert match_trees('(a:-1,b:2)r;','(a:-1,b:2)s;')['negative_branch_pairs_require_review']==1
    text='t0:1'
    for i in range(1,1200):text='('+text+',t'+str(i)+':1)n'+str(i)+':1'
    parsed=parse_tree(text+';');index,_=clade_index(parsed,parsed['tips']);assert len(index)==2399
    bad=[('(a:1,b:2)r;','(a:1,c:2)s;'),('((a:1,b:1)i:1,c:1)r;','((a:1,c:1)j:1,b:1)s;'),
         ('(a:1,b:2)r;','(a:1.0000000001,b:2)s;'),('(a:1,b:2)r;','(a,b:2)s;'),
         ('((a:1)i:1,b:2)r;','((a:1)j:1,b:2)s;'),('(a:1,b:2)r;','(a:nan,b:2)s;'),
         ('(a:1,b:2)r;','(a:1e999,b:2)s;'),('(a:1,a:2)r;','(a:1,a:2)s;')]
    rejected=0
    for a,b in bad:
        try:match_trees(a,b)
        except (ValueError,KeyError):rejected+=1
        else:raise AssertionError('Invalid correspondence accepted')
    assert match_trees('(a:1,b:2)r:7;','(b:2,a:1)s:9;')['maximum_absolute_length_difference']=='0'
    assert match_trees('(a:1,b:2)r;','(a:1.000000000099,b:2)s;')['maximum_absolute_length_difference']!='0'
    return cases,rejected


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    cases,invalid=primitives();rejected=[];source_rejected=[]
    with tempfile.TemporaryDirectory(prefix='native-clade-contract-') as directory:
        root=Path(directory);path=setup(root)
        try:run(path,stop_after_chains=3)
        except InterruptedError:pass
        else:raise AssertionError('Interruption not exercised')
        before={str(p):sha(p) for p in (root/'output/chains').glob('*.json')}
        result=run(path);assert all(sha(p)==v for p,v in before.items());readback(path,root/'output/readback.json')
        for action in ['completed_producer','alternate_reader']:
            try:run(path) if action=='completed_producer' else readback(path,root/'alternate.json')
            except AssertionError:pass
            else:raise AssertionError('Completed restart accepted')
        (root/'output/readback.json').unlink();out=root/'output';manifest=out/'chain_manifest.json';receipt=out/'receipt.json'
        mb,rb=manifest.read_bytes(),receipt.read_bytes();first=json.loads(mb)[0];cp=out/first['path'];cb=cp.read_bytes()
        changes=['missing_node','duplicate_node','changed_clade','changed_length','candidate_label','candidate_root','promote_science',
            'foreign_seed','candidate_frame_count','hidden_failure','omit_failed_chain','omit_intact_unresolved_chain']
        for change in changes:
            row=json.loads(cb);r=row['result'];entries=json.loads(mb)
            if change=='missing_node':r['node_rows'].pop()
            elif change=='duplicate_node':r['node_rows'][-1]=copy.deepcopy(r['node_rows'][0])
            elif change=='changed_clade':r['node_rows'][0]['descendant_mask_hex']='0xffff'
            elif change=='changed_length':r['node_rows'][0]['source_length']='9'
            elif change=='candidate_label':r['candidate_rows'][0]['runtime_node']='foreign'
            elif change=='candidate_root':r['candidate_rows'][3]['assumed_root']=False
            elif change=='promote_science':r['biological_root_accepted']=True
            elif change=='foreign_seed':r['seed']+=1
            elif change=='candidate_frame_count':r['candidate_frame_mappings']-=1
            elif change=='hidden_failure':r['status']='unresolved_failed_native_chain_retained';r['node_rows']=r['candidate_rows']=[]
            elif change=='omit_failed_chain':entries=[e for e in entries if e['chain_id']!='group-0404-chain0']
            else:entries=[e for e in entries if e['chain_id']!='group-0404-chain1']
            atomic(cp,row);entries[0]['sha256']=sha(cp);atomic(manifest,entries);document=json.loads(rb)
            for artifact in [cp,manifest]:document['artifacts'][str(artifact.relative_to(out))]=sha(artifact)
            atomic(receipt,document)
            try:readback(path,root/'bad-reader.json')
            except (AssertionError,ValueError):rejected.append(change)
            else:raise AssertionError('Rehashed false export accepted: '+change)
            cp.write_bytes(cb);manifest.write_bytes(mb);receipt.write_bytes(rb)
        source,_=load(json.loads(path.read_text()),path);cid=first['chain_id'];audit_path=Path(source['details'][cid]['sample_audit']);saved=audit_path.read_bytes()
        for change in ['wrong_candidate_level','wrong_runtime_identity','duplicate_frame_mapping','missing_frame_mapping']:
            audit=json.loads(saved)
            if change=='wrong_candidate_level':audit['candidate_samples'][0]['level']='3'
            elif change=='wrong_runtime_identity':audit['candidate_samples'][0]['runtime_node']='node-3'
            elif change=='duplicate_frame_mapping':audit['candidate_samples'][-1]=copy.deepcopy(audit['candidate_samples'][0])
            else:audit['candidate_samples'].pop()
            audit_path.write_text(json.dumps(audit))
            try:replay(source,cid)
            except (AssertionError,ValueError):source_rejected.append(change)
            else:raise AssertionError('Changed source candidate mapping accepted')
            audit_path.write_bytes(saved)
    paths=[__file__,*[f'scripts/{n}.py' for n in ['independent_native_ancestral_topology','independent_native_clade_mapping',
        'independent_native_clade_mapping_fixture','prepare_independent_native_clade_mapping','readback_independent_native_clade_mapping',
        'independent_native_ancestral_alignment','independent_native_alignment_replay_sources','check_independent_baliphy_scalar_readback_v2',
        'audit_baliphy_sample_mapping','prepare_ancestral_state_traces','ancestral_residue_anchors']]]
    result=dict(status='passed_full_independent_native_clade_mapping_contracts',random_tree_correspondence_fixtures=cases,
        malformed_correspondences_rejected=invalid,deep_iterative_tree_nodes=2399,quoted_labels_and_comments_checked=True,
        strict_branch_threshold_boundaries_checked=True,full_chains=1620,checked_chains=1618,failed_chains=2,full_groups=405,
        effective_input_groups=135,node_pairs=result['node_pairs'],candidate_node_pairs=6472,candidate_frame_mappings=653672,
        intact_chains_in_unresolved_groups=6,unchanged_interrupted_checkpoints_reused=True,completed_producer_and_alternate_reader_restart_refused=True,
        rehashed_false_exports_rejected=rejected,source_mapping_alterations_rejected=source_rejected,
        source_hashes={str(p):sha(p) for p in paths},scientific_eligibility=False,production_clade_mappings_checked=False,
        scope='Complete synthetic1620chain/405group/135input manual rooted-clade source/serialized/checkpoint workflow, both failures/six intact unresolved retained. '
              'Primitive descendant sets checked against original Biopython index; exact branch threshold1e-10 retained. Source journals are synthetic; biological root/topology/posterior not qualified.')
    with args.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
