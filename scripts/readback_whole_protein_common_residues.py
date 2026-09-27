#!/usr/bin/env python3
"""Independently reconstruct every saved whole-protein common-residue mapping."""
import argparse,csv,gzip,hashlib,itertools,json,subprocess,time
import psutil
from collections import Counter
from functools import lru_cache
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha
def common_residues(ab,ar,br):
    assert all(len(set(m.values()))==len(m) for m in [ab,ar,br])
    common=sorted([(aa,bb,ref) for ref,aa in ar.items() for bb in ([br[ref]] if ref in br else [])],key=lambda x:x[2])
    relation=set(ab.items())
    return [list(t) for t in common],[list(t) for t in common if (t[0],t[1]) in relation]



def digest(value):return hashlib.sha256(json.dumps(value,separators=(',',':')).encode()).hexdigest()


def table(path):
    with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--launch',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        assert sha(args.plan)==ph
        for path,h in plan['pins'].items():assert sha(path)==h,path
    verify();launch=json.loads(args.launch.read_text());lh=sha(args.launch)
    assert launch['plan_sha256']==ph
    while True:
        try:
            p=psutil.Process(launch['pid'])
            if p.create_time()!=launch['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:break
        print('waiting_for_complete_whole_protein_residue_maps',launch['pid'],flush=True);time.sleep(30)
    terminal=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert terminal==dict(ActiveState='inactive',Result='success',ExecMainStatus='0'),terminal
    assert sha(args.launch)==lh
    produced=Path(plan['output']);rp=produced/'receipt.json';pr=json.loads(rp.read_text());rh=sha(rp)
    assert pr['status']=='complete_whole_protein_common_residue_inventory_pending_readback' and pr['plan_sha256']==ph
    for name,h in pr['artifacts'].items():assert sha(produced/name)==h
    verify();root=Path(plan['triads']);receipt=json.loads((root/'receipt.json').read_text());audit=json.loads(Path(plan['triad_audit']).read_text())
    assert audit['status']=='passed_full_whole_protein_reference_triad_independent_join' and audit['source_receipt_sha256']==sha(root/'receipt.json')
    for name,h in receipt['artifacts'].items():assert sha(root/name)==h
    triads=table(root/'model_triads.tsv');assert len(triads)==17619
    inputs={}
    for folder in plan['inputs']:
        folder=Path(folder);r=json.loads((folder/'receipt.json').read_text());assert sha(folder/'inputs.jsonl')==r['artifacts']['inputs.jsonl']
        for line in (folder/'inputs.jsonl').open():
            row=json.loads(line);key=row['model_id'],int(row['version']),row['mask'];assert key not in inputs
            assert len(row['sequence'])==len(row['original_positions'])==row['retained_residues']
            assert len(set(row['original_positions']))==len(row['original_positions'])
            inputs[key]=row
    statuses={};hashes={};roots={}
    for label,spec in plan['sources'].items():
        orders=Path(spec['orders']);r=json.loads((orders/'receipt.json').read_text());a=json.loads(Path(spec['audit']).read_text())
        assert a['status']==f'passed_full_{label}_usable_order_summary_readback' and a['source_receipt_sha256']==sha(orders/'receipt.json')
        assert sha(orders/'pair_mask_order_summary.tsv')==r['artifacts']['pair_mask_order_summary.tsv']
        native=Path(spec['native']);nr=json.loads((native/'receipt.json').read_text());assert sha(native/'receipt.json')==r['native_receipt_sha256']
        assert sha(native/'checkpoint_manifest.tsv')==nr['artifacts']['checkpoint_manifest.tsv']
        roots[label]=native;hashes[label]={row['path']:row['sha256'] for row in table(native/'checkpoint_manifest.tsv')}
        for row in table(orders/'pair_mask_order_summary.tsv'):
            for order in [0,1]:
                key=label,row['pair_key'],row['mask'],order;assert key not in statuses
                statuses[key]=(row[f'order{order}_status'],row[f'order{order}_native_status'],row[f'order{order}_numerical_exclusion_reasons'])
    @lru_cache(maxsize=4096)
    def native_map(label,pair,mask,order):
        rel=f'pairs/{pair[:2]}/{pair}-{mask}-{order}.json';path=roots[label]/rel;assert sha(path)==hashes[label][rel]
        r=json.loads(path.read_text());assert (r['pair_key'],r['mask'],r['order'])==(pair,mask,order)
        ends=[(v['model_id'],int(v['version'])) for v in r['inputs']];rows=[inputs[(*end,mask)] for end in ends]
        status,raw,why=statuses[label,pair,mask,order];assert r['status']==raw
        if status=='input_unavailable':
            assert any(v['status']!='ready' for v in rows)
            return ends,[],['input_unavailable'],status
        assert raw=='aligned' and status in ['aligned','excluded_numerically']
        assert (status=='excluded_numerically')==bool(why)
        strings=[r['metrics']['alignment_left'],r['metrics']['alignment_right']]
        assert len(strings[0])==len(strings[1])
        assert all(s.replace('-','')==v['sequence'] for s,v in zip(strings,rows))
        columns=[]
        for text,row in zip(strings,rows):
            positions=iter(row['original_positions'])
            columns.append([None if char=='-' else next(positions) for char in text])
            assert next(positions,None) is None
        pairs=[(a,b) for a,b in zip(*columns) if a is not None and b is not None]
        assert len(pairs)==r['metrics']['aligned_length']
        return ends,pairs,why.split(';') if why else [],status
    out=produced;counts=Counter();n=0;common_count=cycle_count=0;seen=set();nonempty_common=nonempty_cycle=disagree=0
    with gzip.open(out/'common_residue_maps.jsonl.gz','rt') as handle:
        for ix,t in enumerate(triads,1):
            ids=[(t[role+'_model'],int(t[role+'_version'])) for role in ['a','b','reference']]
            assert digest(ids)==t['triad_id'] and len(set(ids))==int(t['distinct_models'])
            for role,model in zip(['a','b','reference'],ids):
                full=inputs[(*model,'full')];assert hashlib.sha256(full['sequence'].encode()).hexdigest()==t[role+'_sequence_sha256']
            edges=[(ids[0],ids[1]),(ids[2],ids[0]),(ids[2],ids[1])]
            labels=['ab','a_reference','b_reference']
            for mask in ['full','plddt70']:
                for orders in itertools.product([0,1],repeat=3):
                    maps=[];exclusions=[];edge_statuses=[]
                    for label,ends,order in zip(labels,edges,orders):
                        pair=t[label+'_pair_key'];assert digest(sorted(ends))==pair
                        source=t[label+'_source']
                        if source=='identical_model':
                            assert ends[0]==ends[1];row=inputs[(*ends[0],mask)]
                            ready=row['status']=='ready';maps.append({p:p for p in row['original_positions']} if ready else {})
                            exclusions.append([] if ready else ['input_unavailable']);edge_statuses.append('identity_map_not_native_fit' if ready else 'identity_input_unavailable')
                        else:
                            source_label={'primary_pair':'primary','additional_reference_pair':'reference'}[source]
                            actual,pairs,why,status=native_map(source_label,pair,mask,order)
                            assert set(actual)==set(ends)
                            mapping=dict(pairs if actual==list(ends) else [(b,a) for a,b in pairs]);assert len(mapping)==len(pairs)
                            maps.append(mapping);exclusions.append(why);edge_statuses.append(status)
                    common,cycle=common_residues(*maps)
                    problems=sorted(set(itertools.chain.from_iterable(exclusions)))
                    if not common:problems.append('no_common_reference_residues')
                    elif not cycle:problems.append('no_cycle_consistent_residues')
                    value=dict(triad_id=t['triad_id'],distinct_models=int(t['distinct_models']),mask=mask,orders=list(orders),edge_order=['ab','reference_to_a','reference_to_b'],edge_statuses=edge_statuses,edge_exclusions=exclusions,exclusions=problems,original_lengths=[inputs[(*model,mask)]['original_length'] for model in ids],reference_common_triples=common,cycle_consistent_triples=cycle,common_reference_count=len(common),cycle_consistent_count=len(cycle))
                    actual=json.loads(next(handle));assert actual==value,(t['triad_id'],mask,orders)
                    key=t['triad_id'],mask,orders;assert key not in seen;seen.add(key)
                    n+=1;common_count+=len(common);cycle_count+=len(cycle)
                    nonempty_common+=bool(common);nonempty_cycle+=bool(cycle);disagree+=common!=cycle
                    counts[mask+':'+(';'.join(problems) if problems else 'recorded_without_edge_exclusions')]+=1
            if ix%500==0:print('Whole-protein triads independently checked',ix,'/ 17619',flush=True)
        assert next(handle,None) is None
    assert n==pr['mask_order_dispositions']==281904
    assert seen=={(t['triad_id'],m,o) for t in triads for m in ['full','plddt70'] for o in itertools.product([0,1],repeat=3)}
    assert dict(counts)==pr['counts']
    assert common_count==pr['common_reference_residue_occurrences'] and cycle_count==pr['cycle_consistent_residue_occurrences']
    for name in ['model_triads.tsv','event_reference_triads.tsv']:
        assert sha(out/name)==receipt['artifacts'][name]
    verify()
    assert sha(rp)==rh and sha(args.launch)==lh
    for name,h in pr['artifacts'].items():assert sha(produced/name)==h
    result=dict(status='passed_full_whole_protein_common_residue_readback',producer_receipt_sha256=rh,plan_sha256=ph,producer_terminal_state=terminal,distinct_oriented_model_triads=len(triads),event_reference_links=receipt['event_reference_links'],mask_order_dispositions=n,common_reference_residue_occurrences=common_count,cycle_consistent_residue_occurrences=cycle_count,dispositions_with_nonempty_reference_intersection=nonempty_common,dispositions_with_nonempty_consistent_mapping=nonempty_cycle,dispositions_with_any_mapping_disagreement=disagree,counts=dict(counts),script_sha256=sha(__file__),scope='Every saved triad/mask/order record reconstructed with iterator-based native position mapping and relation intersections. Exact identities, all residue triples, explicit identity maps, exclusion reasons, complete key universe, and event-link hashes checked. No fit or biological inference.')
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
