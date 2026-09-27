#!/usr/bin/env python3
"""Map all oriented whole-protein triads, masks and orders with explicit exclusions."""
import argparse,csv,gzip,hashlib,itertools,json
from collections import Counter
from functools import lru_cache
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha
from prepare_domain_triad_common_residues import common_residues


def digest(value):return hashlib.sha256(json.dumps(value,separators=(',',':')).encode()).hexdigest()


def table(path):
    with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);args=ap.parse_args()
    plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        assert sha(args.plan)==ph
        for path,h in plan['pins'].items():assert sha(path)==h,path
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
        i=j=0;pairs=[]
        for a,b in zip(*strings):
            if a!='-' and b!='-':pairs.append((rows[0]['original_positions'][i],rows[1]['original_positions'][j]))
            i+=a!='-';j+=b!='-'
        assert len(pairs)==r['metrics']['aligned_length']
        return ends,pairs,why.split(';') if why else [],status
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);counts=Counter();n=0;common_count=cycle_count=0
    with gzip.open(out/'common_residue_maps.jsonl.gz','wt') as handle:
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
                    handle.write(json.dumps(value,separators=(',',':'))+'\n');n+=1;common_count+=len(common);cycle_count+=len(cycle)
                    counts[mask+':'+(';'.join(problems) if problems else 'recorded_without_edge_exclusions')]+=1
            if ix%500==0:print('Whole-protein triads mapped',ix,'/ 17619',flush=True)
    assert n==281904
    for name in ['model_triads.tsv','event_reference_triads.tsv']:
        (out/name).write_bytes((root/name).read_bytes());assert sha(out/name)==receipt['artifacts'][name]
    verify()
    result=dict(status='complete_whole_protein_common_residue_inventory_pending_readback',plan_sha256=ph,distinct_oriented_model_triads=len(triads),event_reference_links=receipt['event_reference_links'],mask_order_dispositions=n,common_reference_residue_occurrences=common_count,cycle_consistent_residue_occurrences=cycle_count,counts=dict(counts),source_inventory_receipt_sha256=sha(root/'receipt.json'),artifacts={p.name:sha(p) for p in out.iterdir()},scope='All 17619 oriented model triads, two masks and eight native-order combinations retained. Original residue positions, common-reference and direct-AB-consistent subsets recorded separately. Identity models use explicit identity correspondence without fabricated native fits. Numerical/missing-input exclusions and 36944 event/reference links preserved. Shared models, tied references and repeated orders are not independent observations. No common-coordinate fit, ancestral state, structural asymmetry or biological inference yet.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
