#!/usr/bin/env python3
"""Independently reconstruct all common-reference and three-edge residue correspondences."""
import argparse,csv,gzip,hashlib,itertools,json,sqlite3
from collections import Counter
from functools import lru_cache
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        if sha(a.plan)!=ph:raise ValueError('Changed plan')
        for p,h in plan['pins'].items():
            if sha(p)!=h:raise ValueError('Changed source: '+p)
    verify();root=Path(plan['output']);rp=root/'receipt.json';r=json.loads(rp.read_text());rh=sha(rp)
    if r['status']!='complete_domain_common_residue_inventory_pending_readback' or r['plan_sha256']!=ph:raise ValueError('Wrong inventory receipt')
    for f,h in r['artifacts'].items():
        if sha(root/f)!=h:raise ValueError('Changed inventory artifact')
    gp=json.loads(Path(plan['geometry_plan']).read_text());dp=json.loads(Path(gp['diagnostic_plan']).read_text());sp=json.loads(Path(dp['source_plan']).read_text());native_root=Path(sp['output'])
    geometry={(x['pair_key'],x['mask'],int(x['order'])):x['geometry_status'] for x in csv.DictReader((Path(gp['output'])/'alignment_geometry.tsv').open(),delimiter='\t')}
    numeric={(x['pair_key'],x['mask'],int(x['order'])):x['rmsd_status'] for x in csv.DictReader((Path(gp['diagnostic'])/'numeric_readback.tsv').open(),delimiter='\t')}
    proofs={x['path']:x['sha256'] for x in csv.DictReader((native_root/'checkpoint_manifest.tsv').open(),delimiter='\t')}
    inputs={}
    for line in (Path(sp['inputs'])/'inputs.jsonl').open():
        x=json.loads(line);inputs[x['interval_id'],x['mask']]={k:x[k] for k in ['original_positions','sequence','status','interval_length']}
    triads={}
    for line in (root/'triads.jsonl').open():
        x=json.loads(line);k=x['triad_key'];ids=[x[f] for f in ['interval_a','interval_b','interval_reference']]
        if k in triads or hashlib.sha256(json.dumps(ids,separators=(',',':')).encode()).hexdigest()!=k:raise ValueError('Invalid triad ID')
        triads[k]=x
    db=sqlite3.connect(f"file:{plan['integration']}/domain_triads.sqlite?mode=ro",uri=True);db.row_factory=sqlite3.Row
    linkfields=['guide','family','gene_node','gene_a','gene_b','reference_gene','policy','boundary','pfam_accession'];expected_links=set();expected_triads=set()
    for x in db.execute('SELECT * FROM matched_domains'):
        ids=[x[f] for f in ['interval_a','interval_b','interval_reference']];k=hashlib.sha256(json.dumps(ids,separators=(',',':')).encode()).hexdigest();t=triads[k]
        for f in ['interval_a','interval_b','interval_reference','duplicate_domain_pair_key','a_reference_domain_pair_key','b_reference_domain_pair_key']:
            if t[f]!=x[f]:raise ValueError('Triad source mismatch')
        expected_links.add((k,)+tuple(x[f] for f in linkfields));expected_triads.add(k)
    db.close();actual_links=set();links=0
    for x in csv.DictReader((root/'event_domain_links.tsv').open(),delimiter='\t'):
        k=(x['triad_key'],)+tuple(x[f] for f in linkfields)
        if k in actual_links:raise ValueError('Repeated event link')
        actual_links.add(k);links+=1
    if actual_links!=expected_links or expected_triads!=set(triads) or links!=r['event_domain_links']:raise ValueError('Incomplete event/triad binding')
    @lru_cache(maxsize=8192)
    def pairs(pair,mask,order):
        rel=f'pairs/{pair[:2]}/{pair}-{mask}-{order}.json';p=native_root/rel
        if sha(p)!=proofs[rel]:raise ValueError('Changed native checkpoint')
        x=json.loads(p.read_text());ends=[v['interval_id'] for v in x['inputs']];key=pair,mask,order
        if x['status']=='input_unavailable':
            if all(inputs[e,mask]['status']=='ready' for e in ends):raise ValueError('Invalid unavailable input')
            return ends,set(),['input_unavailable']
        if x['status']!='aligned':raise ValueError('Unexpected native status')
        problems=[]
        if numeric[key]!='within_printed_rounding':problems.append('numeric_discrepancy')
        if geometry[key]!='unique_at_numeric_tolerance':problems.append('degenerate_geometry')
        strings=[x['metrics']['alignment_left'],x['metrics']['alignment_right']];columns=[]
        for text,end in zip(strings,ends):
            inp=inputs[end,mask]
            if text.replace('-','')!=inp['sequence']:raise ValueError('Sequence mapping differs')
            positions=iter(inp['original_positions']);columns.append([None if letter=='-' else next(positions) for letter in text])
            if next(positions,None) is not None:raise ValueError('Unconsumed positions')
        mapping={(left,right) for left,right in zip(*columns) if left is not None and right is not None}
        if len(mapping)!=x['metrics']['aligned_length']:raise ValueError('Native pair count differs')
        return ends,mapping,problems
    seen=set();counts=Counter();common_total=cycle_total=0;nonempty_common=nonempty_cycle=disagree=0
    with gzip.open(root/'common_residue_maps.jsonl.gz','rt') as f:
        for line in f:
            x=json.loads(line);t=triads[x['triad_key']];mask=x['mask'];orders=tuple(x['orders']);key=x['triad_key'],mask,orders
            if key in seen or mask not in ['full','plddt70'] or orders not in set(itertools.product([0,1],repeat=3)):raise ValueError('Unexpected/duplicate disposition')
            seen.add(key);ids=[t[c] for c in ['interval_a','interval_b','interval_reference']]
            edges=[(ids[0],ids[1]),(ids[0],ids[2]),(ids[1],ids[2])];relations=[];problems=[]
            for edge,field,order in zip(edges,['duplicate_domain_pair_key','a_reference_domain_pair_key','b_reference_domain_pair_key'],orders):
                ends,mapping,exclusions=pairs(t[field],mask,order)
                if set(ends)!=set(edge):raise ValueError('Endpoint mismatch')
                relations.append(mapping if ends==list(edge) else {(b,a) for a,b in mapping});problems.append(exclusions)
            ab,ar,br=relations;ref_a={ref:aa for aa,ref in ar};ref_b={ref:bb for bb,ref in br}
            common=sorted([(aa,ref_b[ref],ref) for ref,aa in ref_a.items() if ref in ref_b],key=lambda z:z[2])
            cycle=[z for z in common if (z[0],z[1]) in ab]
            expected_problem=sorted(set(itertools.chain.from_iterable(problems)))
            if not common:expected_problem.append('no_common_reference_residues')
            elif not cycle:expected_problem.append('no_cycle_consistent_residues')
            wanted=dict(triad_key=x['triad_key'],mask=mask,orders=list(orders),original_interval_lengths=[inputs[e,mask]['interval_length'] for e in ids],edge_exclusions=problems,exclusions=expected_problem,reference_common_triples=[list(z) for z in common],cycle_consistent_triples=[list(z) for z in cycle],common_reference_count=len(common),cycle_consistent_count=len(cycle))
            if x!=wanted:raise ValueError('Common residue record mismatch: '+str(key))
            common_total+=len(common);cycle_total+=len(cycle);nonempty_common+=bool(common);nonempty_cycle+=bool(cycle);disagree+=common!=cycle
            counts[mask+':'+('recorded_without_edge_exclusions' if not expected_problem else ';'.join(expected_problem))]+=1
            if len(seen)%10000==0:print('Common-residue records checked',len(seen),'/',r['mask_order_dispositions'],flush=True)
    expected={(t,m,o) for t in triads for m in ['full','plddt70'] for o in itertools.product([0,1],repeat=3)}
    if seen!=expected or len(seen)!=r['mask_order_dispositions'] or dict(counts)!=r['counts']:raise ValueError('Incomplete dispositions')
    if common_total!=r['common_reference_residue_occurrences'] or cycle_total!=r['cycle_consistent_residue_occurrences']:raise ValueError('Wrong residue totals')
    verify()
    if sha(rp)!=rh:raise ValueError('Producer changed')
    for name,h in r['artifacts'].items():
        if sha(root/name)!=h:raise ValueError('Artifact changed during audit')
    result=dict(status='passed_full_domain_common_residue_readback',producer_receipt_sha256=rh,plan_sha256=ph,checker_sha256=sha(__file__),triads=len(triads),event_links=links,dispositions=len(seen),common_reference_residue_occurrences=common_total,cycle_consistent_residue_occurrences=cycle_total,dispositions_with_nonempty_reference_intersection=nonempty_common,dispositions_with_nonempty_consistent_mapping=nonempty_cycle,dispositions_with_any_mapping_disagreement=disagree,counts=dict(counts),scope='Every event link and triad identity checked against integration source; every residue triple and all exclusions reconstructed from hashed native maps using independent relation intersections. These are alignment hypotheses and repeated order/mask occurrences, not validated evolutionary homology or biological effect estimates.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
