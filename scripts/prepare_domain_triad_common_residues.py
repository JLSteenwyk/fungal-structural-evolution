#!/usr/bin/env python3
"""Build common-reference and cycle-consistent residue maps for all domain triads/orders/masks."""
import argparse,csv,gzip,hashlib,itertools,json,sqlite3
from collections import Counter
from functools import lru_cache
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha


def digest(value):return hashlib.sha256(json.dumps(value,separators=(',',':')).encode()).hexdigest()


def common_residues(ab,ar,br):
    # Inputs: A -> B, reference -> A, reference -> B using ORIGINAL positions.
    if any(len(set(m.values()))!=len(m) for m in [ab,ar,br]):raise ValueError('Non-bijective pair mapping')
    common=[[ar[r],br[r],r] for r in sorted(set(ar)&set(br))]
    consistent=[t for t in common if ab.get(t[0])==t[1]]
    return common,consistent


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True)
    args=parser.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        if sha(args.plan)!=ph:raise ValueError('Changed plan')
        for path,h in plan['pins'].items():
            if sha(path)!=h:raise ValueError('Changed input: '+path)
    verify()
    integration=Path(plan['integration']);ir=json.loads((integration/'receipt.json').read_text());ia=json.loads(Path(plan['integration_readback']).read_text())
    if ia['status']!='passed_full_domain_triad_integration_readback' or ia['producer_receipt_sha256']!=sha(integration/'receipt.json'):raise ValueError('Unbound triad integration')
    if sha(integration/'domain_triads.sqlite')!=ir['artifacts']['domain_triads.sqlite']:raise ValueError('Changed triad database')
    gp=json.loads(Path(plan['geometry_plan']).read_text());geo=Path(gp['output']);gr=json.loads((geo/'receipt.json').read_text());ga=json.loads(Path(plan['geometry_readback']).read_text())
    if ga['status']!='passed_full_domain_geometry_readback' or ga['producer_receipt_sha256']!=sha(geo/'receipt.json'):raise ValueError('Unbound geometry audit')
    if sha(geo/'alignment_geometry.tsv')!=gr['artifacts']['alignment_geometry.tsv']:raise ValueError('Changed geometry table')
    geometry={(r['pair_key'],r['mask'],int(r['order'])):r['geometry_status'] for r in csv.DictReader((geo/'alignment_geometry.tsv').open(),delimiter='\t')}
    diagnostic=Path(gp['diagnostic']);dr=json.loads((diagnostic/'receipt.json').read_text())
    if sha(diagnostic/'receipt.json')!=gr['diagnostic_receipt_sha256']:raise ValueError('Changed diagnostic')
    if sha(diagnostic/'numeric_readback.tsv')!=dr['artifacts']['numeric_readback.tsv']:raise ValueError('Changed numeric diagnostic')
    numeric={(r['pair_key'],r['mask'],int(r['order'])):r['rmsd_status'] for r in csv.DictReader((diagnostic/'numeric_readback.tsv').open(),delimiter='\t')}
    dp=json.loads(Path(gp['diagnostic_plan']).read_text());sp=json.loads(Path(dp['source_plan']).read_text());root=Path(sp['output']);folder=Path(sp['inputs'])
    if sha(root/'receipt.json')!=dr['producer_receipt_sha256']:raise ValueError('Changed native producer')
    nr=json.loads((root/'receipt.json').read_text());manifest=root/'checkpoint_manifest.tsv'
    if sha(manifest)!=nr['artifacts']['checkpoint_manifest.tsv']:raise ValueError('Changed native checkpoint manifest')
    checkpoints={r['path']:r['sha256'] for r in csv.DictReader(manifest.open(),delimiter='\t')}
    inputs={}
    ip=json.loads((folder/'receipt.json').read_text())
    if sha(folder/'inputs.jsonl')!=ip['artifacts']['inputs.jsonl']:raise ValueError('Changed residue positions')
    for line in (folder/'inputs.jsonl').open():
        r=json.loads(line);inputs[r['interval_id'],r['mask']]={k:r[k] for k in ['original_positions','sequence','status','interval_length']}
    @lru_cache(maxsize=8192)
    def native(pair,mask,order):
        rel=f'pairs/{pair[:2]}/{pair}-{mask}-{order}.json';path=root/rel
        if sha(path)!=checkpoints[rel]:raise ValueError('Changed checkpoint')
        record=json.loads(path.read_text());key=pair,mask,order
        if (record['pair_key'],record['mask'],record['order'])!=key:raise ValueError('Wrong native key')
        ends=[r['interval_id'] for r in record['inputs']];rows=[inputs[e,mask] for e in ends]
        if record['status']=='input_unavailable':
            if all(r['status']=='ready' for r in rows) or key in numeric or key in geometry:raise ValueError('Invalid unavailable disposition')
            return ends,[],['input_unavailable']
        if record['status']!='aligned':raise ValueError('Unhandled native outcome')
        problems=[]
        if numeric[key]!='within_printed_rounding':problems.append('numeric_discrepancy')
        if geometry[key]!='unique_at_numeric_tolerance':problems.append('degenerate_geometry')
        strings=[record['metrics']['alignment_left'],record['metrics']['alignment_right']]
        if any(s.replace('-','')!=r['sequence'] for s,r in zip(strings,rows)):raise ValueError('Sequence/position mismatch')
        i=j=0;mapping=[]
        for a,b in zip(*strings):
            if a!='-' and b!='-':mapping.append((rows[0]['original_positions'][i],rows[1]['original_positions'][j]))
            i+=a!='-';j+=b!='-'
        if len(mapping)!=record['metrics']['aligned_length']:raise ValueError('Native aligned length mismatch')
        return ends,mapping,problems
    db=sqlite3.connect(f'file:{integration}/domain_triads.sqlite?mode=ro',uri=True);db.row_factory=sqlite3.Row
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    columns=['interval_a','interval_b','interval_reference','duplicate_domain_pair_key','a_reference_domain_pair_key','b_reference_domain_pair_key']
    triples={};links=0
    with (out/'event_domain_links.tsv').open('w') as f:
        fields=['triad_key','guide','family','gene_node','gene_a','gene_b','reference_gene','policy','boundary','pfam_accession']
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for r in db.execute('SELECT * FROM matched_domains'):
            triple=[r[k] for k in columns[:3]];key=digest(triple);value={c:r[c] for c in columns}
            if key in triples and triples[key]!=value:raise ValueError('Conflicting triad identity')
            triples[key]=value;w.writerow(dict(triad_key=key,**{f:r[f] for f in fields[1:]}));links+=1
    db.close()
    if links!=ir['counts']['matched_domains']:raise ValueError('Incomplete event/domain link inventory')
    counts=Counter();dispositions=0;common_total=cycle_total=0
    with gzip.open(out/'common_residue_maps.jsonl.gz','wt') as f,(out/'triads.jsonl').open('w') as tf:
        for index,(tk,t) in enumerate(sorted(triples.items()),1):
            ids=[t[c] for c in columns[:3]];edges=[(ids[0],ids[1]),(ids[2],ids[0]),(ids[2],ids[1])];pairs=[t[c] for c in columns[3:]]
            if any(digest(sorted(e))!=p for e,p in zip(edges,pairs)):raise ValueError('Domain pair identity mismatch')
            tf.write(json.dumps(dict(triad_key=tk,**t),separators=(',',':'))+'\n')
            for mask in ['full','plddt70']:
                lengths=[inputs[i,mask]['interval_length'] for i in ids]
                for orders in itertools.product([0,1],repeat=3):
                    maps=[];edge_problems=[]
                    for pair,order,endpoints in zip(pairs,orders,edges):
                        ends,mapping,problems=native(pair,mask,order)
                        if set(ends)!=set(endpoints):raise ValueError('Interval endpoints differ')
                        maps.append(dict(mapping if ends==list(endpoints) else [(y,x) for x,y in mapping]));edge_problems.append(problems)
                    common,consistent=common_residues(*maps)
                    problems=sorted(set(itertools.chain.from_iterable(edge_problems)))
                    if not common:problems.append('no_common_reference_residues')
                    elif not consistent:problems.append('no_cycle_consistent_residues')
                    r=dict(triad_key=tk,mask=mask,orders=list(orders),original_interval_lengths=lengths,edge_exclusions=edge_problems,exclusions=problems,reference_common_triples=common,cycle_consistent_triples=consistent,common_reference_count=len(common),cycle_consistent_count=len(consistent))
                    f.write(json.dumps(r,separators=(',',':'))+'\n');dispositions+=1;common_total+=len(common);cycle_total+=len(consistent)
                    counts[mask+':'+('recorded_without_edge_exclusions' if not problems else ';'.join(problems))]+=1
            if index%1000==0:print('Mapped triads',index,'/',len(triples),flush=True)
    verify()
    if dispositions!=len(triples)*16:raise ValueError('Incomplete mask/order combination universe')
    result=dict(status='complete_domain_common_residue_inventory_pending_readback',plan_sha256=ph,distinct_oriented_interval_triads=len(triples),event_domain_links=links,mask_order_dispositions=dispositions,common_reference_residue_occurrences=common_total,cycle_consistent_residue_occurrences=cycle_total,counts=dict(counts),artifacts={p.name:sha(p) for p in out.iterdir() if p.is_file()},scope='Every matched interval triad, both masks, all eight native-order combinations. Original A/B/reference residue positions retained; direct A-B mapping consistency reported separately from reference intersection. All numerical/geometry exclusions propagated. Shared computations retain every guide/event/reference/policy/boundary/Pfam link. No structural contrast, homology validation or biological inference.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
