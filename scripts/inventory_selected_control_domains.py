#!/usr/bin/env python3
"""Link every selected control to shared single-copy domain comparisons."""
import argparse,csv,gzip,json,hashlib
from pathlib import Path
from collections import Counter
from screen_duplication_domain_alignment_coverage import sha


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify();selection=Path(plan['selection']);sr=json.loads((selection/'receipt.json').read_text());audit=json.loads(Path(plan['selection_audit']).read_text())
    assert audit['status']=='passed_full_background_control_selection_readback' and audit['producer_receipt_sha256']==sha(selection/'receipt.json')
    assert sha(selection/'selections.tsv.gz')==sr['artifacts']['selections.tsv.gz']
    graph=Path(plan['graph']);gr=json.loads((graph/'receipt.json').read_text());ga=json.loads(Path(plan['graph_audit']).read_text());assert ga['producer_receipt_sha256']==sha(graph/'receipt.json')
    nodes={}
    for label in ['target','background']:
        path=graph/(label+'_nodes.jsonl');assert sha(path)==gr['artifacts'][path.name]
        node={}
        for line in path.open():
            r=json.loads(line);assert r['node_id'] not in node
            node[r['node_id']]={k:r[k] for k in ['pair_key','same_model','guide','family']}
        nodes[label]=node
    indexes={}
    for label,spec in plan['domains'].items():
        root=Path(spec['root']);receipt=json.loads((root/'receipt.json').read_text());proof=json.loads(Path(spec['audit']).read_text())
        assert proof['producer_receipt_sha256']==sha(root/'receipt.json')
        path=root/'domain_pair_links.tsv';assert sha(path)==receipt['artifacts'][path.name]
        index={}
        for row in csv.DictReader(path.open(),delimiter='\t'):
            if label=='target' and 'primary' not in row['scope'].split('+'):continue
            if label=='background':assert row['scope']=='background'
            k=row['pair_key'],row['policy'];d=index.setdefault(k,{})
            slot=row['pfam_accession'],row['boundary'];assert slot not in d
            d[slot]=row['domain_pair_key']
        indexes[label]=index
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);configs={};counts=Counter();rows=0;seen=set()
    with gzip.open(selection/'selections.tsv.gz','rt') as f,gzip.open(out/'selection_domain_links.tsv.gz','wt') as target:
        reader=csv.DictReader(f,delimiter='\t');w=csv.DictWriter(target,reader.fieldnames+['domain_config_id'],delimiter='\t',lineterminator='\n');w.writeheader()
        for row in reader:
            key=row['target_id'],row['policy'],row['scenario_id'];assert key not in seen;seen.add(key)
            t=nodes['target'][row['target_id']];b=nodes['background'][row['background_id']];assert t['guide']==b['guide'] and t['family']==b['family']
            identity=[t['pair_key'],b['pair_key'],row['policy']];cid=hashlib.sha256(json.dumps(identity,separators=(',',':')).encode()).hexdigest()
            if cid not in configs:
                ti=indexes['target'].get((t['pair_key'],row['policy']),{});bi=indexes['background'].get((b['pair_key'],row['policy']),{})
                tf={x for x,y in ti};bf={x for x,y in bi};shared=sorted(tf&bf)
                for index,families in [(ti,tf),(bi,bf)]:assert set(index)=={(x,y) for x in families for y in ['alignment','envelope']}
                matches=[dict(pfam_accession=pf,boundary=boundary,target_domain_pair=ti[pf,boundary],background_domain_pair=bi[pf,boundary]) for pf in shared for boundary in ['alignment','envelope']]
                status='identical_model_requires_separate_handling' if t['same_model'] or b['same_model'] else 'shared_domain_comparisons' if shared else 'no_shared_domain_comparison'
                configs[cid]=dict(domain_config_id=cid,target_pair_key=t['pair_key'],background_pair_key=b['pair_key'],policy=row['policy'],target_same_model=t['same_model'],background_same_model=b['same_model'],target_accessions=sorted(tf),background_accessions=sorted(bf),shared_accessions=shared,matches=matches,status=status)
            config=configs[cid];assert config['target_same_model']==t['same_model'] and config['background_same_model']==b['same_model']
            counts[t['guide']+'|'+row['policy']+'|'+row['scenario_id']+'|'+config['status']]+=1
            w.writerow(dict(row,domain_config_id=cid));rows+=1
            if rows%500000==0:print('Linked selections',rows,'/',sr['selected_records'],flush=True)
    assert rows==sr['selected_records']==audit['selected_records']
    with (out/'domain_configurations.jsonl').open('w') as f:
        for key in sorted(configs):f.write(json.dumps(configs[key],separators=(',',':'))+'\n')
    verify();r=dict(status='complete_selected_control_domain_inventory_pending_readback',plan_sha256=ph,selected_records=rows,distinct_configurations=len(configs),configuration_status_counts=dict(Counter(x['status'] for x in configs.values())),selected_status_counts=dict(counts),selection_receipt_sha256=sha(selection/'receipt.json'),artifacts={p.name:sha(p) for p in out.iterdir()},scope='All selected metadata controls preserved and linked to shared Pfam accessions in duplicate-pair and background domain inventories. Both boundary definitions retained. Identical models require separate handling; no shared comparisons is not biological domain absence. No structural scores, coverage qualification or effect inference.')
    (out/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k!='selected_status_counts'},indent=2))

if __name__=='__main__':main()
