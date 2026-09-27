#!/usr/bin/env python3
"""Prepare model-pair jobs while retaining every reviewed duplication event link."""
import argparse,csv,hashlib,json,sqlite3
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['review','review-plan','bridge','catalog','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    sources={str(a.review/'receipt.json'):sha(a.review/'receipt.json'),str(a.bridge):sha(a.bridge),str(a.catalog/'receipt.json'):sha(a.catalog/'receipt.json')}
    review=json.loads((a.review/'receipt.json').read_text());catalog=json.loads((a.catalog/'receipt.json').read_text())
    if review['status']!='complete_duplication_candidate_join_and_tree_review':raise ValueError('Incomplete candidate review')
    rp=json.loads(a.review_plan.read_text())
    if sha(a.review_plan)!=review['plan_sha256'] or Path(rp['bridge']).resolve()!=a.bridge.resolve() or sha(a.bridge)!=rp['pins'][rp['bridge']]:raise ValueError('Review bridge binding differs')
    sources[str(a.review_plan)]=sha(a.review_plan)
    bp=a.bridge.parent/'receipt.json';br=json.loads(bp.read_text());sources[str(bp)]=sha(bp)
    if br['artifacts'][a.bridge.name]!=sha(a.bridge) or br['catalog_receipt_sha256']!=sha(a.catalog/'receipt.json'):raise ValueError('Bridge/catalog binding differs')
    rows=[]
    for guide in ['profile','mafft']:
        path=a.review/(guide+'_candidate_tree_checks.tsv')
        if sha(path)!=review['artifacts'][path.name]:raise ValueError('Changed candidate checks')
        sources[str(path)]=sha(path)
        with path.open() as f:
            rows.extend(dict(r,guide=guide) for r in csv.DictReader(f,delimiter='\t'))
    needed={r[k] for r in rows for k in ['gene_a','gene_b']};links={}
    con=sqlite3.connect('file:'+str(a.bridge.resolve())+'?mode=ro',uri=True)
    for t,p,s,m,v,path in con.execute('SELECT taxon_id,protein_id,sequence_sha256,model_id,version,model_path FROM structures'):
        key=t+'_'+p
        if key in needed:
            if key in links:raise ValueError('Ambiguous protein')
            links[key]=(s,m,v,path)
    con.close()
    if set(links)!=needed:raise ValueError('Missing frozen protein link')
    models_needed={(v[1],v[2]) for v in links.values()};models={};path=a.catalog/'models.jsonl'
    if sha(path)!=catalog['artifacts']['models.jsonl']:raise ValueError('Changed frozen catalog')
    sources[str(path)]=sha(path)
    with path.open() as f:
        for line in f:
            m=json.loads(line);key=(m['model_id'],m['version'])
            if key in models_needed:
                if key in models:raise ValueError('Repeated model key')
                models[key]=m
    if set(models)!=models_needed:raise ValueError('Missing catalog model')
    for seq,m,v,path in links.values():
        if (models[m,v]['sequence_sha256'],models[m,v]['path'])!=(seq,path):raise ValueError('Catalog/bridge mismatch')
    a.output.mkdir(parents=True);pairs={};counts={'event_links':len(rows),'same_model_links':0,'unresolved_tree_links':0,'eligible_event_links':0}
    fields=list(rows[0])+['version_a','version_b','pair_key','comparison_status']
    with (a.output/'event_model_pair_links.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for r in rows:
            ma,mb=links[r['gene_a']],links[r['gene_b']]
            if (ma[1],mb[1])!=(r['model_a'],r['model_b']):raise ValueError('Candidate model identity differs')
            keys=sorted([(ma[1],ma[2]),(mb[1],mb[2])]);key=hashlib.sha256(json.dumps(keys,separators=(',',':')).encode()).hexdigest()
            if r['tree_status']!='exact_reported_pair':status='unresolved_tree';counts['unresolved_tree_links']+=1
            elif keys[0]==keys[1]:status='identical_model_no_alignment';counts['same_model_links']+=1
            else:
                status='queued_distinct_models';counts['eligible_event_links']+=1;pairs[key]=keys
            w.writerow(dict(r,version_a=ma[2],version_b=mb[2],pair_key=key,comparison_status=status))
    with (a.output/'model_pairs.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['pair_key','model_a','version_a','model_b','version_b'])
        for key,(ma,mb) in sorted(pairs.items()):w.writerow([key,*ma,*mb])
    with (a.output/'models.jsonl').open('w') as f:
        for key,m in sorted(models.items()):f.write(json.dumps(m,separators=(',',':'))+'\n')
    active_models={tuple(m) for pair in pairs.values() for m in pair};byte_count=sum(Path(models[m]['path']).stat().st_size for m in active_models)
    for path,h in sources.items():
        if sha(path)!=h:raise ValueError('Source changed during join')
    receipt=dict(status='complete_reviewed_duplication_model_pair_queue',source_hashes=sources,script_sha256=sha(__file__),counts=counts,unique_models=len(models),active_models=len(active_models),unique_distinct_model_pairs=len(pairs),directed_alignments_both_orders=2*len(pairs),active_coordinate_bytes=byte_count,artifacts={p.name:sha(p) for p in a.output.iterdir()},scope='Exact model/version and sequence join to frozen catalog with every gene/event link retained. Identical models explicitly excluded from redundant compute, not erased from event universe. Native catalog confidence retained without an eligibility threshold. Coordinate bytes are a size estimate; raw-coordinate validation, alignments, matched controls and duplication/asymmetry inference remain pending.')
    (a.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
