"""Bind the new background pair workload to three audited input collections."""
import csv,hashlib,json
from collections import Counter
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha


def load_handoff(plan):
    inventory=Path(plan['inventory']);qr=json.loads((inventory/'receipt.json').read_text());proof=json.loads(Path(plan['inventory_readback']).read_text())
    if qr['status']!='complete_background_measurement_inventory_pending_readback' or proof['status']!='passed_full_background_measurement_inventory_readback' or proof['producer_receipt_sha256']!=sha(inventory/'receipt.json'):raise ValueError('Unverified background inventory')
    for name in ['active_models.jsonl','model_pairs.tsv']:
        if sha(inventory/name)!=qr['artifacts'][name]:raise ValueError('Changed background workload')
    wanted={}
    with (inventory/'active_models.jsonl').open() as f:
        for line in f:
            m=json.loads(line);key=(m['model_id'],m['version'])
            if key in wanted:raise ValueError('Repeated active model')
            wanted[key]=m['sha256']
    if len(wanted)!=proof['active_models']:raise ValueError('Active model count differs')
    inputs={};bindings={};owners=set()
    for label,spec in plan['input_sources'].items():
        folder=Path(spec['inputs']);rp=folder/'receipt.json';r=json.loads(rp.read_text());materialization=json.loads(Path(spec['input_plan']).read_text());source=Path(spec['model_inventory']);mr=json.loads((source/'receipt.json').read_text())
        if r['status']!=spec['status'] or r['plan_sha256']!=sha(spec['input_plan']) or Path(materialization['output']).resolve()!=folder.resolve():raise ValueError('Unbound input collection')
        if r[spec['inventory_receipt_field']]!=sha(source/'receipt.json'):raise ValueError('Input inventory differs')
        modelpath=source/spec['model_file']
        if sha(modelpath)!=mr['artifacts'][spec['model_file']]:raise ValueError('Changed source model partition')
        expected={}
        with modelpath.open() as f:
            for line in f:
                m=json.loads(line);key=(m['model_id'],m['version'])
                if key in expected or key in owners:raise ValueError('Overlapping input model partitions')
                expected[key]=m['sha256'];owners.add(key)
                if key in wanted and wanted[key]!=m['sha256']:raise ValueError('Active model provenance differs')
        auditpath=Path(materialization['readback'])/'receipt.json';ar=json.loads(auditpath.read_text());coordinatepath=Path(materialization['coordinates'])/'receipt.json'
        if ar['status']!=spec['coordinate_readback_status'] or r['readback_receipt_sha256']!=sha(auditpath) or r['coordinate_receipt_sha256']!=sha(coordinatepath) or ar['producer_receipt_sha256']!=sha(coordinatepath):raise ValueError('Coordinate proof binding differs')
        manifest=folder/'inputs.jsonl'
        if sha(manifest)!=r['artifacts']['inputs.jsonl']:raise ValueError('Changed materialized inputs')
        seen=set();counts=Counter()
        with manifest.open() as f:
            for line in f:
                row=json.loads(line);key=(row['model_id'],row['version'],row['mask'])
                if key[:2] not in expected or key in seen or key[2] not in ['full','plddt70']:raise ValueError('Unexpected or repeated input')
                if row['source_sha256']!=expected[key[:2]]:raise ValueError('Wrong input source bytes')
                if row['status'] not in ['ready','source_rejected','too_few_retained_residues']:raise ValueError('Unknown disposition')
                if row['status']=='ready' and not all(k in row for k in ['path','sha256','sequence']):raise ValueError('Incomplete ready input')
                counts[key[2]+':'+row['status']]+=1;seen.add(key)
                if key[:2] in wanted:inputs[key]={k:row[k] for k in ['status','path','sha256','sequence'] if k in row}
        if len(seen)!=2*len(expected) or r['models']!=len(expected) or r['input_dispositions']!=len(seen) or r['counts']!=dict(counts):raise ValueError('Incomplete input grid')
        for path in [rp,manifest,auditpath,coordinatepath,modelpath,source/'receipt.json']:bindings[str(path)]=sha(path)
    if len(inputs)!=2*len(wanted) or {k[:2] for k in inputs}!=set(wanted):raise ValueError('Incomplete active input set')
    previous={}
    for spec in plan['existing_queues']:
        root=Path(spec['path']);r=json.loads((root/'receipt.json').read_text());path=root/'model_pairs.tsv'
        if sha(path)!=r['artifacts']['model_pairs.tsv']:raise ValueError('Changed existing pairs')
        with path.open() as f:
            for row in csv.DictReader(f,delimiter='\t'):previous.setdefault(row['pair_key'],spec['label'])
        bindings[str(path)]=sha(path);bindings[str(root/'receipt.json')]=sha(root/'receipt.json')
    pairs=[];seen=set();reused=0
    with (inventory/'model_pairs.tsv').open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            ends=sorted([(row['model_a'],int(row['version_a'])),(row['model_b'],int(row['version_b']))]);key=hashlib.sha256(json.dumps(ends,separators=(',',':')).encode()).hexdigest()
            if key!=row['pair_key'] or key in seen or ends[0]==ends[1]:raise ValueError('Invalid pair')
            seen.add(key)
            if any((*end,mask) not in inputs for end in ends for mask in ['full','plddt70']):raise ValueError('Missing pair input')
            if row['work_disposition']!=previous.get(key,'new_model_pair'):raise ValueError('Pair reuse partition differs')
            if key in previous:reused+=1
            else:pairs.append(row)
    if len(seen)!=qr['distinct_eligible_model_pairs'] or len(pairs)!=qr['new_model_pairs'] or len(pairs)!=plan['resources']['new_model_pairs']:raise ValueError('Background workload differs')
    for path in [inventory/'receipt.json',inventory/'model_pairs.tsv',inventory/'active_models.jsonl',Path(plan['inventory_readback'])]:bindings[str(path)]=sha(path)
    return inputs,pairs,bindings,hashlib.sha256(json.dumps(bindings,sort_keys=True).encode()).hexdigest()
