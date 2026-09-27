#!/usr/bin/env python3
"""Materialize all additional background models after full coordinate readback."""
import argparse,gzip,hashlib,json,shutil,time
from collections import Counter
from pathlib import Path
import psutil
from duplication_alignment_inputs import render_ca
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        if sha(a.plan)!=ph:raise ValueError('Changed materialization plan')
        for path,h in plan['pins'].items():
            if sha(path)!=h:raise ValueError('Changed pin '+path)
    verify();dependency=plan['producer']
    while psutil.pid_exists(dependency['pid']):
        try:
            proc=psutil.Process(dependency['pid'])
            if abs(proc.create_time()-dependency['created'])>.01 or proc.status()==psutil.STATUS_ZOMBIE:break
            if proc.cmdline()!=dependency['cmdline']:raise ValueError('Dependency identity changed')
        except psutil.NoSuchProcess:break
        time.sleep(30)
    verify();auditdir=Path(plan['readback']);ap=auditdir/'receipt.json';ar=json.loads(ap.read_text());audit_hash=sha(ap)
    if ar['status']!='passed_background_exported_ca_readback' or ar['plan_sha256']!=sha(plan['readback_plan']):raise ValueError('Full coordinate readback required')
    source=Path(plan['coordinates']);srp=source/'receipt.json';sr=json.loads(srp.read_text());source_hash=sha(srp)
    if sr['status']!='complete_background_coordinate_validation_with_dispositions' or source_hash!=ar['producer_receipt_sha256']:raise ValueError('Source/audit binding differs')
    sourceplan=json.loads(Path(plan['coordinate_plan']).read_text())
    readbackplan=json.loads(Path(plan['readback_plan']).read_text())
    inventory=Path(sourceplan['inventory'])
    if sr['plan_sha256']!=sha(plan['coordinate_plan']):raise ValueError('Coordinate plan differs')
    if sr['source_inventory_receipt_sha256']!=sha(inventory/'receipt.json'):raise ValueError('Inventory binding differs')
    if Path(readbackplan['source_plan']).resolve()!=Path(plan['coordinate_plan']).resolve() or Path(readbackplan['models']).resolve()!=(inventory/'additional_models.jsonl').resolve():raise ValueError('Readback model source differs')
    if Path(readbackplan['source']).resolve()!=source.resolve() or Path(sourceplan['output']).resolve()!=source.resolve():raise ValueError('Coordinate output differs')
    inventory_receipt=json.loads((inventory/'receipt.json').read_text())
    inventory_audit=json.loads(Path(sourceplan['readback']).read_text())
    if inventory_receipt['status']!='complete_background_measurement_inventory_pending_readback' or inventory_audit['status']!='passed_full_background_measurement_inventory_readback' or inventory_audit['producer_receipt_sha256']!=sha(inventory/'receipt.json'):raise ValueError('Unverified background inventory')
    modelpath=inventory/'additional_models.jsonl'
    if sha(modelpath)!=inventory_receipt['artifacts'][modelpath.name]:raise ValueError('Changed background model list')
    with modelpath.open() as handle:models=[json.loads(line) for line in handle]
    if len(models)!=inventory_audit['additional_models'] or len(models)!=inventory_receipt['additional_models']:raise ValueError('Background model partition differs')
    expected={(m['model_id'],m['version']):m for m in models}
    if sr['models']!=len(expected) or ar['models']!=len(expected):raise ValueError('Additional model universe differs')
    out=Path(plan['output'])
    if out.exists():raise FileExistsError('Fresh materialization output required')
    out.mkdir(parents=True);seen=set();counts=Counter();written_bytes=0
    with (out/'inputs.jsonl').open('w') as manifest:
        for index,shard in enumerate(sr['shards']):
            if shutil.disk_usage(out).free<plan['minimum_free_disk_gib']*2**30:raise RuntimeError('Disk reserve reached')
            path=source/shard['output'];proofpath=auditdir/(path.name+'.readback.json');proof=json.loads(proofpath.read_text())
            if sha(proofpath)!=ar['proofs'][proofpath.name] or proof['source_sha256']!=shard['output_sha256'] or sha(path)!=shard['output_sha256']:raise ValueError('Unaudited or changed coordinate shard')
            with gzip.open(path,'rt') as handle:
                for line in handle:
                    r=json.loads(line);key=(r['model_id'],r['version'])
                    if key in seen:raise ValueError('Repeated model')
                    if key not in expected:raise ValueError('Unexpected additional model')
                    m=expected[key]
                    if any(r[k]!=m[v] for k,v in [('source_sha256','sha256'),('source_path','path'),('sequence_sha256','sequence_sha256')]):raise ValueError('Model provenance differs')
                    seen.add(key);name=hashlib.sha256(json.dumps(key,separators=(',',':')).encode()).hexdigest()
                    for label,threshold in [('full',None),('plddt70',70)]:
                        row=dict(model_id=key[0],version=key[1],mask=label,coordinate_shard=path.name,source_sha256=r['source_sha256'])
                        if r['status']!='validated':row.update(status='source_rejected',reason=r['reason'])
                        else:
                            blob,sequence,positions=render_ca(r,threshold);row.update(sequence=sequence,original_positions=positions,retained_residues=len(sequence),original_length=r['length'])
                            if len(sequence)<3:row['status']='too_few_retained_residues'
                            else:
                                target=out/'pdb'/name[:2]/(name+'-'+label+'.pdb');target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(blob);digest=hashlib.sha256(blob).hexdigest()
                                if sha(target)!=digest:raise ValueError('Written coordinate bytes differ')
                                row.update(status='ready',path=str(target),sha256=digest);written_bytes+=len(blob)
                        counts[label+':'+row['status']]+=1;manifest.write(json.dumps(row,separators=(',',':'))+'\n')
            print('Materialized shards',index+1,'/',len(sr['shards']),flush=True)
    if seen!=set(expected):raise ValueError('Materialization model count differs')
    verify()
    if sha(ap)!=audit_hash or sha(srp)!=source_hash:raise ValueError('Upstream receipt changed')
    result=dict(status='complete_background_alignment_input_materialization',plan_sha256=ph,readback_receipt_sha256=audit_hash,coordinate_receipt_sha256=source_hash,inventory_receipt_sha256=sha(inventory/'receipt.json'),models=len(seen),input_dispositions=sum(counts.values()),counts=dict(counts),pdb_bytes=written_bytes,artifacts={'inputs.jsonl':sha(out/'inputs.jsonl')},script_sha256=sha(__file__),scope='Full and pLDDT>=70 C-alpha inputs for the exact additional-model partition of the frozen background inventory. Original residue positions and masked sequence retained; short masks and source rejections explicit. PDB coordinates rounded to 0.001 A with bounds checked, confidence to 0.01. No PAE qualification; masked TM scores normalize to retained length and cannot be interpreted interchangeably with full-model scores. Pair alignments, domain/orientation controls and biological effects pending.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
