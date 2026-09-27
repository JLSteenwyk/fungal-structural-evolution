#!/usr/bin/env python3
"""Extract every verified domain interval after both complete coordinate readbacks."""
import argparse,gzip,hashlib,json,shutil,time
from collections import Counter,defaultdict
from pathlib import Path
import psutil
from duplication_domain_inputs import render_interval
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        if sha(a.plan)!=ph:raise ValueError('Changed extraction plan')
        for path,h in plan['pins'].items():
            if sha(path)!=h:raise ValueError('Changed pin: '+path)
    verify()
    for source in plan['coordinate_sources']:
        dep=source['producer']
        while psutil.pid_exists(dep['pid']):
            try:
                proc=psutil.Process(dep['pid'])
                if abs(proc.create_time()-dep['created'])>.01 or proc.status()==psutil.STATUS_ZOMBIE:break
                if proc.cmdline()!=dep['cmdline']:raise ValueError('Coordinate audit identity changed')
            except psutil.NoSuchProcess:break
            time.sleep(30)
    verify();inventory=Path(plan['inventory']);receipt=json.loads((inventory/'receipt.json').read_text());audit=json.loads(Path(plan['inventory_readback']).read_text())
    if audit['status']!='passed_full_duplication_domain_pair_inventory_readback' or audit['producer_receipt_sha256']!=sha(inventory/'receipt.json'):raise ValueError('Domain inventory not independently verified')
    if sha(inventory/'intervals.jsonl')!=receipt['artifacts']['intervals.jsonl']:raise ValueError('Changed intervals')
    intervals={};by_model=defaultdict(list)
    with (inventory/'intervals.jsonl').open() as f:
        for line in f:
            row=json.loads(line);key=row['interval_id']
            if key in intervals:raise ValueError('Repeated domain interval')
            intervals[key]=row;by_model[(row['model_id'],row['version'])].append(row)
    if len(intervals)!=receipt['unique_intervals']:raise ValueError('Interval universe differs')
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);counts=Counter();seen=set();source_seen=set();bindings={};written=0
    with (out/'inputs.jsonl').open('w') as manifest:
        for source in plan['coordinate_sources']:
            auditdir=Path(source['readback']);coordinates=Path(source['coordinates']);arpath=auditdir/'receipt.json';crpath=coordinates/'receipt.json'
            ar=json.loads(arpath.read_text());cr=json.loads(crpath.read_text());ap=json.loads(Path(source['readback_plan']).read_text())
            if ar['status']!='passed_duplication_exported_ca_readback' or ar['plan_sha256']!=sha(source['readback_plan']) or ar['producer_receipt_sha256']!=sha(crpath):raise ValueError('Incomplete coordinate readback')
            if cr['status']!='complete_duplication_coordinate_validation_with_dispositions' or cr['plan_sha256']!=sha(source['coordinate_plan']) or Path(ap['source']).resolve()!=coordinates.resolve() or Path(ap['source_plan']).resolve()!=Path(source['coordinate_plan']).resolve():raise ValueError('Coordinate source binding differs')
            bindings[str(arpath)]=sha(arpath);bindings[str(crpath)]=sha(crpath)
            for shard in cr['shards']:
                if shutil.disk_usage(out).free<plan['minimum_free_disk_gib']*2**30:raise RuntimeError('Disk reserve reached')
                path=coordinates/shard['output'];proofpath=auditdir/(path.name+'.readback.json');proof=json.loads(proofpath.read_text())
                if sha(proofpath)!=ar['proofs'][proofpath.name] or proof['source_sha256']!=shard['output_sha256'] or sha(path)!=shard['output_sha256']:raise ValueError('Unaudited coordinate shard')
                with gzip.open(path,'rt') as f:
                    for line in f:
                        record=json.loads(line);key=record['model_id'],record['version']
                        if key not in by_model:continue
                        if key in source_seen:raise ValueError('Repeated coordinate model')
                        source_seen.add(key)
                        for interval in by_model[key]:
                            if any(record[k]!=interval[k] for k in ['source_path','source_sha256','sequence_sha256']):raise ValueError('Domain/source provenance differs')
                            iid=interval['interval_id'];seen.add(iid)
                            for mask,threshold in [('full',None),('plddt70',70)]:
                                row=dict(interval_id=iid,model_id=key[0],version=key[1],mask=mask,start=interval['start'],end=interval['end'],interval_length=interval['length'],original_length=interval['original_length'],source_sha256=record['source_sha256'],coordinate_shard=str(path))
                                if record['status']=='rejected_content':row.update(status='source_rejected',reason=record['reason'])
                                else:
                                    if record['length']!=interval['original_length']:raise ValueError('Source length differs')
                                    blob,sequence,positions=render_interval(record,interval['start'],interval['end'],threshold)
                                    row.update(sequence=sequence,original_positions=positions,retained_residues=len(positions))
                                    if len(positions)<3:row['status']='too_few_retained_residues'
                                    else:
                                        dest=out/'pdb'/iid[:2]/(iid+'-'+mask+'.pdb');dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(blob);h=hashlib.sha256(blob).hexdigest()
                                        if sha(dest)!=h:raise ValueError('Written domain PDB differs')
                                        row.update(status='ready',path=str(dest),sha256=h);written+=len(blob)
                                counts[mask+':'+row['status']]+=1;manifest.write(json.dumps(row,separators=(',',':'))+'\n')
                print('Processed source shard',path.name,'intervals',len(seen),'/',len(intervals),flush=True)
    if seen!=set(intervals) or source_seen!=set(by_model):raise ValueError('Incomplete domain coordinate universe')
    verify()
    for path,h in bindings.items():
        if sha(path)!=h:raise ValueError('Coordinate receipts changed')
    result=dict(status='complete_duplication_domain_input_materialization_pending_readback',plan_sha256=ph,inventory_receipt_sha256=sha(inventory/'receipt.json'),coordinate_bindings=bindings,intervals=len(seen),source_models=len(source_seen),input_dispositions=sum(counts.values()),counts=dict(counts),pdb_bytes=written,artifacts={'inputs.jsonl':sha(out/'inputs.jsonl')},scope='All verified alignment/envelope domain intervals extracted from independently audited full-model coordinates. Full and pLDDT70 masks preserve original protein residue positions; short masks and rejected source models explicit. PDB rounding 0.001 A and confidence 0.01. Independent domain serialization readback and alignments pending; no PAE or biological inference.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
