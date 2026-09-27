#!/usr/bin/env python3
"""Independently verify all domain PDBs and dispositions against audited full-model arrays."""
import argparse,gzip,json,math,time
from collections import Counter,defaultdict
from pathlib import Path
import psutil
from Bio.Data.PDBData import protein_letters_3to1
from run_ortholog_pair_guide_comparison import sha


def check_input(row,interval,record):
    for field in ['interval_id','model_id','version','start','end','original_length','source_sha256']:
        if row[field]!=interval[field]:raise ValueError('Domain input identity differs')
    if row['interval_length']!=interval['length']:raise ValueError('Interval length differs')
    if record['status']=='rejected_content':
        if row['status']!='source_rejected' or row['reason']!=record['reason']:raise ValueError('Source rejection differs')
        return 0
    if record['status']!='validated' or row['mask'] not in ['full','plddt70']:raise ValueError('Unknown source/mask')
    positions=[p for p in range(interval['start'],interval['end']+1) if row['mask']=='full' or record['ca_plddt'][p-1]>=70]
    sequence=''.join(record['sequence'][p-1] for p in positions)
    if row['original_positions']!=positions or row['sequence']!=sequence or row['retained_residues']!=len(positions):raise ValueError('Domain residue selection differs')
    if len(positions)<3:
        if row['status']!='too_few_retained_residues' or 'path' in row:raise ValueError('Short-mask disposition differs')
        return 0
    if row['status']!='ready':raise ValueError('Missing ready domain input')
    path=Path(row['path']);blob=path.read_bytes()
    if sha(path)!=row['sha256']:raise ValueError('Changed domain PDB')
    lines=blob.decode().splitlines()
    if lines[-2:]!=['TER','END']:raise ValueError('Invalid PDB terminator')
    atoms=lines[:-2]
    if len(atoms)!=len(positions):raise ValueError('Domain atom count differs')
    for serial,(line,p) in enumerate(zip(atoms,positions),1):
        if len(line)!=80 or not line.startswith('ATOM  ') or line[12:16].strip()!='CA' or line[21]!='A' or line[16]!=' ' or line[26]!=' ':raise ValueError('Invalid domain atom/chain')
        if int(line[6:11])!=serial or int(line[22:26])!=p or protein_letters_3to1[line[17:20]]!=record['sequence'][p-1]:raise ValueError('Domain PDB residue mapping differs')
        for (a,b),value in zip([(30,38),(38,46),(46,54)],record['ca_xyz'][p-1]):
            actual=float(line[a:b])
            if not math.isfinite(actual) or abs(actual-value)>.000501:raise ValueError('Domain coordinate differs beyond rounding')
        confidence=float(line[60:66])
        if not math.isfinite(confidence) or abs(confidence-record['ca_plddt'][p-1])>.005001 or float(line[54:60])!=1.:raise ValueError('Domain confidence/occupancy differs')
    return len(positions)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        if sha(a.plan)!=ph:raise ValueError('Changed readback plan')
        for path,h in plan['pins'].items():
            if sha(path)!=h:raise ValueError('Changed pin: '+path)
    verify();dep=plan['producer']
    while psutil.pid_exists(dep['pid']):
        try:
            p=psutil.Process(dep['pid'])
            if abs(p.create_time()-dep['created'])>.01 or p.status()==psutil.STATUS_ZOMBIE:break
            if p.cmdline()!=dep['cmdline']:raise ValueError('Producer identity changed')
        except psutil.NoSuchProcess:break
        time.sleep(30)
    verify();sp=json.loads(Path(plan['source_plan']).read_text());folder=Path(sp['output']);rp=folder/'receipt.json';r=json.loads(rp.read_text());rh=sha(rp)
    if r['status']!='complete_background_domain_input_materialization_pending_readback' or r['plan_sha256']!=sha(plan['source_plan']):raise ValueError('Incomplete extraction')
    inventory=Path(sp['inventory']);ir=json.loads((inventory/'receipt.json').read_text())
    if r['inventory_receipt_sha256']!=sha(inventory/'receipt.json') or sha(inventory/'intervals.jsonl')!=ir['artifacts']['intervals.jsonl']:raise ValueError('Interval source differs')
    intervals={};by_model=defaultdict(list)
    with (inventory/'intervals.jsonl').open() as f:
        for line in f:
            row=json.loads(line);iid=row['interval_id']
            if iid in intervals:raise ValueError('Repeated interval')
            intervals[iid]=row;by_model[(row['model_id'],row['version'])].append(iid)
    manifest=folder/'inputs.jsonl';mh=sha(manifest)
    if mh!=r['artifacts']['inputs.jsonl']:raise ValueError('Changed input manifest')
    rows={}
    with manifest.open() as f:
        for line in f:
            row=json.loads(line);key=row['interval_id'],row['mask']
            if key in rows:raise ValueError('Duplicate domain input')
            rows[key]=row
    if set(rows)!={(iid,mask) for iid in intervals for mask in ['full','plddt70']}:raise ValueError('Incomplete domain input grid')
    seen=set();counts=Counter();atoms=0;size=0
    for source in sp['coordinate_sources']:
        coords=Path(source['coordinates']);audit=Path(source['readback']);crp=coords/'receipt.json';arp=audit/'receipt.json';cr=json.loads(crp.read_text());ar=json.loads(arp.read_text())
        if r['coordinate_bindings'][str(crp)]!=sha(crp) or r['coordinate_bindings'][str(arp)]!=sha(arp) or ar['producer_receipt_sha256']!=sha(crp):raise ValueError('Coordinate receipts differ')
        for shard in cr['shards']:
            path=coords/shard['output'];proof=audit/(path.name+'.readback.json')
            if sha(path)!=shard['output_sha256'] or sha(proof)!=ar['proofs'][proof.name] or json.loads(proof.read_text())['source_sha256']!=shard['output_sha256']:raise ValueError('Coordinate shard/proof differs')
            with gzip.open(path,'rt') as f:
                for line in f:
                    record=json.loads(line);key=record['model_id'],record['version']
                    if key not in by_model:continue
                    if key in seen:raise ValueError('Repeated coordinate model')
                    seen.add(key)
                    for iid in by_model[key]:
                        interval=intervals[iid]
                        if any(record[k]!=interval[k] for k in ['source_path','source_sha256','sequence_sha256']):raise ValueError('Source model provenance differs')
                        if record['status']=='validated' and record['length']!=interval['original_length']:raise ValueError('Source model length differs')
                        for mask in ['full','plddt70']:
                            row=rows[(iid,mask)]
                            if row['coordinate_shard']!=str(path):raise ValueError('Input coordinate shard differs')
                            atoms+=check_input(row,interval,record);counts[mask+':'+row['status']]+=1
                            if row['status']=='ready':size+=Path(row['path']).stat().st_size
            print('Checked domain source shard',path.name,'models',len(seen),'/',len(by_model),flush=True)
    if seen!=set(by_model) or dict(counts)!=r['counts'] or size!=r['pdb_bytes'] or len(rows)!=r['input_dispositions'] or len(intervals)!=r['intervals']:raise ValueError('Incomplete domain extraction readback')
    verify()
    for path,h in r['coordinate_bindings'].items():
        if sha(path)!=h:raise ValueError('Coordinate receipts changed during readback')
    if len(seen)!=r['source_models']:raise ValueError('Source model count differs')
    if sha(rp)!=rh or sha(manifest)!=mh:raise ValueError('Producer changed during readback')
    result=dict(status='passed_full_background_domain_input_readback',plan_sha256=ph,producer_receipt_sha256=rh,intervals=len(intervals),source_models=len(seen),input_dispositions=len(rows),counts=dict(counts),verified_pdb_ca_atoms=atoms,pdb_bytes=size,scope='All domain/mask dispositions independently reconstructed from audited full-model arrays. Every ready PDB atom identity, original position, XYZ, confidence and occupancy checked to fixed rounding tolerances. Source rejections and short masks preserved. Does not repeat raw-CIF audit, validate physical domain boundaries or establish PAE reliability or biological effects.')
    with Path(plan['output']).open('x') as f:json.dump(result,f,indent=2);f.write('\n')


if __name__=='__main__':main()
