#!/usr/bin/env python3
"""Check every background membership by independent fixed-record binary search."""
import argparse,csv,json,time
from pathlib import Path
import psutil
from run_ortholog_pair_guide_comparison import sha


def record_key(handle,index):
    handle.seek(index*28);raw=handle.read(28)
    assert len(raw)==28 and raw[12:14]==b'0\n' and raw[26:28]==b'1\n' and raw[:12]==raw[14:26]
    key=raw[:12];assert all(c in b'0123456789abcdef' for c in key) and key[:6]<key[6:]
    return key


def membership(handle,n,key):
    lo,hi=0,n
    while lo<hi:
        middle=(lo+hi)//2
        if record_key(handle,middle)<key:lo=middle+1
        else:hi=middle
    if lo<n:
        after=record_key(handle,lo);assert after>=key
        if after==key:return 1
    if lo>0:assert record_key(handle,lo-1)<key
    return 0


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--config',type=Path,required=True);a=ap.parse_args()
    config=json.loads(a.config.read_text());ch=sha(a.config)
    def verify_config():
        assert sha(a.config)==ch
        for path,h in config['pins'].items():assert sha(path)==h,path
    verify_config()
    for dep in config['dependencies']:
        while True:
            try:
                p=psutil.Process(dep['pid'])
                if p.create_time()!=dep['created'] or p.status()==psutil.STATUS_ZOMBIE:break
                assert p.cmdline()==dep['cmdline']
            except psutil.NoSuchProcess:break
            time.sleep(30)
    verify_config();plan=json.loads(Path(config['plan']).read_text());root=Path(plan['output']);source=Path(plan['source'])
    receipt=json.loads((root/'receipt.json').read_text());sr=json.loads((source/'receipt.json').read_text());audit=json.loads(Path(config['guide_readback']).read_text())
    assert receipt['status']=='complete_background_native_ortholog_membership_pending_readback' and receipt['plan_sha256']==sha(config['plan'])
    assert audit['status']=='passed_full_terminal_sister_guide_comparison_readback'
    assert receipt['source_receipt_sha256']==audit['producer_receipt_sha256']==sha(source/'receipt.json')
    bindings={str(root/'receipt.json'):sha(root/'receipt.json'),str(source/'receipt.json'):sha(source/'receipt.json'),config['guide_readback']:sha(config['guide_readback'])}
    bindings.update(plan['pins'])
    bindings.update({str(root/name):h for name,h in receipt['artifacts'].items()})
    source_table=source/'modeled_candidate_union.tsv';bindings[str(source_table)]=sr['artifacts'][source_table.name]
    for stream in plan['streams']:
        r=json.loads(Path(stream['receipt']).read_text());bindings[stream['stream']]=r['sorted_pairs_sha256']
    def verify():
        verify_config()
        for path,h in bindings.items():assert sha(path)==h,path
    verify()
    with source_table.open() as f:source_rows={(r['gene_a'],r['gene_b']):r for r in csv.DictReader(f,delimiter='\t')}
    assert len(source_rows)==receipt['candidate_rows']==audit['modeled_candidate_union']
    species={}
    for line in Path(plan['species_ids']).read_text().splitlines():
        number,name=line.split(': ',1);species[number]=name.rsplit('.',1)[0]
    needed={g for pair in source_rows for g in pair};ids={}
    with Path(plan['sequence_ids']).open() as f:
        for index,line in enumerate(f):
            native,protein=line.rstrip('\n').split(': ',1);gene=species[native.partition('_')[0]]+'_'+protein
            if gene in needed:
                assert gene not in ids;ids[gene]=index
    assert set(ids)==needed
    rows=[];seen=set();keys=[]
    with (root/'candidate_orthology_membership.tsv').open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            pair=(row['gene_a'],row['gene_b']);assert pair in source_rows and pair not in seen;seen.add(pair)
            copied={k:v for k,v in row.items() if k not in ['native_pair_key','profile_native_ortholog','mafft_native_ortholog']}
            assert copied==source_rows[pair]
            x,y=sorted(ids[g] for g in pair);key=f'{x:06x}{y:06x}'
            assert row['native_pair_key']==key and (not keys or keys[-1]<key)
            keys.append(key);rows.append(row)
    assert seen==set(source_rows)
    assert (root/'queries.hex').read_text()==''.join(k+'\n' for k in keys)
    summaries={}
    for stream in plan['streams']:
        guide=stream['guide'];path=Path(stream['stream']);r=json.loads(Path(stream['receipt']).read_text());n=r['unique_unordered_pairs']
        assert path.stat().st_size==n*28
        with (root/(guide+'_membership.tsv')).open() as f:raw=list(csv.reader(f,delimiter='\t'))
        assert len(raw)==len(rows)
        present=0
        with path.open('rb') as handle:
            for i,row in enumerate(rows):
                found=membership(handle,n,row['native_pair_key'].encode());present+=found
                assert row[guide+'_native_ortholog']==str(found)
                assert raw[i]==[row['native_pair_key'],str(found)]
        summary=dict(stream_pairs=n,queries=len(rows),present=present,absent=len(rows)-present)
        assert summary==receipt['guides'][guide];summaries[guide]=summary
        print(guide,json.dumps(summary),flush=True)
    verify()
    result=dict(status='passed_full_background_orthology_membership_readback',config_sha256=ch,producer_receipt_sha256=sha(root/'receipt.json'),guide_readback_sha256=sha(config['guide_readback']),candidate_rows=len(rows),guides=summaries,checker_sha256=sha(__file__),scope='Every candidate source field and protein ordinal reconstructed; all membership flags checked by independent binary search with reciprocal-record validation and bracketing for absences. Full stream hashes bind to prior global ordering/multiplicity audits. Guide union readback required. Native assignment is not independent biological orthology or a matched duplication effect.')
    output=Path(config['output']);assert not output.exists();output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
