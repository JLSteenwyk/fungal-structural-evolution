#!/usr/bin/env python3
"""Check complete retrieval dispositions and every native PAE matrix entry."""
import gzip,hashlib,json,math,subprocess,time
from pathlib import Path
import psutil
from screen_duplication_domain_alignment_coverage import sha


def main():
    launchpath=Path('metadata/whole_domain_case_pae_launch_20260927.json');launch=json.loads(launchpath.read_text());lh=sha(launchpath)
    while True:
        try:
            p=psutil.Process(launch['pid'])
            if p.create_time()!=launch['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:break
        print('Waiting for exact case PAE retrieval',launch['pid'],flush=True);time.sleep(30)
    state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
    assert sha(launchpath)==lh and sha('scripts/retrieve_marker_pae.py')==launch['script_sha256']
    base=Path('results/structural_comparisons');snapshot=base/'whole-domain-case-pae-inputs-20260927-v1';root=base/'whole-domain-case-pae-20260927-v1'
    sr=snapshot/'receipt.json';assert sha(sr)==launch['snapshot_receipt_sha256'];s=json.loads(sr.read_text())
    mp=snapshot/'model_provenance.json';assert sha(mp)==s['artifacts'][mp.name]
    models=json.loads(mp.read_text());expected={(m['model_id'],m['version']):m for m in models};assert len(expected)==39
    rp=root/'receipt.json';r=json.loads(rp.read_text());assert r['mapping_receipt_sha256']==sha(sr)
    manifest=root/'pae_manifest.json';assert sha(manifest)==r['artifacts'][manifest.name]
    records=json.loads(manifest.read_text());assert len(records)==39 and {(m['model_id'],m['version']) for m in records}==set(expected)
    cells=0;verified=0;failures=[];hashes={str(p):sha(p) for p in [launchpath,sr,mp,rp,manifest]}
    for record in records:
        m=expected[record['model_id'],record['version']]
        if record['status']!='verified':
            assert record['status']=='failed' and record['error'];failures.append(record);continue
        verified+=1
        for field in ['model_id','version','length','sequence_sha256']:assert record[field]==m[field]
        assert record['url']==m['pae_url'];path=Path(record['path']);assert sha(path)==record['gzip_sha256']
        raw=gzip.decompress(path.read_bytes());assert hashlib.sha256(raw).hexdigest()==record['json_sha256']
        assert len(raw)==record['json_bytes'] and path.stat().st_size==record['compressed_bytes']
        payload=json.loads(raw);assert isinstance(payload,list) and len(payload)==1
        matrix=payload[0]['predicted_aligned_error'];maximum=float(payload[0]['max_predicted_aligned_error'])
        assert math.isfinite(maximum) and maximum>=0 and len(matrix)==m['length']
        for row in matrix:
            assert len(row)==m['length']
            for value in row:
                assert not isinstance(value,bool) and isinstance(value,(int,float)) and math.isfinite(value) and 0<=value<=maximum+.51
                cells+=1
        hashes[str(path)]=sha(path)
    assert verified==r['models_verified'] and len(failures)==r['models_failed'] and verified+len(failures)==r['models_requested']==39
    assert cells==sum(m['length']**2 for m in models if any(x['model_id']==m['model_id'] and x['version']==m['version'] and x['status']=='verified' for x in records))
    for path,digest in hashes.items():assert sha(path)==digest
    result=dict(status='passed_full_case_pae_retrieval_readback',terminal_state=state,source_hashes=hashes,checker_sha256=sha(__file__),models_requested=39,models_verified=verified,models_failed=len(failures),matrix_entries_checked=cells,failures=failures,scope='Every disposition and matrix entry checked; exact version, sequence-provenance, URL, dimensions, compression and hashes retained. PAE is prediction confidence, not independent experimental validation or calibrated uncertainty for structural contrasts.')
    Path('metadata/whole_domain_case_pae_readback_20260927.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
