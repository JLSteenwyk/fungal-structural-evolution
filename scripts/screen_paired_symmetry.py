#!/usr/bin/env python3
"""Run unfiltered symmetry diagnostics on every ready paired AA/3Di alignment."""
import argparse,csv,fcntl,hashlib,io,json,subprocess,time
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha,read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['inputs','fit_config','output']:p.add_argument('--'+name.replace('_','-'),type=Path,required=True)
    a=p.parse_args();inputs=checked_receipt(a.inputs);original=json.loads(a.fit_config.read_text())
    if original['input_receipt_sha256']!=sha(a.inputs/'receipt.json') or sha(Path(original['executable']))!=original['executable_sha256']:raise ValueError('Source/executable mismatch')
    ready=[r for r in read_table(a.inputs/'marker_summary.tsv') if r['status']=='ready_for_inference']
    if len(ready)!=inputs['ready_markers']:raise ValueError('Marker count differs')
    a.output.mkdir(parents=True,exist_ok=True);lock=(a.output/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    config={'input_receipt_sha256':sha(a.inputs/'receipt.json'),'source_fit_config_sha256':sha(a.fit_config),'script_sha256':sha(Path(__file__)),'executable':original['executable'],'executable_sha256':original['executable_sha256'],'version':original['version'],'workers':4,'threads_per_run':1,'memory_per_run_gb':2,'markers':len(ready),'alignments':2*len(ready),'nominal_pair_threshold':0.05,'interpretation':'Unfiltered IQ-TREE symmetry diagnostics on original paired masks; both alphabets treated as 20-state characters. No marker selection, hypothesis conclusion or multiplicity calibration. Nominal p-values and counts can be affected by sparse contingency tables, missingness and correlated/contextual 3Di characters. Non-rejection does not establish model adequacy.'}
    cp=a.output/'config.json'
    if cp.exists() and json.loads(cp.read_text())!=config:raise ValueError('Configuration changed')
    cp.write_text(json.dumps(config,indent=2)+'\n');ch=sha(cp)
    def run(job):
        m,label=job;marker=m['marker'];out=a.output/marker;out.mkdir(exist_ok=True);rp=out/(label+'.receipt.json');alignment=(a.inputs/marker/(label+'.faa')).resolve();prefix=(out/label).resolve()
        if rp.exists():
            r=json.loads(rp.read_text())
            if r['config_sha256']!=ch or r['alignment_sha256']!=sha(alignment) or any(sha(out/f)!=h for f,h in r['artifacts'].items()):raise ValueError('Completed diagnostic changed')
            return r
        seed=int(hashlib.sha256(marker.encode()).hexdigest()[:8],16)%2147483647
        command=[original['executable'],'-s',str(alignment),'-st','AA','-T','1','--mem','2G','--seed',str(seed),'-keep-ident','--prefix',str(prefix),'--symtest-only','--symtest-pval','0.05']
        start=time.monotonic()
        with (out/(label+'.stdout.log')).open('a') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
        path=out/(label+'.symtest.csv');content='\n'.join(x for x in path.read_text().splitlines() if x.strip() and not x.startswith('#'));rows=list(csv.DictReader(io.StringIO(content)))
        if len(rows)!=1:raise ValueError('Expected one unpartitioned diagnostic row')
        r={'status':'complete_unfiltered_symmetry_diagnostic','marker':marker,'alphabet':label,'config_sha256':ch,'alignment_sha256':sha(alignment),'command':command,'taxa':int(m['eligible_taxa']),'columns':int(m['retained_columns']),'elapsed_seconds':time.monotonic()-start,'raw_result':rows[0],'artifacts':{f:sha(out/f) for f in [label+'.symtest.csv',label+'.log',label+'.stdout.log']}}
        rp.write_text(json.dumps(r,indent=2)+'\n');print(marker,label,flush=True);return r
    results=[]
    with ThreadPoolExecutor(max_workers=4) as pool:
        for f in as_completed([pool.submit(run,(m,label)) for m in ready for label in ['aa','3di']]):results.append(f.result())
    rows=[{'marker':r['marker'],'alphabet':r['alphabet'],'taxa':r['taxa'],'columns':r['columns'],**r['raw_result']} for r in sorted(results,key=lambda r:(r['marker'],r['alphabet']))]
    write_table(a.output/'raw_symmetry_summary.tsv',rows)
    r={'status':'complete_unfiltered_paired_symmetry_screen','config_sha256':ch,'markers':len(ready),'alignments':len(results),'entry_receipts':{x['marker']+'/'+x['alphabet']:sha(a.output/x['marker']/(x['alphabet']+'.receipt.json')) for x in results},'interpretation':config['interpretation'],'artifacts':{'raw_symmetry_summary.tsv':sha(a.output/'raw_symmetry_summary.tsv')}}
    (a.output/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k!='entry_receipts'},indent=2))

if __name__=='__main__':main()
