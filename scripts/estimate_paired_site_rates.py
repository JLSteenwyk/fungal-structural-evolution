#!/usr/bin/env python3
"""Refit existing paired model specifications on fixed AA topologies and export site rates."""
import argparse,fcntl,json,math,subprocess,time
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha,read_table
from run_paired_marker_fits import tree_edges


def rate_rows(path,expected):
    lines=[x.split() for x in path.read_text().splitlines() if x.strip() and not x.startswith('#')]
    header=[{'Cat':'Category','C_Rate':'Categorized_rate'}.get(x,x) for x in lines[0]]
    if len(set(header))!=len(header):raise ValueError('Ambiguous site-rate columns')
    rows=[dict(zip(header,x)) for x in lines[1:]]
    if any(len(x)!=len(header) for x in lines[1:]) or not {'Site','Rate','Category','Categorized_rate'}<=set(header):raise ValueError('Unexpected site-rate schema')
    if len(rows)!=expected or [int(x['Site']) for x in rows]!=list(range(1,expected+1)):raise ValueError('Incomplete rate site grid')
    for r in rows:
        if not 1<=int(r['Category'])<=4 or any(not math.isfinite(float(r[k])) or float(r[k])<0 for k in ['Rate','Categorized_rate']):raise ValueError('Invalid site rate')
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['inputs','fits','audit','output']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--heterogeneity', choices=['G4','R4'], default='G4')
    a=p.parse_args();inputs=checked_receipt(a.inputs);audit=checked_receipt(a.audit)
    source=json.loads((a.fits/'receipt.json').read_text());parent=json.loads((a.fits/'config.json').read_text())
    if audit['fit_receipt_sha256']!=sha(a.fits/'receipt.json') or source['config_sha256']!=sha(a.fits/'config.json') or parent['input_receipt_sha256']!=sha(a.inputs/'receipt.json'):raise ValueError('Source lineage differs')
    if sha(Path(parent['executable']))!=parent['executable_sha256']:raise ValueError('IQ-TREE executable changed')
    ready=[r for r in read_table(a.inputs/'marker_summary.tsv') if r['status']=='ready_for_inference']
    if {r['marker'] for r in ready}!={r['marker'] for r in source['results']} or len(ready)!=inputs['ready_markers']:raise ValueError('Marker universe differs')
    labels=['aa','3di_af','3di_af_empirical','3di_llm'];pins={}
    for m in ready:
        marker=m['marker']
        for label in labels:
            cp=a.fits/marker/(label+'.config.json');fr=json.loads((a.fits/marker/(label+'.receipt.json')).read_text())
            if sha(cp)!=fr['config_sha256']:raise ValueError('Source fit configuration changed')
            pins[str(cp)]=sha(cp)
            command=json.loads(cp.read_text())['command'];model=command[command.index('-m')+1].split('+')[0]
            if model!='LG':pins[model]=sha(Path(model))
        tree=a.fits/marker/'aa.treefile';fr=json.loads((a.fits/marker/'aa.receipt.json').read_text())
        if sha(tree)!=fr['artifacts'][tree.name]:raise ValueError('Source topology changed')
        pins[str(tree)]=sha(tree)
    a.output.mkdir(parents=True,exist_ok=True);lock=(a.output/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    config={'source_receipts':{k:sha(getattr(a,k)/'receipt.json') for k in ['inputs','fits','audit']},'script_sha256':sha(Path(__file__)),'tree_helper_sha256':sha(Path(__file__).with_name('run_paired_marker_fits.py')),'executable':parent['executable'],'executable_sha256':parent['executable_sha256'],'pinned_files':pins,'workers':4,'threads_per_fit':1,'memory_per_fit_gb':2,'markers':len(ready),'fits':len(ready)*4,'interpretation':'Existing AA LG+F+G4 and three 3Di G4 specifications, refitted on fixed original AA topologies with original paired masks. Branch lengths and model parameters reestimated; not an exact replay of previous parameter values. Empirical-Bayes posterior mean site-rate multipliers are conditional on these models and trees, not rates/year, physical displacement, selection or calibrated uncertainty. Model adequacy and Gamma-versus-other heterogeneity sensitivity remain pending.'}
    config['heterogeneity']=a.heterogeneity
    if a.heterogeneity=='R4':
        config['interpretation']='Four-category FreeRate sensitivity for all four existing matrix/frequency specifications on the same paired observations and original fixed AA topology. Branches, category rates and weights reestimated. Comparison with Gamma4 is conditional model sensitivity, not model adequacy, calibrated uncertainty, physical displacement or selection.'
    cp=a.output/'config.json'
    if cp.exists() and json.loads(cp.read_text())!=config:raise ValueError('Run configuration changed')
    cp.write_text(json.dumps(config,indent=2)+'\n');config_hash=sha(cp)
    def run(m):
        marker=m['marker'];folder=a.output/marker;folder.mkdir(exist_ok=True);results=[]
        taxa={x.id for x in SeqIO.parse(a.inputs/marker/'aa.faa','fasta')};topology=(a.fits/marker/'aa.treefile').resolve()
        for label in labels:
            prefix=(folder/label).resolve();rp=folder/(label+'.receipt.json')
            if rp.exists():
                r=json.loads(rp.read_text())
                if r['parent_config_sha256']!=config_hash or any(sha(folder/f)!=h for f,h in r['artifacts'].items()):raise ValueError('Completed output changed')
                results.append(r);continue
            original=json.loads((a.fits/marker/(label+'.config.json')).read_text())['command'];command=[];i=0
            while i<len(original):
                value=original[i]
                if value in ['--alrt','-B','-te']:i+=2;continue
                if value in ['--bnni','--boot-trees']:i+=1;continue
                if value=='--prefix':command.extend([value,str(prefix)]);i+=2;continue
                command.append(value);i+=1
            command.extend(['-te',str(topology),'--rate','--sitelh'])
            mi=command.index('-m')+1
            if not command[mi].endswith('+G4'):raise ValueError('Expected Gamma4 source specification')
            command[mi]=command[mi][:-2]+a.heterogeneity
            align=Path(command[command.index('-s')+1]);expected_align=json.loads((a.fits/marker/(label+'.config.json')).read_text())['alignment_sha256']
            if sha(align)!=expected_align:raise ValueError('Alignment changed')
            request={'command':command,'alignment_sha256':sha(align),'topology_sha256':sha(topology),'parent_config_sha256':config_hash}
            request_path=folder/(label+'.config.json')
            if request_path.exists() and json.loads(request_path.read_text())!=request:raise ValueError('Fit request changed')
            request_path.write_text(json.dumps(request,indent=2)+'\n');start=time.monotonic()
            with (folder/(label+'.stdout.log')).open('a') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
            rows=rate_rows(folder/(label+'.rate'),int(m['retained_columns']))
            if set(tree_edges(folder/(label+'.treefile'),taxa))!=set(tree_edges(topology,taxa)):raise ValueError('Fixed topology changed')
            required=[label+x for x in ['.rate','.sitelh','.treefile','.iqtree','.log','.config.json']]
            r={'status':'complete_fixed_topology_site_rates','marker':marker,'label':label,'parent_config_sha256':config_hash,'sites':len(rows),'elapsed_seconds':time.monotonic()-start,'artifacts':{f:sha(folder/f) for f in required}}
            rp.write_text(json.dumps(r,indent=2)+'\n');results.append(r)
        print(marker,'four site-rate fits complete',flush=True);return results
    results=[]
    with ThreadPoolExecutor(max_workers=4) as pool:
        for f in as_completed([pool.submit(run,m) for m in ready]):results.extend(f.result())
    if any(sha(Path(p))!=h for p,h in pins.items()):raise ValueError('Pinned source changed during run')
    out={'status':'complete_paired_site_rate_exports','config_sha256':config_hash,'markers':len(ready),'fits':len(results),'site_rate_rows':sum(r['sites'] for r in results),'fit_receipts':{r['marker']+'/'+r['label']:sha(a.output/r['marker']/(r['label']+'.receipt.json')) for r in results},'interpretation':config['interpretation']}
    (a.output/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='fit_receipts'},indent=2))

if __name__=='__main__':main()
