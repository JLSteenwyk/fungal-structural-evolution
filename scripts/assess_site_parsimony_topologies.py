#!/usr/bin/env python3
"""Evaluate original paired characters on every saved AA bootstrap topology."""
import argparse,fcntl,json
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import numpy as np
from Bio import Phylo,SeqIO
from audit_busco_gene_copies import sha,read_table
from assess_pae_sensitivity import checked_receipt
from prepare_paired_phylogenetic_inputs import write_table
from summarize_site_parsimony_exposure import minimum_changes


def run_marker(job):
    marker,inputs,fits,output,config_hash,reference,expected=job
    out=Path(output)/marker;out.mkdir(parents=True,exist_ok=True);rp=out/'receipt.json'
    if rp.exists():
        r=checked_receipt(out)
        if r['config_sha256']!=config_hash:raise ValueError('Marker configuration changed')
        return r
    folder=Path(inputs)/marker;fit=Path(fits)/marker;boot=fit/'aa.ufboot'
    r=json.loads((fit/'aa.receipt.json').read_text())
    if sha(boot)!=r['artifacts'][boot.name]:raise ValueError('Bootstrap tree file changed')
    seqs=[{x.id:str(x.seq) for x in SeqIO.parse(folder/file,'fasta')} for file in ['aa.faa','3di.faa']]
    if len(seqs[0])>65000:raise ValueError('Count storage range too small')
    aa=[];di=[]
    for tree in Phylo.parse(boot,'newick'):
        aa.append(minimum_changes(tree,seqs[0]));di.append(minimum_changes(tree,seqs[1]))
    if len(aa)!=expected:raise ValueError('Incomplete bootstrap tree count')
    arrays={'aa':np.asarray(aa,dtype=np.uint16),'3di':np.asarray(di,dtype=np.uint16)}
    if any(v.shape!=(expected,len(reference)) for v in arrays.values()):raise ValueError('Bootstrap site grid differs')
    rows=[dict(x) for x in reference]
    for label,values in arrays.items():
        quantiles=np.quantile(values,[.025,.5,.975],axis=0)
        for j,row in enumerate(rows):
            original=int(row[label+'_minimum_changes']);v=values[:,j]
            row.update({label+'_topology_min':int(v.min()),label+'_topology_p025':float(quantiles[0,j]),label+'_topology_median':float(quantiles[1,j]),label+'_topology_p975':float(quantiles[2,j]),label+'_topology_max':int(v.max()),label+'_topology_fraction_equal_reference':float(np.mean(v==original)),label+'_topology_distinct_scores':len(np.unique(v))})
    np.savez_compressed(out/'bootstrap_site_scores.npz',**arrays)
    write_table(out/'site_topology_sensitivity.tsv',rows)
    # Read stored arrays back before treating this marker as complete.
    with np.load(out/'bootstrap_site_scores.npz') as saved:
        if any(not np.array_equal(saved[k],v) for k,v in arrays.items()):raise ValueError('Stored scores differ')
    result={'status':'complete_marker_topology_sensitivity','marker':marker,'config_sha256':config_hash,'bootstrap_tree_sha256':sha(boot),'bootstrap_trees':expected,'sites':len(reference),'aa_sites_with_variable_scores':int(np.sum(np.ptp(arrays['aa'],axis=0)>0)),'3di_sites_with_variable_scores':int(np.sum(np.ptp(arrays['3di'],axis=0)>0)),'artifacts':{f.name:sha(f) for f in out.iterdir() if f.name!='receipt.json'}}
    rp.write_text(json.dumps(result,indent=2)+'\n');return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['inputs','fits','diagnostic','output']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--workers',type=int,default=8);p.add_argument('--bootstrap-trees',type=int,default=1000);a=p.parse_args()
    if a.workers<1 or a.bootstrap_trees<1:raise ValueError('Invalid resources')
    d=checked_receipt(a.diagnostic);checked_receipt(a.inputs)
    if d['source_receipts']['inputs']!=sha(a.inputs/'receipt.json') or d['source_receipts']['fits']!=sha(a.fits/'receipt.json'):raise ValueError('Source lineage differs')
    grouped={}
    for r in read_table(a.diagnostic/'site_parsimony_exposure.tsv'):grouped.setdefault(r['marker'],[]).append(r)
    if len(grouped)!=d['markers']:raise ValueError('Marker count differs')
    for rows in grouped.values():
        if [int(r['paired_column_1based']) for r in rows]!=list(range(1,len(rows)+1)):raise ValueError('Reference column grid differs')
    a.output.mkdir(parents=True,exist_ok=True);lock=(a.output/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    config={'source_receipts':{k:sha(getattr(a,k)/'receipt.json') for k in ['inputs','fits','diagnostic']},'script_sha256':sha(Path(__file__)),'parsimony_helper_sha256':sha(Path(__file__).with_name('summarize_site_parsimony_exposure.py')),'workers':a.workers,'bootstrap_trees_per_marker':a.bootstrap_trees,'markers':len(grouped),'interpretation':'All saved AA UFBoot topologies evaluated using the unchanged original paired AA/3Di characters; isolates topology sensitivity conditional on the AA inference procedure. Quantiles are not calibrated confidence intervals or posterior probabilities. Bootstrap trees may repeat and are all retained. No new topology fits, sites, branch assignments or ancestral reconstructions are generated.'}
    cp=a.output/'config.json'
    if cp.exists() and json.loads(cp.read_text())!=config:raise ValueError('Configuration changed')
    cp.write_text(json.dumps(config,indent=2)+'\n');results=[]
    jobs=[(m,str(a.inputs),str(a.fits),str(a.output),sha(cp),rows,a.bootstrap_trees) for m,rows in sorted(grouped.items())]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futures=[pool.submit(run_marker,j) for j in jobs]
        for f in as_completed(futures):
            r=f.result();results.append(r);print(r['marker'],r['sites'],r['aa_sites_with_variable_scores'],r['3di_sites_with_variable_scores'],flush=True)
    allrows=[]
    for r in sorted(results,key=lambda r:r['marker']):allrows.extend(read_table(a.output/r['marker']/'site_topology_sensitivity.tsv'))
    if len(allrows)!=d['sites']:raise ValueError('Full site count differs')
    write_table(a.output/'site_topology_sensitivity.tsv',allrows)
    r={'status':'complete_all_marker_topology_sensitivity','config_sha256':sha(cp),'markers':len(results),'sites':len(allrows),'bootstrap_topologies_evaluated':sum(x['bootstrap_trees'] for x in results),'aa_sites_with_variable_scores':sum(x['aa_sites_with_variable_scores'] for x in results),'3di_sites_with_variable_scores':sum(x['3di_sites_with_variable_scores'] for x in results),'marker_receipts':{x['marker']:sha(a.output/x['marker']/'receipt.json') for x in results},'interpretation':config['interpretation'],'artifacts':{'site_topology_sensitivity.tsv':sha(a.output/'site_topology_sensitivity.tsv')}}
    (a.output/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k!='marker_receipts'},indent=2))

if __name__=='__main__':main()
