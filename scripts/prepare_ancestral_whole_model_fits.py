#!/usr/bin/env python3
"""Freeze the complete whole-protein evolutionary-model fit grid."""
import json
from pathlib import Path
from Bio import SeqIO
from prepare_case_ancestral_neighborhoods import sha,read

def main():
    root=Path('results/ancestral/case-alignments-20260927-v1');audit=Path('results/ancestral/case-alignment-readback-20260927-v1');trees=Path('results/ancestral/case-local-trees-20260927-v1')
    receipt=json.loads((root/'receipt.json').read_text());ar=json.loads((audit/'receipt.json').read_text());assert ar['source_receipt_sha256']==sha(root/'receipt.json')
    closure=Path('metadata/ancestral_case_alignments_completed_20260927.json');assert json.loads(closure.read_text())['source_receipt_sha256']==sha(audit/'receipt.json')
    tr=json.loads((trees/'receipt.json').read_text());tree_rows=read(trees/'tree_inputs.tsv');pins={str(p):sha(p) for p in [root/'receipt.json',audit/'receipt.json',trees/'receipt.json',closure]};jobs=[]
    for r in receipt['results']:
        family,method=r['family'],r['method'];alignment=root/(family+'-'+method)/'alignment.faa';assert sha(alignment)==r['artifacts']['alignment.faa']
        tree=trees/('profile-'+family+'-whole.nwk');alias=trees/('mafft-'+family+'-whole.nwk')
        assert sha(tree)==tr['artifacts'][tree.name] and sha(alias)==tr['artifacts'][alias.name] and sha(tree)==sha(alias)
        pins[str(tree)]=sha(tree);pins[str(alias)]=sha(alias);pins[str(alignment)]=sha(alignment)
        records=list(SeqIO.parse(alignment,'fasta'));chars=''.join(str(x.seq) for x in records);assert set(chars)<=set('ARNDCQEGHILKMFPSTWYV-X')
        for model in ['LG','WAG','JTT']:
            jobs.append(dict(job_id=family+'-whole-'+method+'-'+model,family=family,boundary='whole',method=method,model=model+'+F+G4',alignment=str(alignment),alignment_sha256=sha(alignment),tree=str(tree),tree_sha256=sha(tree),proteins=r['proteins'],columns=r['columns'],unknown_X_residues=chars.count('X'),guide_aliases=['profile','mafft']))
    assert len(jobs)==78 and sum(j['unknown_X_residues'] for j in jobs)==18
    old=json.loads(Path('metadata/ancestral_domain_model_fit_plan_20260927.json').read_text())
    for f in ['scripts/prepare_ancestral_whole_model_fits.py','scripts/run_ancestral_whole_model_fits.py','scripts/prepare_case_ancestral_neighborhoods.py',old['iqtree']]:pins[f]=sha(f)
    p=dict(iqtree=old['iqtree'],output='results/ancestral/whole-protein-model-fits-20260927-v1',jobs=jobs,pins=pins,resources=dict(concurrent_fits=2,threads_per_fit=4,aggregate_cpus=8,memory_gib=12,swap_gib=0,output_gib=4,planning_hours=[1,48]),scope='26 complete untrimmed whole-protein alignments x3 models, all1025 input proteins. X and gaps retained as unknown observations. Shared guide aliases fitted once. Full model/likelihood audit and ancestral inference remain downstream.')
    Path('metadata/ancestral_whole_model_fit_plan_20260927.json').write_text(json.dumps(p,indent=2)+'\n');print('Prepared78 whole-protein fits')

if __name__=='__main__':main()
