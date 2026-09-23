#!/usr/bin/env python3
"""Rebuild cluster composition from source rows with Python sets and counters."""
import argparse,csv,json,sqlite3,time
from pathlib import Path
from collections import Counter,defaultdict
import psutil
from annotate_whole_proteome_cluster_composition import sha

def connect(path):return sqlite3.connect('file:'+str(Path(path).resolve())+'?mode=ro',uri=True)
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        if sha(args.plan)!=ph:raise ValueError('Audit plan changed')
        for f,h in plan['pins'].items():
            if sha(f)!=h:raise ValueError('Changed pinned input: '+f)
    verify();pred=plan['predecessor']
    while True:
        try:
            proc=psutil.Process(pred['pid'])
            live=proc.create_time()==pred['create_time'] and proc.status()!=psutil.STATUS_ZOMBIE
        except psutil.NoSuchProcess:live=False
        if not live:break
        time.sleep(20)
    verify();p=json.loads(Path(plan['producer_plan']).read_text());root=Path(p['output'])
    r=json.loads((root/'receipt.json').read_text())
    if r['status']!='complete_cluster_composition_pending_independent_readback' or r['plan_sha256']!=sha(plan['producer_plan']):raise ValueError('Producer incomplete or mismatched')
    for f,h in p['pins'].items():
        if sha(f)!=h:raise ValueError('Changed producer input: '+f)
    for f,h in r['artifacts'].items():
        if sha(root/f)!=h:raise ValueError('Changed output: '+f)
    out=Path(plan['output'])
    if out.exists():raise FileExistsError(out)
    out.mkdir();db=connect(root/'cluster_composition.sqlite');members={}
    for rep,model in csv.reader(Path(p['members']).open(),delimiter='\t'):
        if model in members:raise ValueError('Duplicate source model')
        members[model]=rep
    remaining=dict(members)
    for model,rep in db.execute('SELECT model,representative FROM members'):
        if remaining.pop(model,None)!=rep:raise ValueError('Membership differs')
    if remaining:raise ValueError('Missing model rows')
    model_counts=Counter(members.values());genes={};taxa=defaultdict(set);protein_counts=Counter();covered=set()
    source=connect(p['structures'])
    output_rows=iter(db.execute('SELECT native_gene_id,taxon_id,protein_id,model,representative FROM links ORDER BY native_gene_id'))
    for gene,taxon,protein,mid,version,path in source.execute('SELECT native_gene_id,taxon_id,protein_id,model_id,version,model_path FROM structures ORDER BY native_gene_id'):
        model=Path(path).stem
        if model!=f'{mid}-v{version}' or model not in members:raise ValueError('Source model identity differs')
        rep=members[model];expected=(gene,taxon,protein,model,rep)
        if next(output_rows,None)!=expected or gene in genes:raise ValueError('Protein linkage differs')
        genes[gene]=rep;covered.add(model);taxa[rep].add(taxon);protein_counts[rep]+=1
    if next(output_rows,None) is not None or covered!=set(members):raise ValueError('Unexpected links or missing models')
    family=connect(p['families']);family_counts={}
    for guide in ['profile','mafft']:
        output_rows=iter(db.execute('SELECT native_gene_id,family FROM family_links WHERE guide=? ORDER BY native_gene_id',(guide,)))
        groups=defaultdict(set);matched=0
        for gene,fam in family.execute('SELECT native_gene_id,family FROM assignments WHERE guide=? ORDER BY native_gene_id',(guide,)):
            if gene not in genes:continue
            if next(output_rows,None)!=(gene,fam):raise ValueError('Family mapping differs')
            groups[genes[gene]].add(fam);matched+=1
        if next(output_rows,None) is not None or matched!=len(genes):raise ValueError('Family coverage differs')
        family_counts[guide]={k:len(v) for k,v in groups.items()}
    seen=set()
    for row in csv.DictReader((root/'cluster_composition.tsv').open(),delimiter='\t'):
        rep=row['representative']
        if rep in seen or rep not in model_counts:raise ValueError('Unexpected summary cluster')
        expected=[model_counts[rep],protein_counts[rep],len(taxa[rep]),family_counts['profile'][rep],family_counts['mafft'][rep]]
        if [int(row[k]) for k in ['models','proteins','taxa','profile_families','mafft_families']]!=expected:raise ValueError('Composition counts differ')
        seen.add(rep)
    if seen!=set(model_counts) or len(genes)!=r['protein_links'] or len(members)!=r['models'] or len(seen)!=r['clusters']:raise ValueError('Incomplete universe')
    verify()
    result=dict(status='passed_full_cluster_composition_source_readback',plan_sha256=ph,producer_receipt_sha256=sha(root/'receipt.json'),models=len(members),protein_links=len(genes),clusters=len(seen),family_assignments=2*len(genes),scope='Every cluster membership and source protein identity, both complete family-assignment joins and every composition-table count rebuilt with Python sets/counters. Does not validate structural alignment thresholds, homology, orthology, confidence or evolutionary events.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
