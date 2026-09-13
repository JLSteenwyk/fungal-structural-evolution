#!/usr/bin/env python3
"""Realign all directed representative/member edges without a search prefilter."""
import argparse,csv,json,math,shutil,subprocess
from collections import Counter,defaultdict
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha,read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--clusters',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use new immutable edge validation output')
    r=checked_receipt(a.clusters);c=json.loads((a.clusters/'config.json').read_text())
    if r['config_sha256']!=sha(a.clusters/'config.json'):raise ValueError('Cluster config differs')
    db=(a.clusters/'tmp/latest/input').resolve();lookup={line.split('\t')[1]:int(line.split('\t')[0]) for line in (db.with_suffix('.lookup')).read_text().splitlines()}
    members=read_table(a.clusters/'model_cluster_membership.tsv');models={r['cluster_input_id']:r for r in members}
    if set(models)!=set(lookup):raise ValueError('Database model universe differs')
    expected=set();edges=[]
    for row in members:
        rep=row['representative_input_id'];member=row['cluster_input_id']
        if rep!=member:edges.append((rep,member));expected.update([(rep,member),(member,rep)])
    if len(expected)!=2*len(edges):raise ValueError('Duplicate representative/member edges')
    a.output.mkdir(parents=True);pairs=a.output/'requested_pairs.tsv'
    with pairs.open('w') as f:
        for q,t in sorted(expected,key=lambda pair:(lookup[pair[0]],lookup[pair[1]])):f.write(f'{lookup[q]}\t{lookup[t]}\t0\t0\n')
    exe=Path(shutil.which('foldseek'))
    if sha(exe)!=c['executable_sha256']:raise ValueError('Foldseek executable differs')
    database_pins={str(p):sha(p) for p in sorted(db.parent.glob('input*')) if p.is_file()}
    fields=['query','target','evalue','qstart','qend','qlen','tstart','tend','tlen','alnlen','qcov','tcov','alntmscore','qtmscore','ttmscore','lddt','rmsd']
    pref=a.output/'pairs';align=a.output/'alignments';table=a.output/'alignments.tsv'
    commands=[[str(exe),'tsv2db',str(pairs),str(pref),'--output-dbtype','7'],[str(exe),'structurealign',str(db),str(db),str(pref),str(align),'--threads','8','--alignment-type','2','--comp-bias-corr','0','--sort-by-structure-bits','0','-a','1','--alignment-mode','3','-e','1000000','-c','0','--tmscore-threshold','0','--exact-tmscore','0'],[str(exe),'convertalis',str(db),str(db),str(align),str(table),'--threads','8','--format-output',','.join(fields),'--exact-tmscore','0']]
    cfg={'cluster_receipt_sha256':sha(a.clusters/'receipt.json'),'database_file_sha256':database_pins,'script_sha256':sha(Path(__file__)),'executable_sha256':sha(exe),'commands':commands,'requested_pair_table_sha256':sha(pairs),'directed_pairs':len(expected),'criterion':'Reported evalue<=1e-5, qcov>=0.8, tcov>=0.8, alignment-normalized approximate TM score>=0.5. Reported values have finite output precision; near-boundary results require review.','interpretation':'Same Foldseek aligner, prefilter bypassed, scoring settings matched to clustering log. Both directions preserved. Independent of candidate-search screening, not an independent geometry implementation or all-pairs within-cluster proof.'}
    cp=a.output/'config.json';cp.write_text(json.dumps(cfg,indent=2)+'\n');print('Realigning',len(expected),'directed edges',flush=True)
    with (a.output/'foldseek.log').open('w') as log:
        for cmd in commands:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
    found={}
    with table.open() as f:
        for values in csv.reader(f,delimiter='\t'):
            row=dict(zip(fields,values));key=(row['query'],row['target'])
            if len(values)!=len(fields) or key not in expected or key in found:raise ValueError('Unexpected/duplicate alignment')
            numbers={k:float(row[k]) for k in fields[2:]}
            if not all(math.isfinite(v) for v in numbers.values()):raise ValueError('Nonfinite result')
            if int(row['qlen'])!=int(models[key[0]]['length']) or int(row['tlen'])!=int(models[key[1]]['length']):raise ValueError('Alignment length identity differs')
            row['passes_reported_cluster_criteria']=numbers['evalue']<=1e-5 and numbers['qcov']>=.8 and numbers['tcov']>=.8 and numbers['alntmscore']>=.5
            found[key]=row
    statuses=[]
    for rep,member in sorted(edges):
        forward=found.get((rep,member));reverse=found.get((member,rep));passed=[x is not None and x['passes_reported_cluster_criteria'] for x in [forward,reverse]]
        status='both_directions_pass' if all(passed) else 'one_direction_passes' if any(passed) else 'neither_direction_passes'
        statuses.append({'representative_input_id':rep,'member_input_id':member,'forward_alignment_returned':forward is not None,'reverse_alignment_returned':reverse is not None,'forward_passes':passed[0],'reverse_passes':passed[1],'status':status})
    write_table(a.output/'edge_review.tsv',statuses)
    if found:write_table(a.output/'directed_alignment_review.tsv',list(found.values()))
    result={'status':'complete_representative_member_edge_realignment','config_sha256':sha(cp),'representative_member_pairs':len(edges),'requested_directed_pairs':len(expected),'returned_directed_pairs':len(found),'missing_directed_pairs':len(expected)-len(found),'edge_status_counts':dict(Counter(r['status'] for r in statuses)),'interpretation':cfg['interpretation']+' Missing alignments are retained as failures of this validation; they do not prove biological dissimilarity. Original clusters are immutable and require this review before interpretation.','artifacts':{p.name:sha(p) for p in [pairs,table,a.output/'edge_review.tsv',a.output/'directed_alignment_review.tsv'] if p.exists()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
