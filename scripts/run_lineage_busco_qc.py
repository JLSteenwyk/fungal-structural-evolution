#!/usr/bin/env python3
"""Run explicit lineage-specific protein QC across the full fungal ingroup."""
import argparse
import csv
import fcntl
import json
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from prepare_lineage_busco_datasets import ROOT, sha

SPECIFIC={x:x.lower()+'_odb12.2' for x in ['Ascomycota','Basidiomycota','Chytridiomycota','Microsporidia','Mucoromycota']}


def validate_summary(path,dataset):
    r=json.loads(path.read_text());v=r['results'];lineage=r['lineage_dataset']
    if lineage['name']!=dataset['dataset'] or lineage['creation_date']!=dataset['config']['creation_date']:raise ValueError('BUSCO dataset identity differs')
    n=int(dataset['config']['number_of_BUSCOs'])
    if int(v['n_markers'])!=n or int(lineage['number_of_buscos'])!=n:raise ValueError('Marker denominator differs')
    fields=['Single copy BUSCOs','Multi copy BUSCOs','Fragmented BUSCOs','Missing BUSCOs']
    if any(int(v[k])<0 for k in fields) or sum(int(v[k]) for k in fields)!=n or int(v['Complete BUSCOs'])!=int(v[fields[0]])+int(v[fields[1]]):raise ValueError('Inconsistent BUSCO category counts')
    return {k:v[k] for k in ['n_markers','Complete BUSCOs',*fields]}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--datasets',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    collection=json.loads((a.datasets/'receipt.json').read_text())
    if collection['status']!='complete_pinned_lineage_busco_dataset_collection':raise ValueError('Completed pinned dataset collection required')
    datasets={r['dataset']:r for r in collection['datasets']}
    for dataset in datasets.values():
        for name,h in dataset['files_sha256'].items():
            if sha(ROOT/dataset['path']/name)!=h:raise ValueError('Changed lineage dataset')
    manifest=ROOT/'metadata/analysis_manifest.tsv';inputs=ROOT/'metadata/qc_input_receipts.json'
    with manifest.open() as f:taxa=list(csv.DictReader(f,delimiter='\t'))
    sequences={r['taxon_id']:r for r in json.loads(inputs.read_text())}
    if len(taxa)!=len({r['taxon_id'] for r in taxa}) or set(sequences)!={r['taxon_id'] for r in taxa}:raise ValueError('QC and manifest identities differ')
    assignments=[]
    for row in taxa:
        group=row['lineage'].split(';')[0];ingroup=row['study_role']=='ingroup'
        name=SPECIFIC.get(group,'fungi_odb12.2') if ingroup else ''
        assignments.append({'taxon_id':row['taxon_id'],'species_name':row['species_name'],'study_role':row['study_role'],'lineage_group':group,'dataset':name,'selection_reason':('Exact manifest-group panel' if group in SPECIFIC else 'Fungal parent panel; no selected narrower panel') if ingroup else 'Retain existing eukaryotic QC; fungal panel not applied'})
    a.output.mkdir(parents=True,exist_ok=True);lock=(a.output/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    exe=ROOT/'.cache/envs/busco/bin/busco';version=subprocess.run(['conda','run','--prefix',str(exe.parent.parent),'busco','--version'],capture_output=True,text=True,check=True).stdout.strip()
    config={'dataset_receipt_sha256':sha(a.datasets/'receipt.json'),'manifest_sha256':sha(manifest),'input_receipts_sha256':sha(inputs),'script_sha256':sha(Path(__file__)),'busco_executable_sha256':sha(exe),'version':version,'workers':4,'threads_per_job':4,'mode':'proteins','offline':True,'selection':SPECIFIC,'fallback':'fungi_odb12.2','taxa':len(taxa),'ingroup_jobs':sum(bool(r['dataset']) for r in assignments),'interpretation':'Complementary lineage-panel QC. Broad scores retained, outgroups not tested against fungal panels; no automatic filtering or contamination/ploidy inference.'}
    cp=a.output/'config.json'
    if cp.exists() and json.loads(cp.read_text())!=config:raise ValueError('Changed lineage QC configuration')
    cp.write_text(json.dumps(config,indent=2)+'\n')
    with (a.output/'taxon_dataset_assignments.tsv').open('w') as f:
        writer=csv.DictWriter(f,list(assignments[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(assignments)
    (a.output/'logs').mkdir(exist_ok=True);(a.output/'jobs').mkdir(exist_ok=True)
    def run(assignment):
        taxon=assignment['taxon_id'];row=sequences[taxon];ds=datasets[assignment['dataset']]
        if sha(ROOT/row['input_path'])!=row['sha256']:raise ValueError('Changed taxon proteome')
        rp=a.output/'jobs'/(taxon+'.json');folder=a.output/taxon
        if rp.exists():
            old=json.loads(rp.read_text())
            if old['status']=='complete_lineage_protein_qc' and old['config_sha256']==sha(cp) and old['input_sha256']==row['sha256']:
                summary=ROOT/old['summary_path']
                if sha(summary)!=old['summary_sha256']:raise ValueError('Changed completed taxon summary')
                validate_summary(summary,ds);return old
            raise ValueError('Existing failed taxon needs explicit recovery')
        if folder.exists():raise FileExistsError('Partial taxon output needs explicit recovery')
        cmd=['conda','run','--prefix',str(exe.parent.parent),'busco','-i',str(ROOT/row['input_path']),'-m','proteins','-l',str(ROOT/ds['path']),'--offline','-c','4','-o',taxon,'--out_path',str(a.output.resolve())]
        start=time.monotonic()
        with (a.output/'logs'/(taxon+'.log')).open('w') as f:done=subprocess.run(cmd,stdout=f,stderr=f,cwd=ROOT)
        result={'taxon_id':taxon,'dataset':ds['dataset'],'input_sha256':row['sha256'],'config_sha256':sha(cp),'command':cmd,'returncode':done.returncode,'elapsed_seconds':time.monotonic()-start,'status':'failed_requires_review'}
        summaries=list(folder.glob('short_summary.specific.*.json'))
        if done.returncode==0 and len(summaries)==1:
            counts=validate_summary(summaries[0],ds)
            result.update(status='complete_lineage_protein_qc',counts=counts,summary_path=str(summaries[0].resolve().relative_to(ROOT)),summary_sha256=sha(summaries[0]))
        rp.write_text(json.dumps(result,indent=2)+'\n');return result
    results=[]
    with ThreadPoolExecutor(max_workers=4) as pool:
        for future in as_completed([pool.submit(run,r) for r in assignments if r['dataset']]):
            r=future.result();results.append(r);print(len(results),r['taxon_id'],r['dataset'],r['status'],flush=True)
    result={'status':'complete_full_ingroup_lineage_qc' if all(r['status']=='complete_lineage_protein_qc' for r in results) else 'incomplete_lineage_qc_requires_review','config_sha256':sha(cp),'jobs':len(results),'completed':sum(r['status']=='complete_lineage_protein_qc' for r in results),'outgroups_retaining_broad_qc':sum(not r['dataset'] for r in assignments),'artifacts':{'taxon_dataset_assignments.tsv':sha(a.output/'taxon_dataset_assignments.tsv'),**{str(p.relative_to(a.output)):sha(p) for p in (a.output/'jobs').glob('*.json')}}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='artifacts'},indent=2))


if __name__=='__main__':main()
