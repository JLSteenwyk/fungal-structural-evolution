#!/usr/bin/env python3
"""Create exploratory whole-chain similarity clusters while retaining model provenance."""
import argparse,csv,fcntl,json,shutil,subprocess,time
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT,sha
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a fresh immutable cluster output; inspect partial failures before recovery')
    sources={'AlphaFold':ROOT/'results/structural_markers/gdm-expanded-v1','ESMFold':ROOT/'results/structural_markers/esmfold-partial-v1'};rows=[];pins={}
    for source,base in sources.items():
        checked_receipt(base);pins[source]=sha(base/'receipt.json')
        for m in sorted(json.loads((base/'model_provenance.json').read_text()),key=lambda r:r['model_id']):
            if sha(ROOT/m['path'])!=m['sha256']:raise ValueError('Changed source structure')
            rows.append({'cluster_input_id':'M'+str(len(rows)+1).zfill(6),'source':source,'model_id':m['model_id'],'model_version':m['version'],'sequence_sha256':m['sequence_sha256'],'length':m['length'],'mean_ca_plddt':m['mean_ca_plddt'],'model_path':m['path'],'model_sha256':m['sha256']})
    if len({(r['source'],r['model_id']) for r in rows})!=len(rows):raise ValueError('Duplicate source model')
    a.output.mkdir(parents=True);inputs=a.output/'inputs';inputs.mkdir();temp=a.output/'tmp';prefix=a.output/'clusters'
    for r in rows:(inputs/(r['cluster_input_id']+'.cif')).symlink_to((ROOT/r['model_path']).resolve())
    write_table(a.output/'input_models.tsv',rows)
    exe=Path(shutil.which('foldseek'));command=[str(exe),'easy-cluster',str(inputs),str(prefix),str(temp),'--threads','8','--split-memory-limit','32G','-s','7.5','--max-seqs','2000','-e','1e-5','-c','0.8','--cov-mode','0','--alignment-type','2','--tmscore-threshold','0.5','--tmscore-threshold-mode','0','--cluster-mode','2','--single-step-clustering','1','--cluster-reassign','1','--mask-bfactor-threshold','70','--remove-tmp-files','0']
    version=subprocess.run([str(exe),'version'],capture_output=True,text=True,check=True).stdout.strip()
    config={'source_receipt_sha256':pins,'input_table_sha256':sha(a.output/'input_models.tsv'),'script_sha256':sha(Path(__file__)),'executable_sha256':sha(exe),'version':version,'command':command,'models':len(rows),'residues':sum(r['length'] for r in rows),'interpretation':'Exploratory whole-chain 3Di+AA similarity clustering. Aligned coverage >=0.8 for both query and target, E<=1e-5, approximate alignment-normalized TM score threshold 0.5; greedy length-based single-step clustering with reassignment. pLDDT70 masks seeding only, not every scored residue. No orthology, homology, novelty or biological-function inference. Prefilter hit cap may limit recall; pairwise/member criteria and threshold sensitivity require later audits.'}
    cp=a.output/'config.json';cp.write_text(json.dumps(config,indent=2)+'\n');print('Starting',len(rows),'models',flush=True)
    with (a.output/'foldseek.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    assignments=a.output/'clusters_cluster.tsv';seen=set();mapping={r['cluster_input_id']:r for r in rows};clusters={};annotated=[]
    with assignments.open() as f:
        for line in f:
            representative,member=line.rstrip('\n').split('\t')
            if representative not in mapping or member not in mapping or member in seen:raise ValueError('Unexpected or duplicate cluster model ID')
            seen.add(member);clusters.setdefault(representative,[]).append(member);annotated.append({'representative_input_id':representative,**mapping[member]})
    if seen!=set(mapping) or any(rep not in members for rep,members in clusters.items()):raise ValueError('Incomplete cluster membership')
    write_table(a.output/'model_cluster_membership.tsv',annotated)
    result={'status':'complete_exploratory_frozen_model_clustering','config_sha256':sha(cp),'models':len(rows),'clusters':len(clusters),'singletons':sum(len(v)==1 for v in clusters.values()),'mixed_prediction_source_clusters':sum(len({mapping[m]['source'] for m in members})>1 for members in clusters.values()),'interpretation':config['interpretation'],'artifacts':{p.name:sha(p) for p in [a.output/'input_models.tsv',assignments,a.output/'model_cluster_membership.tsv']}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
