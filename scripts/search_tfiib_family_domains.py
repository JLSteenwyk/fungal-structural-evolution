#!/usr/bin/env python3
"""Search selected Pfam domains across all additional representative proteins."""
import argparse,json,re,shutil,subprocess,time
from pathlib import Path
from audit_busco_gene_copies import ROOT,sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable focused family search')
    pfam_path=ROOT/'metadata/pfam_release_receipt.json';pfam=json.loads(pfam_path.read_text())
    hmm_info=next(x for x in pfam['files'] if x['file']=='Pfam-A.hmm.gz');hmm=ROOT/'data/pfam/38.2'/hmm_info['uncompressed_file']
    if sha(hmm)!=hmm_info['uncompressed_sha256']:raise ValueError('Changed Pfam HMM library')
    inputs=ROOT/'data/domains/full-inputs-v1';inp=json.loads((inputs/'receipt.json').read_text())
    if inp['status']!='complete_search_input_preparation' or inp['taxa']!=526:raise ValueError('Full representative inputs required')
    for name,h in inp['artifacts'].items():
        if sha(inputs/name)!=h:raise ValueError('Changed full protein inputs')
    marker=ROOT/'results/domains/marker-annotations-v1';mp=marker/'receipt.json';mr=json.loads(mp.read_text())
    for name,h in mr['artifacts'].items():
        if sha(marker/name)!=h:raise ValueError('Changed reusable marker hits')
    ms=ROOT/'results/domains/marker-search-v1';mc=json.loads((ms/'config.json').read_text())
    if mc['pfam_receipt_sha256']!=sha(pfam_path) or mc['input_receipt_sha256']!=inp['marker_input_receipt_sha256']:raise ValueError('Marker/full search provenance differs')
    if mr['search_receipt_sha256']!=sha(ms/'receipt.json'):raise ValueError('Changed marker search receipt')
    search=shutil.which('hmmsearch');fetch=shutil.which('hmmfetch')
    if sha(Path(search))!=mc['executable_sha256']:raise ValueError('Different HMMER search binary')
    version=subprocess.run([search,'-h'],capture_output=True,text=True,check=True).stdout.splitlines()[1]
    a.output.mkdir(parents=True);ids=a.output/'profile_accessions.txt';requested=['PF00382.25','PF07741.19','PF08271.18'];ids.write_text('\n'.join(requested)+'\n')
    profiles=a.output/'selected_profiles.hmm';fetch_command=[fetch,'-f',str(hmm),str(ids)]
    with profiles.open('w') as f:subprocess.run(fetch_command,stdout=f,check=True)
    text=profiles.read_text();accessions=re.findall(r'^ACC\s+(\S+)',text,re.M);names=re.findall(r'^NAME\s+(\S+)',text,re.M)
    if set(accessions)!=set(requested) or len(accessions)!=3 or text.count('\n//')!=3:raise ValueError('Incomplete profile extraction')
    table=a.output/'hits.domtblout';report=a.output/'search.txt'
    command=[search,'--cut_ga','--cpu','2','--seed','42','--noali','--domtblout',str(table),'-o',str(report),str(profiles),str(inputs/'additional_sequences.faa')]
    config={'profiles':requested,'profile_names':names,'pfam_receipt_sha256':sha(pfam_path),'full_input_receipt_sha256':sha(inputs/'receipt.json'),'marker_annotation_receipt_sha256':sha(mp),'selected_profile_sha256':sha(profiles),'fetch_command':fetch_command,'fetch_executable_sha256':sha(Path(fetch)),'search_command':command,'search_executable_sha256':sha(Path(search)),'version':version,'additional_unique_sequences':inp['additional_unique_sequences'],'representative_proteins':inp['representative_proteins'],'taxa':inp['taxa'],'interpretation':'Focused acceleration of three domain searches across all additional representative sequences; existing marker hits are available for later combination. Gathering thresholds match the primary Pfam searches, but E-values use this additional database. Full-family membership/phylogeny/reconciliation still pending.'}
    cp=a.output/'config.json';cp.write_text(json.dumps(config,indent=2)+'\n');print('Starting',inp['additional_unique_sequences'],'additional sequences',requested,flush=True)
    start=time.monotonic()
    with (a.output/'stderr.log').open('w') as err:subprocess.run(command,stdout=err,stderr=err,check=True)
    if '[ok]' not in report.read_text()[-100:]:raise ValueError('Search report lacks completion marker')
    queries=re.findall(r'^Query:\s+(\S+)',report.read_text(),re.M)
    if queries!=names:raise ValueError('Search did not finish the exact profile set')
    hits=0
    with table.open() as f:
        for line in f:
            if not line.strip() or line.startswith('#'):continue
            fields=line.split(maxsplit=22)
            if len(fields)<22 or fields[4] not in requested or not fields[0].startswith('S'):raise ValueError('Unexpected hit identity/orientation')
            hits+=1
    result={'status':'complete_focused_tfiib_domain_search','config_sha256':sha(cp),'profiles':requested,'additional_unique_sequences':inp['additional_unique_sequences'],'domain_hit_rows':hits,'elapsed_seconds':time.monotonic()-start,'script_sha256':sha(Path(__file__)),'interpretation':config['interpretation'],'artifacts':{p.name:sha(p) for p in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
