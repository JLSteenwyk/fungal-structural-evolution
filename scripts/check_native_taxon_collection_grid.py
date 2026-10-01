#!/usr/bin/env python3
"""Full 16-job collection I/O fixture; all inference/provenance is synthetic."""
import argparse
import csv
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import sys
from check_native_taxon_pmsf_reader_cases import make_case, write, rehash
from run_ortholog_pair_guide_comparison import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    root=args.output.resolve();root.mkdir(parents=True,exist_ok=False)
    inputs=root/'inputs';inputs.mkdir();native=root/'native';native.mkdir();matrices=[];jobs=[]
    policies=['sparse','occupancy','hybrids','labels'];counts=[5,6,8,11]
    for policy,n in zip(policies,counts):
        for alignment in ['profile','mafft']:
            base=root/'fixtures'/(policy+'-'+alignment);make_case(base,n)
            prepared=inputs/(policy+'-'+alignment);prepared.mkdir();shutil.copyfile(base/'matrix.faa',prepared/'matrix.faa')
            with (prepared/'taxa.tsv').open('w') as handle:
                writer=csv.DictWriter(handle,['taxon_id','study_role'],delimiter='\t',lineterminator='\n');writer.writeheader()
                writer.writerows(dict(taxon_id='T_'+str(i),study_role='outgroup' if i>=n-2 else 'ingroup') for i in range(n))
            matrices.append(dict(alignment=alignment,policy=policy,path=str(prepared),taxa=n,columns=17,
                                 roles=dict(ingroup=n-2,outgroup=2),matrix_sha256=sha(prepared/'matrix.faa')))
            guide=native/'guides'/(alignment+'-'+policy);guide.mkdir(parents=True)
            shutil.copyfile(base/'run'/'input_guide.treefile',guide/'guide.treefile')
            write(guide/'receipt.json',dict(status='checked_full_taxon_sensitivity_conditioning_guide',taxa=n,
                matrix_sha256=sha(prepared/'matrix.faa'),artifacts={'guide.treefile':sha(guide/'guide.treefile')}))
        for alignment,gl in itertools.product(['profile','mafft'],repeat=2):
            label=policy+'-'+alignment+'-'+gl;folder=native/label;shutil.copytree(root/'fixtures'/(policy+'-'+alignment)/'run',folder)
            jobs.append(dict(label=label,alignment=alignment,guide_alignment=gl,policy=policy))
    write(inputs/'receipt.json',dict(matrices=matrices))
    closed=root/'closed_inputs.json';write(closed,dict(status='complete_verified_native_species_taxon_refit_inputs',
        services=[{'synthetic_not_a_journal':1},{'synthetic_not_a_journal':2}],summary=dict(matrices=8),source_hashes={}))
    np=root/'native_plan.json';audit_impl='scripts/audit_species_taxon_pmsf.py'
    write(np,dict(output=str(native),inputs=str(inputs),input_completion=str(closed),jobs=jobs,pins={audit_impl:sha(audit_impl)}))
    records=[]
    for job in jobs:
        folder=native/job['label'];spec=next(m for m in matrices if m['alignment']==job['alignment'] and m['policy']==job['policy'])
        c=json.loads((folder/'config.json').read_text());c.update(batch_plan_sha256=sha(np),**{k:job[k] for k in ['policy','alignment','guide_alignment']})
        c['command'][c['command'].index('-s')+1]=str(Path(spec['path'])/'matrix.faa')
        c['command'][c['command'].index('--tree-freq')+1]=str(folder/'input_guide.treefile')
        c['pinned_files']={c['command'][c['command'].index('-s')+1]:spec['matrix_sha256']};write(folder/'config.json',c)
        r=json.loads((folder/'receipt.json').read_text());r['config_sha256']=sha(folder/'config.json')
        r['artifacts']={p.name:sha(p) for p in folder.iterdir() if p.is_file() and p.name!='receipt.json'};write(folder/'receipt.json',r)
        a=json.loads((folder/'audit'/'receipt.json').read_text());a.update(status='passed_taxon_sensitivity_pmsf_profile_tree_and_bootstrap_readback',
            source_receipt_sha256=sha(folder/'receipt.json'),script_sha256=sha(audit_impl));write(folder/'audit'/'receipt.json',a)
        records.append(dict(**job,run=str(folder),audit=str(folder/'audit'),taxa=spec['taxa'],columns=17,
            run_receipt_sha256=sha(folder/'receipt.json'),audit_receipt_sha256=sha(folder/'audit'/'receipt.json')))
    for m in matrices:
        gp=native/'guides'/(m['alignment']+'-'+m['policy'])/'receipt.json';gr=json.loads(gp.read_text());gr['plan_sha256']=sha(np);write(gp,gr)
    write(native/'receipt.json',dict(status='complete_native_taxon_pmsf_batch_pending_full_independent_collection_readback',
        plan_sha256=sha(np),runs=records,source_hashes={}))
    plan=root/'collection_plan.json';write(plan,dict(native_plan=str(np),output=str(root/'collection'),pins={},
        scope='Entire16-job collection software I/O fixture. Every tree/profile/likelihood/inference status and prior-journal contract is synthetic; no biological run or completion evidence.'))
    result=subprocess.run([sys.executable,'scripts/readback_native_taxon_pmsf_collection.py','--plan',str(plan)],capture_output=True,text=True)
    (root/'reader.stdout.log').write_text(result.stdout+result.stderr)
    assert result.returncode==0,result.stderr
    proof=json.loads((root/'collection'/'receipt.json').read_text())
    assert proof['run_count']==16 and proof['tree_views']==proof['role_boundary_rows']==proof['boundary_views_with_role_split']==32
    assert proof['raw_bootstrap_trees']==16000 and proof['site_profiles']==272
    write(root/'receipt.json',dict(status='passed_complete_native_taxon_collection_grid_software_fixture',run_count=16,tree_views=32,
        bootstrap_trees=16000,site_profiles=272,collection_receipt_sha256=sha(root/'collection'/'receipt.json'),
        reader_script_sha256=sha('scripts/readback_native_taxon_pmsf_collection.py'),fixture_script_sha256=sha(Path(__file__)),
        scope='Complete16 synthetic job identities/matrix-role grids/conditioning-guide links/raw outputs/32 boundary rows tested. Prior statuses and journal records are synthetic contracts, not completed original services or production data. Full-sized production readback and journal closure remain required.'))
    print(json.dumps(dict(status='passed_complete_native_taxon_collection_grid_software_fixture',run_count=16,tree_views=32)),flush=True)


if __name__=='__main__':main()
